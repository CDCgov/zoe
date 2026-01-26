//! Implementations and methods for iterating over CIGAR strings.

use crate::{
    alignment::{AlignmentStates, StatesSequence},
    data::cigar::{Cigar, CigarError, CigarView, CigarViewMut, Ciglet, USIZE_WIDTH, is_valid_op},
    unwrap_or_return_some_err,
};

/// An iterator over a CIGAR string.
///
/// This iterator has the following unchecked behavior which is worth noting:
///
/// - If any increment is missing, a [`Ciglet`] with an increment of 0 is
///   yielded.
/// - If an invalid operation is encountered, an increment is bigger than
///   [`usize::MAX`], or an increment has no following operation, then the
///   iterator ends early.
/// - If adjacent operations are the same, this iterator will yield them
///   separately as two [`Ciglet`] values.
#[derive(Clone, Debug)]
pub struct CigletIterator<'a> {
    buffer: &'a [u8],
    valid:  bool,
}

impl<'a> CigletIterator<'a> {
    /// Constructs an iterator from a CIGAR byte buffer.
    #[inline]
    pub(crate) fn new(buffer: &'a [u8]) -> Self {
        CigletIterator { buffer, valid: true }
    }

    /// For internal `eq()` only
    #[inline]
    pub(crate) fn valid(&self) -> bool {
        self.valid
    }
}

impl Iterator for CigletIterator<'_> {
    type Item = Ciglet;

    fn next(&mut self) -> Option<Self::Item> {
        let mut num = 0;
        let mut index = 0;
        let limit = std::cmp::min(self.buffer.len(), USIZE_WIDTH - 1);
        while index != limit {
            let b = self.buffer[index];

            if b.is_ascii_digit() {
                num *= 10;
                num += usize::from(b - b'0');
            } else if is_valid_op(b) && num > 0 {
                self.buffer = &self.buffer[index + 1..];
                return Some(Ciglet { inc: num, op: b });
            } else {
                self.valid = false;
                return None;
            }
            index += 1;
        }

        while index != self.buffer.len() {
            let b = self.buffer[index];

            if b.is_ascii_digit() {
                num = num.checked_mul(10)?;
                num = num.checked_add(usize::from(b - b'0'))?;
            } else if is_valid_op(b) && num > 0 {
                self.buffer = &self.buffer[index + 1..];
                return Some(Ciglet { inc: num, op: b });
            } else {
                self.valid = false;
                return None;
            }
            index += 1;
        }

        None
    }
}

impl DoubleEndedIterator for CigletIterator<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let op = match self.buffer {
            [rest @ .., op] if is_valid_op(*op) => {
                self.buffer = rest;
                *op
            }
            _ => return None,
        };

        let mut num = 0;

        if let [rest @ .., b] = self.buffer {
            if b.is_ascii_digit() {
                num += usize::from(b - b'0');
            } else {
                return (num > 0).then_some(Ciglet { inc: num, op });
            }
            self.buffer = rest;
        } else {
            return None;
        }

        let mut multiplier = 1;
        let mut index = 1;

        while let [rest @ .., b] = self.buffer
            && index < USIZE_WIDTH - 1
        {
            if b.is_ascii_digit() {
                multiplier *= 10;
                num += usize::from(b - b'0') * multiplier;
            } else {
                return (num > 0).then_some(Ciglet { inc: num, op });
            }
            self.buffer = rest;
            index += 1;
        }
        if self.buffer.is_empty() {
            return (num > 0).then_some(Ciglet { inc: num, op });
        }

        let mut num_zeros = 0;

        while let [rest @ .., b] = self.buffer {
            if b.is_ascii_digit() {
                if *b == b'0' {
                    num_zeros += 1;
                } else {
                    for _ in 0..=num_zeros {
                        multiplier = multiplier.checked_mul(10)?;
                    }
                    num = usize::from(b - b'0').checked_mul(multiplier)?.checked_add(num)?;
                }
                self.buffer = rest;
            } else {
                return (num > 0).then_some(Ciglet { inc: num, op });
            }
        }

        Some(Ciglet { inc: num, op })
    }
}

impl StatesSequence for CigletIterator<'_> {
    #[inline]
    fn peek_op(&self) -> Option<u8> {
        self.buffer.iter().find(|x| !x.is_ascii_digit()).copied()
    }

    #[inline]
    fn peek_back_op(&self) -> Option<u8> {
        self.buffer.last().copied()
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }

    #[inline]
    fn next_ciglet(&mut self) -> Option<Ciglet> {
        self.next()
    }

    #[inline]
    fn next_ciglet_back(&mut self) -> Option<Ciglet> {
        self.next_back()
    }
}

impl<'a> IntoIterator for &'a Cigar {
    type Item = Ciglet;
    type IntoIter = CigletIterator<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        CigletIterator::new(&self.0)
    }
}

impl<'a> IntoIterator for CigarView<'a> {
    type Item = Ciglet;
    type IntoIter = CigletIterator<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        CigletIterator::new(self.as_bytes())
    }
}

impl<'a> IntoIterator for &CigarView<'a> {
    type Item = Ciglet;
    type IntoIter = CigletIterator<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        CigletIterator::new(self.as_bytes())
    }
}

impl<'a> IntoIterator for &'a CigarViewMut<'a> {
    type Item = Ciglet;
    type IntoIter = CigletIterator<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        CigletIterator::new(self.as_bytes())
    }
}

/// A trait for types from which an iterator of [`Ciglet`] values can be made.
///
/// This is useful for trait bounds, where a generic implementation should
/// function over anything that "looks like" a CIGAR string. In particular, this
/// is stronger than `IntoIterator<Item = Ciglet>` since it is implemented on
/// both owned types and references (calling [`copied`] on the iterator for the
/// latter).
///
/// [`copied`]: Iterator::copied
pub trait ToCigletIterator {
    /// The type of the iterator, which must yield [`Ciglet`] items.
    type Iter<'a>: Iterator<Item = Ciglet>
    where
        Self: 'a;

    /// Creates an iterator over the [`Ciglet`] values.
    fn to_ciglet_iterator(&self) -> Self::Iter<'_>;
}

impl<'a> ToCigletIterator for CigletIterator<'a> {
    type Iter<'b>
        = CigletIterator<'a>
    where
        Self: 'b;

    #[inline]
    fn to_ciglet_iterator(&self) -> Self::Iter<'_> {
        self.clone()
    }
}

impl ToCigletIterator for Cigar {
    type Iter<'a> = CigletIterator<'a>;

    #[inline]
    fn to_ciglet_iterator(&self) -> CigletIterator<'_> {
        self.iter()
    }
}

impl ToCigletIterator for &Cigar {
    type Iter<'a>
        = CigletIterator<'a>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> CigletIterator<'_> {
        self.iter()
    }
}

impl ToCigletIterator for &mut Cigar {
    type Iter<'a>
        = CigletIterator<'a>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> CigletIterator<'_> {
        self.iter()
    }
}

impl<'b> ToCigletIterator for CigarView<'b> {
    type Iter<'a>
        = CigletIterator<'b>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> CigletIterator<'b> {
        self.iter()
    }
}

impl<'b> ToCigletIterator for &CigarView<'b> {
    type Iter<'a>
        = CigletIterator<'b>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> CigletIterator<'b> {
        self.iter()
    }
}

impl ToCigletIterator for &CigarViewMut<'_> {
    type Iter<'a>
        = CigletIterator<'a>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> CigletIterator<'_> {
        self.iter()
    }
}

impl ToCigletIterator for AlignmentStates {
    type Iter<'a> = std::iter::Copied<std::slice::Iter<'a, Ciglet>>;

    #[inline]
    fn to_ciglet_iterator(&self) -> Self::Iter<'_> {
        self.iter().copied()
    }
}

impl ToCigletIterator for &AlignmentStates {
    type Iter<'a>
        = std::iter::Copied<std::slice::Iter<'a, Ciglet>>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> Self::Iter<'_> {
        self.iter().copied()
    }
}

impl ToCigletIterator for &mut AlignmentStates {
    type Iter<'a>
        = std::iter::Copied<std::slice::Iter<'a, Ciglet>>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> Self::Iter<'_> {
        self.iter().copied()
    }
}

impl ToCigletIterator for Vec<Ciglet> {
    type Iter<'a> = std::iter::Copied<std::slice::Iter<'a, Ciglet>>;

    #[inline]
    fn to_ciglet_iterator(&self) -> Self::Iter<'_> {
        self.iter().copied()
    }
}

impl ToCigletIterator for &[Ciglet] {
    type Iter<'a>
        = std::iter::Copied<std::slice::Iter<'a, Ciglet>>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> Self::Iter<'_> {
        self.iter().copied()
    }
}

impl ToCigletIterator for &mut [Ciglet] {
    type Iter<'a>
        = std::iter::Copied<std::slice::Iter<'a, Ciglet>>
    where
        Self: 'a;

    #[inline]
    fn to_ciglet_iterator(&self) -> Self::Iter<'_> {
        self.iter().copied()
    }
}

/// Similar to [`CigletIterator`], but performs checking to ensure each ciglet
/// is valid during iteration.
///
/// Exhausting the iterator and checking for errors is equivalent to
/// [`Cigar::is_valid`], except that `*` will not be considered valid by this
/// method.
///
/// The following assumptions are validated by this iterator:
///
/// - The CIGAR operations must be among `M, I, D, N, S, H, P, X, =`
/// - Every operation in the CIGAR string must have a preceding increment
/// - Every increment must be followed by an operation
/// - The increment for each operation must be non-zero and less than or equal
///   to [`usize::MAX`]
///
/// If any of these are not valid, a [`CigarError`] is returned. Note that this
/// iterator does not confirm that adjacent operations are distinct.
pub(crate) struct CigletIteratorChecked<'a> {
    /// The bytes remaining to decode
    buffer: &'a [u8],
}

impl<'a> CigletIteratorChecked<'a> {
    /// Creates a new [`CigletIteratorChecked`] from the provided bytes.
    #[inline]
    pub(crate) fn new(bytes: &'a [u8]) -> CigletIteratorChecked<'a> {
        Self { buffer: bytes }
    }
}

impl Iterator for CigletIteratorChecked<'_> {
    type Item = Result<Ciglet, CigarError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut num = 0;
        let mut has_number = false;
        let mut index = 0;
        let limit = std::cmp::min(self.buffer.len(), USIZE_WIDTH - 1);
        while index != limit {
            let b = self.buffer[index];

            if b.is_ascii_digit() {
                num *= 10;
                num += usize::from(b - b'0');
                has_number = true;
            } else if is_valid_op(b) {
                if !has_number {
                    return Some(Err(CigarError::MissingInc));
                } else if num == 0 {
                    return Some(Err(CigarError::IncZero));
                }
                self.buffer = &self.buffer[index + 1..];
                return Some(Ok(Ciglet { inc: num, op: b }));
            } else {
                return Some(Err(CigarError::InvalidOperation));
            }
            index += 1;
        }

        while index != self.buffer.len() {
            let b = self.buffer[index];

            if b.is_ascii_digit() {
                num = unwrap_or_return_some_err!(
                    num.checked_mul(10)
                        .and_then(|n| n.checked_add(usize::from(b - b'0')))
                        .ok_or(CigarError::IncOverflow)
                );
                has_number = true;
            } else if is_valid_op(b) {
                if !has_number {
                    return Some(Err(CigarError::MissingInc));
                } else if num == 0 {
                    return Some(Err(CigarError::IncZero));
                }
                self.buffer = &self.buffer[index + 1..];
                return Some(Ok(Ciglet { inc: num, op: b }));
            } else {
                return Some(Err(CigarError::InvalidOperation));
            }
            index += 1;
        }

        if has_number {
            return Some(Err(CigarError::MissingOp));
        }

        None
    }
}
