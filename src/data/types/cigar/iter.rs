//! Implementations and methods for iterating over CIGAR strings.

use crate::{
    alignment::{AlignmentStates, StatesSequence},
    data::cigar::{Cigar, CigarError, CigarView, CigarViewMut, Ciglet, USIZE_WIDTH, is_valid_op},
    unwrap_or_return_some_err,
};
use std::hint::cold_path;

/// An iterator over a CIGAR string.
///
/// This iterator has the following unchecked behavior which is worth noting:
///
/// - If any increment is zero, a [`Ciglet`] with an increment of 0 is yielded.
///   If the increment field is empty, the iterator ends early.
/// - If an invalid operation is encountered, an increment is bigger than
///   [`usize::MAX`], or an increment has no following operation, then the
///   iterator ends early.
/// - If adjacent operations are the same, this iterator will yield them
///   separately as two [`Ciglet`] values.
#[derive(Clone, Debug)]
pub struct CigletIterator<'a> {
    buffer: &'a [u8],
    // TODO: Shouldn't this handle overflow?
    /// Whether or not an invalid state has been reached via [`next`] (namely,
    /// an empty increment field or an invalid operation). This is not updated
    /// with [`next_back`]. This is used for the [`Eq`] implementation of
    /// [`AlignmentStates`] with [`Cigar`].
    ///
    /// [`next`]: CigletIterator::next
    /// [`next_back`]: CigletIterator::next_back
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
            } else if is_valid_op(b) {
                if index == 0 {
                    // The increment field was completely missing, which is
                    // invalid
                    self.valid = false;
                    return None;
                }

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
            } else if is_valid_op(b) {
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
        let [buffer @ .., op] = self.buffer else {
            return None;
        };
        let mut buffer = buffer;
        let op = *op;

        if !is_valid_op(op) {
            return None;
        }

        let mut num = 0;

        if let [rest @ .., b] = buffer {
            if b.is_ascii_digit() {
                num += usize::from(b - b'0');
            } else {
                return None;
            }
            buffer = rest;
        } else {
            return None;
        }

        // At least one literal digit was present in the increment field, so
        // even if `num` is 0, it is still valid to yield a Ciglet

        let mut multiplier = 1;
        let mut index = 1;

        while let [rest @ .., b] = buffer
            && index < USIZE_WIDTH - 1
        {
            if b.is_ascii_digit() {
                multiplier *= 10;
                num += usize::from(b - b'0') * multiplier;
            } else {
                self.buffer = buffer;
                return Some(Ciglet { inc: num, op });
            }
            buffer = rest;
            index += 1;
        }

        if buffer.is_empty() {
            self.buffer = buffer;
            return Some(Ciglet { inc: num, op });
        }

        let mut num_zeros = 0;

        while let [rest @ .., b] = buffer {
            if b.is_ascii_digit() {
                if *b == b'0' {
                    num_zeros += 1;
                } else {
                    for _ in 0..=num_zeros {
                        multiplier = multiplier.checked_mul(10)?;
                    }
                    num = usize::from(b - b'0').checked_mul(multiplier)?.checked_add(num)?;
                }
                buffer = rest;
            } else {
                self.buffer = buffer;
                return Some(Ciglet { inc: num, op });
            }
        }

        self.buffer = buffer;
        Some(Ciglet { inc: num, op })
    }
}

impl StatesSequence for CigletIterator<'_> {
    #[inline]
    fn peek_op(&mut self) -> Option<u8> {
        loop {
            let op_pos = self.buffer.iter().position(|x| !x.is_ascii_digit())?;
            let (inc, rest) = self.buffer.split_at(op_pos);
            let op = rest[0];

            if inc.iter().all(|b| *b == b'0') {
                cold_path();
                self.buffer = &self.buffer[op_pos + 1..];
            } else {
                return Some(op);
            }
        }
    }

    #[inline]
    fn peek_back_op(&mut self) -> Option<u8> {
        loop {
            let (op, rest) = self.buffer.split_last()?;

            let inc_start_idx = rest
                .iter()
                .rposition(|x| !x.is_ascii_digit())
                .map_or(0, |next_op_idx| next_op_idx + 1);
            let inc = &rest[inc_start_idx..];

            if inc.iter().all(|b| *b == b'0') {
                cold_path();
                self.buffer = &self.buffer[..inc_start_idx];
            } else {
                return Some(*op);
            }
        }
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }

    #[inline]
    fn next_ciglet(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = self.next()?;

            if ciglet.inc == 0 {
                cold_path();
            } else {
                return Some(ciglet);
            }
        }
    }

    #[inline]
    fn next_ciglet_back(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = self.next_back()?;

            if ciglet.inc == 0 {
                cold_path();
            } else {
                return Some(ciglet);
            }
        }
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
