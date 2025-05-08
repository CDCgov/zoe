use crate::data::mappings::IS_CIGAR;

#[cfg(test)]
mod test;

#[cfg(test)]
mod benches;

mod error;

pub use error::*;

/// A [CIGAR string] of length-opcode pairs used in sequence alignment.
///
/// [CIGAR string]: https://en.wikipedia.org/wiki/Sequence_alignment#CIGAR_Format
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Cigar(pub(crate) Vec<u8>);

#[cfg(target_pointer_width = "16")]
const USIZE_WIDTH: usize = 5;
#[cfg(target_pointer_width = "32")]
const USIZE_WIDTH: usize = 10;
#[cfg(target_pointer_width = "64")]
const USIZE_WIDTH: usize = 20;

impl Cigar {
    /// Creates a CIGAR string from a Vec of bytes without checking for validity
    #[inline]
    #[must_use]
    pub fn from_vec_unchecked(v: Vec<u8>) -> Self {
        Cigar(v)
    }

    /// Creates a CIGAR string from a slice of bytes without checking for validity.
    #[inline]
    #[must_use]
    pub fn from_slice_unchecked<T: AsRef<[u8]>>(v: T) -> Self {
        Cigar(v.as_ref().to_vec())
    }

    /// Creates a new empty CIGAR string
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Cigar(Vec::new())
    }

    /// Expands a CIGAR to an opcode-only byte string, e.g. "3M" → "MMM".
    #[must_use]
    pub(crate) fn expand_cigar(&self) -> ExpandedCigar {
        let mut expanded = Vec::new();

        for Ciglet { inc, op } in self {
            expanded.extend(std::iter::repeat_n(op, inc));
        }

        ExpandedCigar(expanded)
    }

    /// Sums lengths for operations `M`, `D`, `N`, `X`, and `=`.
    #[must_use]
    pub fn match_length(&self) -> usize {
        self.into_iter()
            .filter(|Ciglet { inc: _, op }| matches!(op, b'M' | b'D' | b'N' | b'X' | b'='))
            .map(|Ciglet { inc, .. }| inc)
            .sum()
    }

    /// Returns an iterator of [`Ciglet`] for the CIGAR.
    #[inline]
    #[must_use]
    pub fn iter(&self) -> CigletIterator<'_> {
        CigletIterator::new(&self.0)
    }

    /// Returns an iterator over the opcodes in the CIGAR string, with each
    /// repeated the number of times specified in the CIGAR string.
    #[inline]
    pub fn expanded_cigar_iter(&self) -> impl Iterator<Item = u8> {
        self.iter().flat_map(|Ciglet { inc, op }| std::iter::repeat_n(op, inc))
    }

    /// Validates that a CIGAR only contains valid (inc, op) pairs
    #[inline]
    #[must_use]
    pub fn is_valid(&self) -> bool {
        Cigar::try_from(&self.0[..]).is_ok()
    }
}

impl Default for Cigar {
    #[inline]
    fn default() -> Self {
        Self::new()
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

impl std::fmt::Display for Cigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.0.is_empty() || self.0 == b"*" {
            write!(f, "*")
        } else if self.0.is_ascii() {
            // SAFETY: we just checked it is ASCII and ASCII is valid UTF8.
            write!(f, "{}", unsafe { std::str::from_utf8_unchecked(&self.0) })
        } else {
            write!(f, "{}", String::from_utf8_lossy(&self.0))
        }
    }
}

impl std::fmt::Debug for Cigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Cigar({})", String::from_utf8_lossy(&self.0))
    }
}

impl From<ExpandedCigar> for Cigar {
    #[inline]
    fn from(c: ExpandedCigar) -> Cigar {
        c.condense_to_cigar()
    }
}

impl TryFrom<String> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(s: String) -> Result<Self, Self::Error> {
        Cigar::try_from(s.into_bytes())
    }
}

impl TryFrom<Vec<u8>> for Cigar {
    type Error = CigarError;

    fn try_from(bytes: Vec<u8>) -> Result<Self, Self::Error> {
        let mut num: usize = 0;
        let mut has_number = false;

        if bytes == b"*" || bytes.is_empty() {
            return Ok(Cigar::new());
        }

        for b in &bytes {
            if b.is_ascii_digit() {
                num = num
                    .checked_mul(10)
                    .and_then(|n| n.checked_add((b - b'0') as usize))
                    .ok_or(CigarError::IncOverflow)?;
                has_number = true;
            } else if IS_CIGAR[*b as usize] {
                if !has_number {
                    return Err(CigarError::MissingInc);
                }
                if num == 0 {
                    return Err(CigarError::IncZero);
                }
                num = 0;
                has_number = false;
            } else {
                return Err(CigarError::InvalidOperation);
            }
        }

        if has_number {
            return Err(CigarError::MissingOp);
        }

        Ok(Cigar(bytes))
    }
}

impl<const N: usize> TryFrom<[u8; N]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: [u8; N]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_vec())
    }
}

impl TryFrom<&[u8]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &[u8]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_owned())
    }
}

impl TryFrom<&mut [u8]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_vec())
    }
}

impl<const N: usize> TryFrom<&[u8; N]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &[u8; N]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_vec())
    }
}

impl<const N: usize> TryFrom<&mut [u8; N]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8; N]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_vec())
    }
}

impl TryFrom<&str> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        Cigar::try_from(s.as_bytes())
    }
}

/// A single increment-opcode pair.
///
/// This is yielded when iterating over a [`Cigar`] string. One can also collect
/// an iterator of [`Ciglet`]s into a `Result<Cigar, CigarError>`.
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Ciglet {
    /// Increment or repetition count.
    pub inc: usize,
    /// Alignment operation code.
    pub op:  u8,
}

/// An iterator over a CIGAR string.
#[derive(Debug)]
pub struct CigletIterator<'a> {
    buffer: &'a [u8],
}

impl<'a> CigletIterator<'a> {
    /// Constructs an iterator from a CIGAR byte buffer.
    #[inline]
    fn new(buffer: &'a [u8]) -> Self {
        CigletIterator { buffer }
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
            index += 1;
        }

        Some(Ciglet { inc: num, op })
    }
}

impl std::fmt::Display for Ciglet {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}{}", self.inc, self.op as char)
    }
}

impl std::fmt::Debug for Ciglet {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({},{})", self.inc, self.op as char)
    }
}

impl FromIterator<Ciglet> for Result<Cigar, CigarError> {
    /// Converts an iterator of [`Ciglet`]s into a [`Cigar`] string, and checks
    /// that the result is valid.
    #[inline]
    fn from_iter<T: IntoIterator<Item = Ciglet>>(data: T) -> Self {
        let mut byte_string = Vec::new();
        let mut format_buffer = itoa::Buffer::new();

        for Ciglet { inc, op } in data {
            if inc == 0 {
                return Err(CigarError::IncZero);
            } else if !is_valid_op(op) {
                return Err(CigarError::InvalidOperation);
            }

            byte_string.extend_from_slice(format_buffer.format(inc).as_bytes());
            byte_string.push(op);
        }

        Ok(Cigar(byte_string))
    }
}

/// A CIGAR string as an "expanded" opcode-only byte string, e.g. `3M` → `MMM`.
#[derive(Clone, PartialEq)]
pub(crate) struct ExpandedCigar(Vec<u8>);

impl ExpandedCigar {
    /// Condenses the [`ExpandedCigar`] back to its standard form,  e.g. `MMM` → `3M`.
    #[inline]
    #[must_use]
    pub(crate) fn condense_to_cigar(self) -> Cigar {
        let mut condensed: Vec<u8> = Vec::new();
        let mut format_buffer = itoa::Buffer::new();

        let mut cigars = self.0.iter().copied().filter(|op| is_valid_op(*op));

        let Some(mut previous) = cigars.next() else {
            return Cigar(condensed);
        };

        let mut count: usize = 1;

        for op in cigars {
            if previous == op {
                count += 1;
            } else {
                condensed.extend_from_slice(format_buffer.format(count).as_bytes());
                condensed.push(previous);
                previous = op;
                count = 1;
            }
        }

        condensed.extend_from_slice(format_buffer.format(count).as_bytes());
        condensed.push(previous);

        Cigar(condensed)
    }
}

impl std::fmt::Display for ExpandedCigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}

impl std::fmt::Debug for ExpandedCigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}

impl From<Cigar> for ExpandedCigar {
    #[inline]
    fn from(c: Cigar) -> ExpandedCigar {
        c.expand_cigar()
    }
}

impl From<Vec<u8>> for ExpandedCigar {
    #[inline]
    fn from(vec: Vec<u8>) -> Self {
        ExpandedCigar(vec)
    }
}

impl From<&str> for ExpandedCigar {
    #[inline]
    fn from(s: &str) -> Self {
        ExpandedCigar(s.as_bytes().to_owned())
    }
}

impl From<&[u8]> for ExpandedCigar {
    #[inline]
    fn from(v: &[u8]) -> Self {
        ExpandedCigar(v.to_vec())
    }
}

impl<const N: usize> From<&[u8; N]> for ExpandedCigar {
    #[inline]
    fn from(v: &[u8; N]) -> Self {
        ExpandedCigar(v.to_vec())
    }
}

#[inline]
const fn is_valid_op(op: u8) -> bool {
    matches!(op, b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'=')
}
