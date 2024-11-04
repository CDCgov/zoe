/// A CIGAR string.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Cigar(pub(crate) Vec<u8>);

#[cfg(target_pointer_width = "16")]
const USIZE_WIDTH: usize = 5;
#[cfg(target_pointer_width = "32")]
const USIZE_WIDTH: usize = 10;
#[cfg(target_pointer_width = "64")]
const USIZE_WIDTH: usize = 20;

impl Cigar {
    #[must_use]
    pub(crate) fn expand_cigar(&self) -> ExpandedCigar {
        let mut expanded = Vec::new();

        for Ciglet { inc, op } in self.into_iter() {
            expanded.extend(std::iter::repeat(op).take(inc));
        }

        ExpandedCigar(expanded)
    }

    #[must_use]
    pub fn match_length(&self) -> usize {
        self.into_iter()
            .filter(|Ciglet { inc: _, op }| matches!(op, b'M' | b'D' | b'N' | b'X' | b'='))
            .map(|Ciglet { inc, .. }| inc)
            .sum()
    }

    pub fn into_iter(&self) -> impl Iterator<Item = Ciglet> + '_ {
        CigletIterator::new(&self.0)
    }
}

use core::str;
use std::fmt;
impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))?;
        Ok(())
    }
}

impl fmt::Debug for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))?;
        Ok(())
    }
}

impl From<&str> for Cigar {
    fn from(s: &str) -> Self {
        Cigar(s.as_bytes().to_owned())
    }
}

impl From<&[u8]> for Cigar {
    fn from(v: &[u8]) -> Self {
        Cigar(v.to_vec())
    }
}

impl From<ExpandedCigar> for Cigar {
    fn from(c: ExpandedCigar) -> Cigar {
        c.condense_to_cigar()
    }
}

/// A single increment quantifier-operation pair.
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Ciglet {
    pub inc: usize,
    pub op:  u8,
}

/// An iterator over a CIGAR string.
#[derive(Debug)]
pub struct CigletIterator<'a> {
    buffer: &'a [u8],
}

impl<'a> CigletIterator<'a> {
    #[inline]
    fn new(buffer: &'a [u8]) -> Self {
        CigletIterator {
            buffer: buffer.trim_ascii_start(),
        }
    }
}

impl<'a> Iterator for CigletIterator<'a> {
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
            } else if matches!(b, b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'=') && num > 0 {
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
            } else if matches!(b, b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'=') && num > 0 {
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

impl fmt::Display for Ciglet {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}", self.inc, self.op as char)?;
        Ok(())
    }
}

impl fmt::Debug for Ciglet {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({},{})", self.inc, self.op as char)?;
        Ok(())
    }
}

#[derive(Clone, PartialEq)]
pub(crate) struct ExpandedCigar(Vec<u8>);

impl ExpandedCigar {
    #[inline]
    #[must_use]
    pub(crate) fn condense_to_cigar(self) -> Cigar {
        let mut condensed: Vec<u8> = Vec::new();
        let mut format_buffer = itoa::Buffer::new();

        let mut cigars = self
            .0
            .iter()
            .copied()
            .filter(|op| matches!(op, b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'='));

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

impl fmt::Display for ExpandedCigar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))?;
        Ok(())
    }
}

impl fmt::Debug for ExpandedCigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))?;
        Ok(())
    }
}

impl From<Cigar> for ExpandedCigar {
    fn from(c: Cigar) -> ExpandedCigar {
        c.expand_cigar()
    }
}

impl From<Vec<u8>> for ExpandedCigar {
    fn from(vec: Vec<u8>) -> Self {
        ExpandedCigar(vec)
    }
}

impl From<&str> for ExpandedCigar {
    fn from(s: &str) -> Self {
        ExpandedCigar(s.as_bytes().to_owned())
    }
}

impl From<&[u8]> for ExpandedCigar {
    fn from(v: &[u8]) -> Self {
        ExpandedCigar(v.to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::AlignmentStates;

    #[test]
    fn test_expand() {
        let cigar: Cigar = " 4S10M2I2D3M4H4P".into();
        let expanded: ExpandedCigar = "SSSSMMMMMMMMMMIIDDMMMHHHHPPPP".into();
        assert_eq!(cigar.expand_cigar(), expanded);
    }

    #[test]
    fn test_iter() {
        let cigar: Cigar = "  1M22I333D4444N55555S6666H777P88X9=  ".into();
        let mut cigar = cigar.into_iter();

        assert_eq!(cigar.next(), Some(Ciglet { op: b'M', inc: 1 }));
        assert_eq!(cigar.next(), Some(Ciglet { op: b'I', inc: 22 }));
        assert_eq!(cigar.next(), Some(Ciglet { op: b'D', inc: 333 }));
        assert_eq!(cigar.next(), Some(Ciglet { op: b'N', inc: 4444 }));
        assert_eq!(cigar.next(), Some(Ciglet { op: b'S', inc: 55555 }));
        assert_eq!(cigar.next(), Some(Ciglet { op: b'H', inc: 6666 }));
        assert_eq!(cigar.next(), Some(Ciglet { op: b'P', inc: 777 }));
        assert_eq!(cigar.next(), Some(Ciglet { op: b'X', inc: 88 }));
        assert_eq!(cigar.next(), Some(Ciglet { op: b'=', inc: 9 }));
        assert_eq!(cigar.next(), None);

        // Illegal cigar op
        let cigar: Cigar = "8K".into();
        assert_eq!(cigar.into_iter().next(), None);

        // Leading zeroes don't matter
        let cigar: Cigar = "000000000000000000000000000000155M".into();
        assert_eq!(cigar.into_iter().next(), Some(Ciglet { op: b'M', inc: 155 }));

        // Overflows
        let cigar: Cigar = "100000000000000000000000000000155M".into();
        assert_eq!(cigar.into_iter().next(), None);

        // Bad order
        let cigar: Cigar = "M155M".into();
        assert_eq!(cigar.into_iter().next(), None);

        // usize == u64
        if USIZE_WIDTH == 20 {
            let cigar: Cigar = "18446744073709551615M".into();
            assert_eq!(
                cigar.into_iter().next(),
                Some(Ciglet {
                    op:  b'M',
                    inc: 18_446_744_073_709_551_615,
                })
            );

            let cigar: Cigar = "001234567890123456789M".into();
            assert_eq!(
                cigar.into_iter().next(),
                Some(Ciglet {
                    op:  b'M',
                    inc: 1_234_567_890_123_456_789,
                })
            );
        }
    }

    #[test]
    fn test_match_length() {
        let cigars = [
            ("4S10M2I2D3M4H4P", 15),
            ("", 0),
            ("3M2D1M", 6),
            ("255M", 255),
            ("3M1D4I8X9=4M", 25),
            ("M", 0),
        ];

        for (c, l) in cigars {
            let cigar = Cigar::from(c);
            assert_eq!(cigar.match_length(), l);
        }
    }

    #[test]
    fn test_condense_cigar() {
        let cigars = ["4S10M2I2D3M4H4P", "", "3M2D1M", "255M", "3M1D4I8X9=4M"];

        for c in cigars {
            let cigar = Cigar::from(c);
            assert_eq!(cigar.expand_cigar().condense_to_cigar(), cigar);
        }
    }

    #[test]
    fn test_add_state() {
        let cigar: Cigar = "4S10M2I2D3M4H4P".into();
        let mut states = AlignmentStates::new();
        for Ciglet { inc, op } in cigar.into_iter() {
            for _ in 0..inc {
                states.add_state(op);
            }
        }

        assert_eq!(cigar, states.to_cigar());
    }
}

#[cfg(test)]
mod benches {
    use super::*;
    use test::Bencher;
    extern crate test;

    #[bench]
    fn expand_cigar(b: &mut Bencher) {
        let cigar: Cigar = "3S10M2I2D3M4H4P".into();

        b.iter(|| cigar.expand_cigar());
    }

    #[bench]
    fn condense_cigar(b: &mut Bencher) {
        let cigar = Cigar::from("4S10M2I2D3M4H4P").expand_cigar();
        b.iter(|| cigar.clone().condense_to_cigar());
    }

    #[bench]
    fn match_length(b: &mut Bencher) {
        let cigar = Cigar::from("4S10M2I2D3M4H4P");
        b.iter(|| cigar.match_length());
    }
}
