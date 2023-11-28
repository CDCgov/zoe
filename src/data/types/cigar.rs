use atoi::FromRadix10Checked;

#[derive(Clone, PartialEq)]
pub struct Cigar(Vec<u8>);

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
#[derive(Clone, Copy)]
pub struct Ciglet {
    pub inc: usize,
    pub op:  u8,
}

pub struct CigletIterator<'a> {
    buffer: &'a [u8],
    offset: usize,
}

impl<'a> CigletIterator<'a> {
    #[inline]
    fn new(buffer: &'a [u8]) -> Self {
        CigletIterator { buffer, offset: 0 }
    }
}

impl<'a> Iterator for CigletIterator<'a> {
    type Item = Ciglet;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.buffer.get(self.offset..)
            && let (Some(count), width) = usize::from_radix_10_checked(s)
            && let Some(&state) = self.buffer.get(self.offset + width)
        {
            self.offset += width + 1;
            if matches!(state, b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'=') {
                Some(Ciglet { inc: count, op: state })
            } else {
                self.next()
            }
        } else {
            None
        }
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

#[derive(Clone, Debug)]
pub(crate) struct AlignmentStates(Vec<Ciglet>);

#[allow(dead_code)]
impl AlignmentStates {
    #[must_use]
    #[inline]
    pub(crate) fn new() -> Self {
        AlignmentStates(Vec::new())
    }

    pub(crate) fn add_state(&mut self, op: u8) {
        if let Some(c) = self.0.last_mut()
            && c.op == op
        {
            c.inc += 1;
        } else {
            self.0.push(Ciglet { inc: 1, op });
        }
    }

    pub(crate) fn to_cigar(&self) -> Cigar {
        let mut condensed = Vec::new();
        let mut format_buffer = itoa::Buffer::new();

        for Ciglet { inc, op } in self.0.iter().copied() {
            condensed.extend_from_slice(format_buffer.format(inc).as_bytes());
            condensed.push(op);
        }

        Cigar(condensed)
    }
}

#[derive(PartialEq, Clone)]
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

    #[test]
    fn test_expand() {
        let cigar: Cigar = " 4S10M2I2D3M4H4P".into();
        let expanded: ExpandedCigar = "SSSSMMMMMMMMMMIIDDMMMHHHHPPPP".into();
        assert_eq!(cigar.expand_cigar(), expanded);
    }

    #[test]
    fn test_match_length() {
        let cigars = [
            ("4S10M2I2D3M4H4P", 15),
            ("", 0),
            ("3M2D1M", 6),
            ("255M", 255),
            ("3M1D4I8X9=4M", 25),
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
