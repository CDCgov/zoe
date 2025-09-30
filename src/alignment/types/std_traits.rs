use super::AlignmentStates;
use crate::data::types::cigar::{Cigar, CigarError, Ciglet, CigletIteratorChecked};
use std::fmt::Write;

impl AsRef<[Ciglet]> for AlignmentStates {
    #[inline]
    fn as_ref(&self) -> &[Ciglet] {
        &self.0
    }
}

impl Default for AlignmentStates {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl PartialEq<Cigar> for AlignmentStates {
    #[inline]
    fn eq(&self, other: &Cigar) -> bool {
        let mut o = other.iter();
        let matches = self.iter().copied().eq(o.by_ref());
        matches && o.valid()
    }
}

impl<'a> IntoIterator for &'a AlignmentStates {
    type Item = Ciglet;
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, Ciglet>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter().copied()
    }
}

impl IntoIterator for AlignmentStates {
    type Item = Ciglet;
    type IntoIter = <Vec<Ciglet> as IntoIterator>::IntoIter;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl FromIterator<Ciglet> for AlignmentStates {
    #[inline]
    fn from_iter<T: IntoIterator<Item = Ciglet>>(iter: T) -> Self {
        AlignmentStates(iter.into_iter().collect::<Vec<_>>())
    }
}

impl std::fmt::Debug for AlignmentStates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.0.iter().map(|c| (c.inc, c.op as char))).finish()
    }
}

impl std::fmt::Display for AlignmentStates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buff = itoa::Buffer::new();
        for Ciglet { inc, op } in self {
            f.write_str(buff.format(inc))?;
            f.write_char(op as char)?;
        }
        Ok(())
    }
}

impl TryFrom<&[u8]> for AlignmentStates {
    type Error = CigarError;

    fn try_from(v: &[u8]) -> Result<Self, Self::Error> {
        if v == b"*" {
            Ok(AlignmentStates::new())
        } else {
            CigletIteratorChecked::new(v)
                .collect::<Result<Vec<_>, _>>()
                .map(AlignmentStates)
        }
    }
}

impl TryFrom<String> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(s: String) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(s.as_bytes())
    }
}

impl TryFrom<Vec<u8>> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(bytes: Vec<u8>) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(bytes.as_slice())
    }
}

impl<const N: usize> TryFrom<[u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: [u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl TryFrom<&mut [u8]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(&*v)
    }
}

impl<const N: usize> TryFrom<&[u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &[u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl<const N: usize> TryFrom<&mut [u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl TryFrom<&str> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(s.as_bytes())
    }
}
