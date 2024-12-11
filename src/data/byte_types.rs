#![allow(dead_code)]

use super::mappings::DNA_PROFILE_MAP;

pub(crate) trait ByteMappings {
    fn to_dna_index(self) -> usize;
}

impl ByteMappings for u8 {
    #[inline]
    fn to_dna_index(self) -> usize {
        DNA_PROFILE_MAP[self] as usize
    }
}

pub(crate) trait ByteAvg {
    fn avg(self, other: u8) -> u8;
}

#[allow(clippy::cast_possible_truncation)]
impl ByteAvg for u8 {
    fn avg(self, other: u8) -> u8 {
        ((u16::from(self) + u16::from(other)) / 2u16) as u8
    }
}

#[allow(clippy::wrong_self_convention)]
pub(crate) trait IsBase {
    fn is_base_acgt(self) -> bool;
    fn is_not_base_acgt(self) -> bool;
}

impl IsBase for u8 {
    #[inline]
    fn is_base_acgt(self) -> bool {
        self == b'A' || self == b'G' || self == b'T' || self == b'C'
    }

    #[inline]
    fn is_not_base_acgt(self) -> bool {
        !self.is_base_acgt()
    }
}
