#![allow(dead_code)]

use crate::data::DNA_PROFILE_MAP;

pub(crate) trait ByteMappings {
    fn to_dna_index(self) -> usize;
}

impl ByteMappings for u8 {
    #[inline]
    fn to_dna_index(self) -> usize {
        DNA_PROFILE_MAP[self] as usize
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
