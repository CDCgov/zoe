//! A private module for helper functions on [`u8`].

use crate::{
    data::DNA_PROFILE_MAP,
    prelude::{IsValidDNA, RecodeDNAStrat, RefineDNAStrat},
};

/// An extension trait allowing for bytes to be converted to indices with
/// [`DNA_PROFILE_MAP`].
pub(crate) trait ByteMappings {
    /// Converts the byte to an index using [`DNA_PROFILE_MAP`].
    fn to_dna_index(self) -> usize;
}

impl ByteMappings for u8 {
    #[inline]
    fn to_dna_index(self) -> usize {
        DNA_PROFILE_MAP[self] as usize
    }
}

/// Extension trait to verify whether a byte is uppercase `ACGT` or not.
pub(crate) trait IsBase: Copy {
    /// Returns `true` if the base is uppercase `ACGT`.
    fn is_base_acgt(self) -> bool;

    /// Returns `true` if the base is not uppercase `ACGT`.
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

/// A trait allowing for validation, recoding, and refinement of a single DNA
/// base.
///
/// Similar functionality for entire sequences is in [`CheckNucleotides`],
/// [`RecodeNucleotides`], and [`RetainNucleotides`]. This trait can be useful
/// for more customized scenarios, such as recoding an array/codon via `map`.
///
/// [`CheckNucleotides`]: crate::data::types::nucleotides::CheckNucleotides
/// [`RecodeNucleotides`]: crate::data::types::nucleotides::RecodeNucleotides
/// [`RetainNucleotides`]: crate::data::types::nucleotides::RetainNucleotides
pub trait SanitizeBase: Sized + Copy {
    /// Checks whether a single byte is valid under the given validation
    /// strategy
    #[must_use]
    fn is_valid(self, strategy: IsValidDNA) -> bool;

    /// Recodes a single byte using the given strategy.
    #[must_use]
    fn recode_base(self, strategy: RecodeDNAStrat) -> u8;

    /// Refines and recodes a single byte using the given strategy.
    #[must_use]
    fn refine_base(self, strategy: RefineDNAStrat) -> Option<u8>;
}

impl SanitizeBase for u8 {
    #[inline]
    fn is_valid(self, strategy: IsValidDNA) -> bool {
        strategy.mapping()[self as usize]
    }

    #[inline]
    fn recode_base(self, strategy: RecodeDNAStrat) -> u8 {
        strategy.mapping()[self]
    }

    #[inline]
    fn refine_base(self, strategy: RefineDNAStrat) -> Option<u8> {
        let recoded = strategy.mapping()[self];
        if recoded > 0 { Some(recoded) } else { None }
    }
}
