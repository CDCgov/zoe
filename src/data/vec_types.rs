use crate::alignment::pairwise_align_with_cigar;
use crate::data::types::cigar::Cigar;
use crate::prelude::{AminoAcids, Nucleotides};

/// Utility trait for sequence validation using byte-wise mappings.
pub trait ValidateSequence {
    fn retain_by_validation(&mut self, validation_mapping: [bool; 256]);
    fn retain_by_recoding(&mut self, transformation_mapping: [u8; 256]);
    fn recode(&mut self, transformation_mapping: [u8; 256]);
}

impl ValidateSequence for Vec<u8> {
    /// Allows for filtering of biological sequences using a validation mapping.
    #[inline]
    fn retain_by_validation(&mut self, validation_mapping: [bool; 256]) {
        self.retain(|b| validation_mapping[*b as usize]);
    }

    /// Allows for filtering of biological sequences using byte re-encoding. The
    /// 0-byte is assumed to be an invalid pattern.
    #[inline]
    fn retain_by_recoding(&mut self, transformation_mapping: [u8; 256]) {
        self.retain_mut(|b| {
            *b = transformation_mapping[*b as usize];
            *b > 0
        });
    }

    /// Allows for re-coding the sequence to the allowed values in-place.
    fn recode(&mut self, transformation_mapping: [u8; 256]) {
        for b in self {
            *b = transformation_mapping[*b as usize];
        }
    }
}

/// A trait for expanding a sequence to its aligned state.
pub trait PairwiseSequence {
    type Output;
    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output);
}

impl PairwiseSequence for Vec<u8> {
    type Output = Self;

    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        pairwise_align_with_cigar(self, query, cigar, position)
    }
}

impl PairwiseSequence for &[u8] {
    type Output = Vec<u8>;

    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        pairwise_align_with_cigar(self, query, cigar, position)
    }
}

impl PairwiseSequence for Nucleotides {
    type Output = Self;

    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        let (r, q) = pairwise_align_with_cigar(self.as_bytes(), query.as_bytes(), cigar, position);
        (r.into(), q.into())
    }
}

impl PairwiseSequence for AminoAcids {
    type Output = Self;

    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        let (r, q) = pairwise_align_with_cigar(self.as_bytes(), query.as_bytes(), cigar, position);
        (r.into(), q.into())
    }
}
