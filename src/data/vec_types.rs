use crate::alignment::pairwise_align_with_cigar;
use crate::data::types::cigar::Cigar;
use crate::prelude::{AminoAcids, Nucleotides};
use crate::simd::SimdByteFunctions;
use std::simd::{prelude::*, LaneCount, SupportedLaneCount};

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

pub trait CheckSequence {
    fn is_ascii_simd<const N: usize>(&self) -> bool
    where
        LaneCount<N>: SupportedLaneCount;

    fn is_graphic_simd<const N: usize>(&self) -> bool
    where
        LaneCount<N>: SupportedLaneCount;
}

impl<T> CheckSequence for T
where
    T: AsRef<[u8]>,
{
    #[inline]
    fn is_ascii_simd<const N: usize>(&self) -> bool
    where
        LaneCount<N>: SupportedLaneCount, {
        let (pre, mid, suffix) = self.as_ref().as_simd();
        pre.is_ascii() && suffix.is_ascii() && mid.iter().fold(Mask::splat(true), |acc, b| acc & b.is_ascii()).all()
    }

    #[inline]
    fn is_graphic_simd<const N: usize>(&self) -> bool
    where
        LaneCount<N>: SupportedLaneCount, {
        let (pre, mid, suffix) = self.as_ref().as_simd();
        pre.iter().fold(true, |acc, b| acc & b.is_ascii_graphic())
            && suffix.iter().fold(true, |acc, b| acc & b.is_ascii_graphic())
            && mid.iter().fold(Mask::splat(true), |acc, b| acc & b.is_ascii_graphic()).all()
    }
}

#[cfg(test)]
mod test {
    use super::CheckSequence;
    use crate::data::alphas::AMINO_ACIDS;

    #[test]
    fn is_ascii() {
        let s = crate::generate::rand_sequence(AMINO_ACIDS, 151, 42);
        assert_eq!(s.is_ascii(), s.is_ascii_simd::<16>());
    }
}

#[cfg(test)]
mod bench {
    use super::CheckSequence;
    use crate::data::alphas::AMINO_ACIDS;
    use std::sync::LazyLock;
    use test::Bencher;
    extern crate test;

    const N: usize = 151;
    const SEED: u64 = 99;

    static SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(AMINO_ACIDS, N, SEED));

    #[bench]
    fn is_ascii_std(b: &mut Bencher) {
        b.iter(|| SEQ.is_ascii());
    }

    #[bench]
    fn is_ascii_zoe(b: &mut Bencher) {
        let (p, m, s) = SEQ.as_simd::<16>();
        eprintln!("{p} {m} {s}", p = p.len(), m = m.len(), s = s.len());
        b.iter(|| SEQ.is_ascii_simd::<16>());
    }
}
