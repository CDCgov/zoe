use crate::{
    alignment::sw::sw_simd_score,
    data::{mappings::to_dna_profile_index, matrices::BiasedWeightMatrix, types::Uint},
};
use std::{
    convert::Into,
    simd::{prelude::*, LaneCount, SimdElement, SupportedLaneCount},
    vec,
};

/// A striped profile for DNA sequence alignment using SIMD operations.
///
/// The profile contains pre-computed scoring vectors for each canonical DNA
/// base (`A`,`C`,`G`,`T`) + `N` arranged in a striped pattern to optimize SIMD
/// operations during alignment.
///
/// # Type Parameters
/// * `T` - The numeric type used for scores (u8, u16, u32, or u64)
/// * `N` - The number of SIMD lanes (usually 16, 32 or 64)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct StripedDNAProfile<T, const N: usize>
where
    T: SimdElement,
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) profile: Vec<Simd<T, N>>,
    pub(crate) bias:    T,
}

impl<T, const N: usize> StripedDNAProfile<T, N>
where
    T: Uint,
    LaneCount<N>: SupportedLaneCount,
{
    const ALPHA: &'static [u8] = b"ACGTN";

    /// Creates a new striped DNA profile from a sequence and scoring matrix.
    ///
    /// See: [`BiasedWeightMatrix`]
    #[must_use]
    pub fn new<const ALPHA_LEN: usize>(seq: &[u8], matrix: &BiasedWeightMatrix<ALPHA_LEN>) -> Self
    where
        T: From<u8>, {
        // SupportedLaneCount cannot presently be zero.
        let number_vectors = seq.len().div_ceil(N);
        let total_lanes = N * number_vectors;

        let bias: T = matrix.bias.into();
        let biases = Simd::splat(bias);
        let mut profile = vec![biases; Self::ALPHA.len() * number_vectors];

        for ref_index in Self::ALPHA.iter().copied().map(to_dna_profile_index) {
            for v in 0..number_vectors {
                let mut vector = biases;
                for (i, q) in (v..total_lanes).step_by(number_vectors).enumerate() {
                    if q < seq.len() {
                        let query_index = to_dna_profile_index(seq[q]);
                        vector[i] = matrix.mapping[ref_index][query_index].into();
                    }
                }
                profile[ref_index * number_vectors + v] = vector;
            }
        }

        StripedDNAProfile { profile, bias }
    }

    /// Returns the number of SIMD vectors in the profile.
    #[must_use]
    #[inline]
    pub fn number_vectors(&self) -> usize {
        self.profile.len() / Self::ALPHA.len()
    }
}

impl<const N: usize> StripedDNAProfile<u8, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Computes the Smith-Waterman alignment score between the `u8` profile and
    /// a sequence. Returns [`None`] if the score overflowed.
    ///
    /// # Arguments
    /// * `seq` - The sequence to align against the profile
    /// * `gap_open` - The gap opening penalty (the unsigned integer will be subtracted)
    /// * `gap_extend` - The gap extension penalty (the unsigned integer will be subtracted)
    ///
    /// For more info, see: [`sw_simd_score`].
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8], gap_open: u8, gap_extend: u8) -> Option<u64> {
        sw_simd_score::<u8, N, _>(seq, self, gap_open, gap_extend).map(Into::into)
    }
}

impl<const N: usize> StripedDNAProfile<u16, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Computes the Smith-Waterman alignment score between the `u16` profile and
    /// a sequence. Returns [`None`] if the score overflowed.
    ///
    /// # Arguments
    /// * `seq` - The sequence to align against the profile
    /// * `gap_open` - The gap opening penalty (the unsigned integer will be subtracted)
    /// * `gap_extend` - The gap extension penalty (the unsigned integer will be subtracted)
    ///
    /// For more info, see: [`sw_simd_score`].
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8], gap_open: u8, gap_extend: u8) -> Option<u64> {
        sw_simd_score::<u16, N, _>(seq, self, gap_open.into(), gap_extend.into()).map(Into::into)
    }
}

impl<const N: usize> StripedDNAProfile<u32, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Computes the Smith-Waterman alignment score between the `u32` profile and
    /// a sequence. Returns [`None`] if the score overflowed.
    ///
    /// # Arguments
    /// * `seq` - The sequence to align against the profile
    /// * `gap_open` - The gap opening penalty (the unsigned integer will be subtracted)
    /// * `gap_extend` - The gap extension penalty (the unsigned integer will be subtracted)
    ///
    /// For more info, see: [`sw_simd_score`].
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8], gap_open: u8, gap_extend: u8) -> Option<u64> {
        sw_simd_score::<u32, N, _>(seq, self, gap_open.into(), gap_extend.into()).map(Into::into)
    }
}

impl<const N: usize> StripedDNAProfile<u64, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Computes the Smith-Waterman alignment score between the `u64` profile and
    /// a sequence. Returns [`None`] if the score overflowed.
    ///
    /// # Arguments
    /// * `seq` - The sequence to align against the profile
    /// * `gap_open` - The gap opening penalty (the unsigned integer will be subtracted)
    /// * `gap_extend` - The gap extension penalty (the unsigned integer will be subtracted)
    ///
    /// For more info, see: [`sw_simd_score`].
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8], gap_open: u8, gap_extend: u8) -> Option<u64> {
        sw_simd_score::<u64, N, _>(seq, self, gap_open.into(), gap_extend.into())
    }
}

#[cfg(test)]
mod bench {
    use test::Bencher;
    extern crate test;
    use super::*;
    use crate::data::constants::matrices::SimpleWeightMatrix;
    pub(crate) static DATA: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/KJ907631.1.txt")); // H5 HA, complete CDS
    pub(crate) static MATRIX: BiasedWeightMatrix<5> =
        SimpleWeightMatrix::new(b"ACGTN", 2, -5, Some(b'N')).into_biased_matrix();

    #[bench]
    fn build_profile_u8(b: &mut Bencher) {
        b.iter(|| StripedDNAProfile::<u8, 32>::new(DATA, &MATRIX));
    }
}
