#![allow(dead_code)]

use std::simd::{prelude::*, LaneCount, SimdElement, SupportedLaneCount};

use crate::data::{mappings::to_dna_profile_index, matrices::BiasedWeightMatrix, types::Uint};

pub struct QueryProfileDNA<T>
where
    T: Uint + From<u8>, {
    pub(crate) profile: Vec<T>,
    bias:               T,
}

impl<T> QueryProfileDNA<T>
where
    T: Uint + From<u8>,
{
    const ALPHA: &'static [u8] = b"ACGTN";
    #[must_use]
    pub fn new<const ALPHA_LEN: usize>(query: &[u8], matrix: &BiasedWeightMatrix<ALPHA_LEN>) -> Self {
        let mut profile = Vec::with_capacity(5 * query.len());

        for ref_index in Self::ALPHA.iter().copied().map(to_dna_profile_index) {
            for query_index in query.iter().copied().map(to_dna_profile_index) {
                profile.push(matrix.mapping[ref_index][query_index].into());
            }
        }

        let bias = matrix.bias.into();
        QueryProfileDNA { profile, bias }
    }
}

pub struct QueryProfileStripedDNA<T, const N: usize>
where
    T: SimdElement,
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) profile: Vec<Simd<T, N>>,
    pub(crate) bias:    T,
}

impl<T, const N: usize> QueryProfileStripedDNA<T, N>
where
    T: Uint + From<u8>,
    LaneCount<N>: SupportedLaneCount,
{
    const ALPHA: &'static [u8] = b"ACGTN";

    #[must_use]
    pub fn new<const ALPHA_LEN: usize>(query: &[u8], matrix: &BiasedWeightMatrix<ALPHA_LEN>) -> Self {
        // SupportedLaneCount cannot presently be zero.
        let number_vectors = (query.len() + (N - 1)) / N;
        let total_lanes = N * number_vectors;

        let zeroes = Simd::from_array([T::zero(); N]);
        let mut profile = vec![zeroes; 5 * number_vectors];

        for ref_index in Self::ALPHA.iter().copied().map(to_dna_profile_index) {
            for v in 0..number_vectors {
                let mut vector = [T::zero(); N];
                for (i, q) in (v..total_lanes).step_by(number_vectors).enumerate() {
                    if q < query.len() {
                        let query_index = to_dna_profile_index(query[q]);
                        vector[i] = matrix.mapping[ref_index][query_index].into();
                    }
                }
                profile[ref_index * number_vectors + v] = Simd::from_array(vector);
            }
        }

        let bias = matrix.bias.into();
        QueryProfileStripedDNA { profile, bias }
    }
}
