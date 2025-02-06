use crate::{
    alignment::sw::{sw_scalar_score, sw_simd_score},
    data::{WeightMatrix, err::QueryProfileError, mappings::ByteIndexMap, types::cigar::Cigar},
    math::{AnyInt, FromSameSignedness},
};
use std::{
    convert::Into,
    simd::{LaneCount, SimdElement, SupportedLaneCount, prelude::*},
    vec,
};

use super::sw::sw_scalar_alignment;

#[inline]
pub(crate) fn validate_profile_args<Q: AsRef<[u8]>, U: AnyInt>(
    query: Q, gap_open: U, gap_extend: U,
) -> Result<(), QueryProfileError> {
    if query.as_ref().is_empty() {
        Err(QueryProfileError::EmptyQuery)
    } else if gap_extend > gap_open {
        Err(QueryProfileError::BadGapWeights)
    } else {
        Ok(())
    }
}

/// A profile for DNA sequence alignment using scalar operations.
///
/// The profile stores the query, weight matrix, and gap penalties to be used in
/// an alignment. The API mirrors that of [`StripedProfile`], but this profile
/// does not use SIMD or create a striped layout.
///
/// # Type Parameters
/// * `S` - The size of the alphabet (usually 5 for DNA including *N*)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScalarProfile<'a, const S: usize> {
    pub(crate) query:      &'a [u8],
    pub(crate) matrix:     WeightMatrix<i8, S>,
    pub(crate) gap_open:   i32,
    pub(crate) gap_extend: i32,
}

impl<'a, const S: usize> ScalarProfile<'a, S> {
    /// Create a new profile for use with the scalar alignment algorithm.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    pub fn new<Q: AsRef<[u8]> + ?Sized>(
        query: &'a Q, matrix: WeightMatrix<i8, S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, QueryProfileError> {
        validate_profile_args(query, gap_open, gap_extend)?;

        Ok(ScalarProfile {
            query: query.as_ref(),
            matrix,
            gap_open: -i32::from(gap_open),
            gap_extend: -i32::from(gap_extend),
        })
    }

    /// Computes the Smith-Waterman local alignment score between the profile
    /// and a passed sequence.
    ///
    /// For more info, see: [`sw_scalar_score`].
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{ScalarProfile, sw::sw_scalar_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"GGCCACAGGATTGAG";
    /// let query: &[u8] = b"CTCAGATTG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = ScalarProfile::<5>::new(query, WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score(reference);
    /// assert_eq!(score, 27);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8]) -> i32 {
        sw_scalar_score(seq, self)
    }

    /// Computes the Smith-Waterman local alignment between the profile and a
    /// passed sequence, yielding the reference starting position (indexed from
    /// 1), cigar, and optimal score.
    ///
    /// For more info, see: [`sw_scalar_alignment`].
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{ScalarProfile, sw::sw_scalar_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"GGCCACAGGATTGAG";
    /// let query: &[u8] = b"CTCAGATTG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = ScalarProfile::<5>::new(query, WEIGHTS, 3, 1).unwrap();
    /// let (start, cigar, score) = profile.smith_waterman_alignment(reference);
    /// assert_eq!(start, 4);
    /// assert_eq!(cigar, b"5M1D4M".into());
    /// assert_eq!(score, 27);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_alignment(&self, seq: &[u8]) -> (usize, Cigar, i32) {
        sw_scalar_alignment(seq, self)
    }
}

/// A striped profile for DNA sequence alignment using SIMD operations.
///
/// The profile contains pre-computed scoring vectors for each residue in the
/// specified mapping, arranged in a striped pattern to optimize SIMD operations
/// during alignment.
///
/// # Type Parameters
/// * `T` - The numeric type used for scores (u8, u16, u32, or u64)
/// * `N` - The number of SIMD lanes (usually 16, 32 or 64)
/// * `S` - The size of the alphabet (usually 5 for DNA including *N*)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct StripedProfile<T, const N: usize, const S: usize>
where
    T: SimdElement,
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) profile:    Vec<Simd<T, N>>,
    pub(crate) gap_open:   T,
    pub(crate) gap_extend: T,
    pub(crate) bias:       T,
    pub(crate) mapping:    &'static ByteIndexMap<S>,
}

impl<T, const N: usize, const S: usize> StripedProfile<T, N, S>
where
    T: AnyInt + SimdElement,
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates a new striped profile from a sequence and scoring matrix.
    ///
    /// See: [`WeightMatrix`]
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    pub fn new<U>(query: &[u8], matrix: &WeightMatrix<U, S>, gap_open: U, gap_extend: U) -> Result<Self, QueryProfileError>
    where
        T: FromSameSignedness<U>,
        U: AnyInt, {
        validate_profile_args(query, gap_open, gap_extend)?;
        Ok(Self::new_unchecked(query, matrix, gap_open, gap_extend))
    }

    /// Creates a new striped profile from a sequence and scoring matrix.
    ///
    /// See: [`WeightMatrix`]
    pub(crate) fn new_unchecked<U>(query: &[u8], matrix: &WeightMatrix<U, S>, gap_open: U, gap_extend: U) -> Self
    where
        T: From<U>,
        U: AnyInt, {
        // SupportedLaneCount cannot presently be zero.
        let number_vectors = query.len().div_ceil(N);
        let total_lanes = N * number_vectors;

        let bias = matrix.bias.into();
        let biases = Simd::splat(bias);
        let mut profile = vec![biases; S * number_vectors];

        for ref_index in 0..matrix.mapping.len() {
            for v in 0..number_vectors {
                let mut vector = biases;
                for (i, q) in (v..total_lanes).step_by(number_vectors).enumerate() {
                    if q < query.len() {
                        let query_index = matrix.mapping.to_index(query[q]);
                        vector[i] = matrix.weights[ref_index][query_index].into();
                    }
                }
                profile[ref_index * number_vectors + v] = vector;
            }
        }

        StripedProfile {
            profile,
            gap_open: gap_open.into(),
            gap_extend: gap_extend.into(),
            bias,
            mapping: matrix.mapping,
        }
    }

    /// Returns the number of SIMD vectors in the profile.
    #[must_use]
    #[inline]
    pub fn number_vectors(&self) -> usize {
        self.profile.len() / S
    }
}

impl<const N: usize, const S: usize> StripedProfile<u8, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Computes the Smith-Waterman local alignment score between the `u8`
    /// profile and a passed sequence. Returns [`None`] if the score overflowed.
    ///
    /// For more info, see: [`sw_simd_score`].
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{StripedProfile, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = StripedProfile::<u8, 32, 5>::new(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8]) -> Option<u64> {
        sw_simd_score::<u8, N, S>(seq, self).map(Into::into)
    }
}

impl<const N: usize, const S: usize> StripedProfile<u16, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Computes the Smith-Waterman local alignment score between the `u16`
    /// profile and a passed sequence. Returns [`None`] if the score overflowed.
    ///
    /// For more info, see: [`sw_simd_score`].
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{StripedProfile, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = StripedProfile::<u16, 32, 5>::new(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8]) -> Option<u64> {
        sw_simd_score::<u16, N, S>(seq, self).map(Into::into)
    }
}

impl<const N: usize, const S: usize> StripedProfile<u32, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Computes the Smith-Waterman local alignment score between the `u32`
    /// profile and a passed sequence. Returns [`None`] if the score overflowed.
    ///
    /// For more info, see: [`sw_simd_score`].
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{StripedProfile, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = StripedProfile::<u32, 32, 5>::new(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8]) -> Option<u64> {
        sw_simd_score::<u32, N, S>(seq, self).map(Into::into)
    }
}

impl<const N: usize, const S: usize> StripedProfile<u64, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Computes the Smith-Waterman local alignment score between the `u64`
    /// profile and a passed sequence. Returns [`None`] if the score overflowed.
    ///
    /// For more info, see: [`sw_simd_score`].
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{StripedProfile, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = StripedProfile::<u64, 32, 5>::new(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, seq: &[u8]) -> Option<u64> {
        sw_simd_score::<u64, N, S>(seq, self)
    }
}

#[cfg(test)]
mod bench {
    use test::Bencher;
    extern crate test;
    use super::*;
    use crate::alignment::sw::test_data::{GAP_EXTEND, GAP_OPEN};
    use crate::data::{constants::matrices::WeightMatrix, mappings::DNA_PROFILE_MAP};
    pub(crate) static DATA: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/KJ907631.1.txt")); // H5 HA, complete CDS
    pub(crate) static MATRIX: WeightMatrix<u8, 5> =
        WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N')).into_biased_matrix();

    #[bench]
    fn build_profile_u8(b: &mut Bencher) {
        b.iter(|| StripedProfile::<u8, 32, 5>::new(DATA, &MATRIX, GAP_OPEN, GAP_EXTEND));
    }
}
