#[cfg(feature = "dev-3pass")]
use crate::alignment::sw::{sw_banded_alignment, sw_simd_score_ranges};
use crate::{
    alignment::{
        Alignment, MaybeAligned,
        sw::{sw_scalar_alignment, sw_scalar_score, sw_simd_alignment, sw_simd_score, sw_simd_score_ends},
    },
    data::{err::QueryProfileError, mappings::ByteIndexMap, matrices::WeightMatrix},
    math::{AlignableIntWidth, AnyInt, FromSameSignedness},
    simd::SimdAnyInt,
};
use std::{
    convert::Into,
    simd::{LaneCount, SimdElement, SupportedLaneCount, prelude::*},
    vec,
};

/// Validate the arguments for [`ScalarProfile`] or [`StripedProfile`].
///
/// ## Errors
///
/// The following errors are possible:
/// - [`QueryProfileError::EmptyQuery`] if `query` is empty
/// - [`QueryProfileError::GapOpenOutOfRange`] if `gap_open` is not between -127
///   and 0, inclusive
/// - [`QueryProfileError::GapExtendOutOfRange`] if `gap_extend` is not between
///   -127 and 0, inclusive
/// - [`QueryProfileError::BadGapWeights`] if `gap_extend` is less than
///   `gap_open`
#[inline]
pub(crate) fn validate_profile_args<Q: AsRef<[u8]>>(
    query: Q, gap_open: i8, gap_extend: i8,
) -> Result<(), QueryProfileError> {
    if query.as_ref().is_empty() {
        Err(QueryProfileError::EmptyQuery)
    } else if !(-127..=0).contains(&gap_open) {
        Err(QueryProfileError::GapOpenOutOfRange { gap_open })
    } else if !(-127..=0).contains(&gap_extend) {
        Err(QueryProfileError::GapExtendOutOfRange { gap_extend })
    } else if gap_extend < gap_open {
        Err(QueryProfileError::BadGapWeights { gap_open, gap_extend })
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
/// ## Type Parameters
///
/// - `S` - The size of the alphabet (usually 5 for DNA including *N*)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScalarProfile<'a, const S: usize> {
    pub(crate) seq:        &'a [u8],
    pub(crate) matrix:     &'a WeightMatrix<'a, i8, S>,
    pub(crate) gap_open:   i32,
    pub(crate) gap_extend: i32,
}

impl<'a, const S: usize> ScalarProfile<'a, S> {
    /// Creates a new profile for use with the scalar alignment algorithm.
    ///
    /// ## Errors
    ///
    /// The following errors are possible:
    ///
    /// - [`QueryProfileError::EmptyQuery`] if `query` is empty
    /// - [`QueryProfileError::GapOpenOutOfRange`] if `gap_open` is not between
    ///   -127 and 0, inclusive
    /// - [`QueryProfileError::GapExtendOutOfRange`] if `gap_extend` is not
    ///   between -127 and 0, inclusive
    /// - [`QueryProfileError::BadGapWeights`] if `gap_extend` is less than
    ///   `gap_open`
    pub fn new<Q: AsRef<[u8]> + ?Sized>(
        query: &'a Q, matrix: &'a WeightMatrix<'a, i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, QueryProfileError> {
        validate_profile_args(query, gap_open, gap_extend)?;
        // Validity: validate_profile_args has checked our assumptions
        Ok(Self::new_unchecked(query, matrix, gap_open, gap_extend))
    }

    /// Creates a new profile for use with the scalar alignment algorithm
    /// without checking for validity of the arguments.
    ///
    /// ## Validity
    ///
    /// The query must be non-empty, `gap_open` and `gap_extend` must be in the
    /// range `-127..=0`, and `gap_extend` must be greater than or equal to
    /// `gap_open` (`gap_open` must be a worse penalty).
    pub fn new_unchecked<Q: AsRef<[u8]> + ?Sized>(
        query: &'a Q, matrix: &'a WeightMatrix<'a, i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Self {
        ScalarProfile {
            seq: query.as_ref(),
            matrix,
            gap_open: i32::from(gap_open),
            gap_extend: i32::from(gap_extend),
        }
    }

    /// Computes the Smith-Waterman local alignment score between the profile
    /// and a passed sequence.
    ///
    /// For more info, see: [`sw_scalar_score`].
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{alignment::{ScalarProfile, sw::sw_scalar_score}, data::matrices::WeightMatrix};
    /// let reference: &[u8] = b"GGCCACAGGATTGAG";
    /// let query: &[u8] = b"CTCAGATTG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = ScalarProfile::<5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.smith_waterman_score(reference).unwrap();
    /// assert_eq!(score, 27);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, reference: &[u8]) -> MaybeAligned<u32> {
        sw_scalar_score(reference, self)
    }

    /// Computes the Smith-Waterman local alignment between the profile and a
    /// passed sequence, yielding the reference starting position (indexed from
    /// 1), cigar, and optimal score.
    ///
    /// For more info, see: [`sw_scalar_alignment`].
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     alignment::{Alignment, ScalarProfile, sw::sw_scalar_score},
    /// #     data::{matrices::WeightMatrix, cigar::Cigar}
    /// # };
    /// let reference: &[u8] = b"GGCCACAGGATTGAG";
    /// let query: &[u8] = b"CTCAGATTG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = ScalarProfile::<5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let alignment = profile.smith_waterman_alignment(reference).unwrap();
    /// assert_eq!(alignment.ref_range.start, 3);
    /// assert_eq!(alignment.states, Cigar::from_slice_unchecked("5M1D4M"));
    /// assert_eq!(alignment.score, 27);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_alignment(&self, reference: &[u8]) -> MaybeAligned<Alignment<u32>> {
        sw_scalar_alignment(reference, self)
    }
}

/// A striped profile for DNA sequence alignment using SIMD operations.
///
/// The profile contains pre-computed scoring vectors for each residue in the
/// specified mapping, arranged in a striped pattern to optimize SIMD operations
/// during alignment.
///
/// ## Type Parameters
///
/// - `T` - The numeric type used for scores. i8, i16, i32, and i64 use the
///   signed algorithm, which is the most common. u8, u16, u32, and u64 use the
///   unsigned algorithm.
/// - `N` - The number of SIMD lanes (usually 16, 32 or 64)
/// - `S` - The size of the alphabet (usually 5 for DNA including *N*)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct StripedProfile<'a, T, const N: usize, const S: usize>
where
    T: SimdElement,
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) profile:    Vec<Simd<T, N>>,
    pub(crate) gap_open:   T,
    pub(crate) gap_extend: T,
    pub(crate) bias:       T,
    pub(crate) mapping:    &'a ByteIndexMap<S>,
    pub(crate) seq_len:    usize,
}

impl<'a, T, const N: usize, const S: usize> StripedProfile<'a, T, N, S>
where
    T: AnyInt + SimdElement,
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates a new striped profile from a sequence and scoring matrix.
    ///
    /// See: [`WeightMatrix`]
    ///
    /// ## Errors
    ///
    /// The following errors are possible:
    ///
    /// - [`QueryProfileError::EmptyQuery`] if `query` is empty
    /// - [`QueryProfileError::GapOpenOutOfRange`] if `gap_open` is not between
    ///   -127 and 0, inclusive
    /// - [`QueryProfileError::GapExtendOutOfRange`] if `gap_extend` is not
    ///   between -127 and 0, inclusive
    /// - [`QueryProfileError::BadGapWeights`] if `gap_extend` is less than
    ///   `gap_open`
    pub fn new<U>(
        query: &[u8], matrix: &WeightMatrix<'a, U, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, QueryProfileError>
    where
        T: FromSameSignedness<U> + AlignableIntWidth,
        U: AnyInt, {
        validate_profile_args(query, gap_open, gap_extend)?;
        // Validity: we validated the required assumptions with
        // validate_profile_args
        Ok(Self::new_unchecked(query, matrix, gap_open, gap_extend))
    }

    /// Creates a new striped profile from a sequence and scoring matrix.
    ///
    /// See: [`WeightMatrix`]
    ///
    /// ## Validity
    ///
    /// The query must be non-empty, `gap_open` and `gap_extend` must be in the
    /// range `-127..=0`, and `gap_extend` must be greater than or equal to
    /// `gap_open` (`gap_open` must be a worse penalty).
    pub fn new_unchecked<U>(query: &[u8], matrix: &WeightMatrix<'a, U, S>, gap_open: i8, gap_extend: i8) -> Self
    where
        T: From<U>,
        U: AnyInt, {
        // SupportedLaneCount cannot presently be zero.
        let number_vectors = query.len().div_ceil(N);
        let total_lanes = N * number_vectors;

        let bias = matrix.bias.into();
        let biases = Simd::splat(bias);
        let mut profile = vec![biases; S * number_vectors];

        for v in 0..number_vectors {
            for ref_index in 0..matrix.mapping.len() {
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
            gap_open: T::from_literal(-gap_open),
            gap_extend: T::from_literal(-gap_extend),
            bias,
            mapping: matrix.mapping,
            seq_len: query.len(),
        }
    }

    /// Given a profile for a forward sequence, returns a reversed profile of a
    /// subsequence ending at the given index.
    ///
    /// `seq_end` represents the 0-based exclusive-end index. `None` is returned
    /// if the bound is invalid.
    #[must_use]
    #[cfg(feature = "dev-3pass")]
    pub fn reverse_from_forward(&self, seq_end: usize) -> Option<Self> {
        if seq_end == 0 || seq_end > self.seq_len {
            return None;
        }

        let number_vectors = seq_end.div_ceil(N);
        let number_vecs_old = self.number_vectors();
        let total_lanes = N * number_vectors;

        let biases = Simd::splat(self.bias);
        let mut profile = vec![biases; S * number_vectors];

        for v in 0..number_vectors {
            for ref_index in 0..self.mapping.len() {
                let mut vector = biases;
                for (i, q) in (v..total_lanes).step_by(number_vectors).enumerate() {
                    if q < seq_end {
                        let q_old = seq_end - 1 - q;
                        let v_old = q_old % number_vecs_old;
                        let lane_old = q_old / number_vecs_old;

                        vector[i] = self.profile[ref_index * number_vecs_old + v_old][lane_old];
                    }
                }
                profile[ref_index * number_vectors + v] = vector;
            }
        }

        Some(StripedProfile {
            profile,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            bias: self.bias,
            mapping: self.mapping,
            seq_len: seq_end,
        })
    }

    /// Returns a sub-profile for a specific range of the query sequence.
    ///
    /// `None` is returned if the bound is invalid.
    #[must_use]
    #[cfg(feature = "dev-3pass")]
    pub fn new_with_range(&self, range: std::ops::Range<usize>) -> Option<Self> {
        if range.is_empty() || range.end > self.seq_len {
            return None;
        }

        let new_len = range.end - range.start;
        let number_vectors = new_len.div_ceil(N);
        let number_vecs_old = self.number_vectors();
        let total_lanes = N * number_vectors;

        let biases = Simd::splat(self.bias);
        let mut profile = vec![biases; S * number_vectors];

        for v in 0..number_vectors {
            for ref_index in 0..self.mapping.len() {
                let mut vector = biases;
                for (i, q) in (v..total_lanes).step_by(number_vectors).enumerate() {
                    if q < new_len {
                        let q_old = range.start + q;
                        let v_old = q_old % number_vecs_old;
                        let lane_old = q_old / number_vecs_old;

                        vector[i] = self.profile[ref_index * number_vecs_old + v_old][lane_old];
                    }
                }
                profile[ref_index * number_vectors + v] = vector;
            }
        }

        Some(StripedProfile {
            profile,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            bias: self.bias,
            mapping: self.mapping,
            seq_len: new_len,
        })
    }

    /// Returns the number of SIMD vectors in the profile.
    #[inline]
    #[must_use]
    pub fn number_vectors(&self) -> usize {
        self.profile.len() / S
    }
}

impl<T, const N: usize, const S: usize> StripedProfile<'_, T, N, S>
where
    LaneCount<N>: SupportedLaneCount,
    T: AlignableIntWidth,
    Simd<T, N>: SimdAnyInt<T, N>,
{
    /// Computes the Smith-Waterman local alignment score between the query
    /// profile and the passed reference.
    ///
    /// Returns [`Overflowed`] if the score overflowed and [`Unmapped`] if no
    /// portion of the query mapped to the reference.
    ///
    /// For more info, see: [`sw_simd_score`].
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{alignment::{StripedProfile, sw::sw_simd_score}, data::matrices::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = StripedProfile::<u8, 32, 5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.smith_waterman_score(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    ///
    /// [`Overflowed`]: MaybeAligned::Overflowed
    /// [`Unmapped`]: MaybeAligned::Unmapped
    #[inline]
    #[must_use]
    pub fn smith_waterman_score(&self, reference: &[u8]) -> MaybeAligned<u32> {
        sw_simd_score::<T, N, S>(reference, self)
    }

    /// Similar to [`smith_waterman_score`] but includes reference and query
    /// 0-based, exclusive end indices.
    ///
    /// Note: these coordinates are equivalent to the 1-based end positions.
    ///
    /// [`smith_waterman_score`]: Self::smith_waterman_score
    #[inline]
    #[must_use]
    pub fn smith_waterman_score_ends(&self, reference: &[u8]) -> MaybeAligned<(u32, usize, usize)> {
        sw_simd_score_ends::<T, N, S>(reference, self)
    }

    /// Similar to [`sw_simd_score`], but includes the reference and query
    /// inclusive start index (0-based).
    ///
    /// <div class="warning important">
    ///
    /// **Important**
    ///
    /// The query profile should also be in reverse orientation, such as using
    /// [`reverse_from_forward`].
    ///
    /// </div>
    ///
    /// [`reverse_from_forward`]: StripedProfile::reverse_from_forward
    #[inline]
    #[must_use]
    #[cfg(feature = "dev-3pass")]
    pub(crate) fn smith_waterman_score_ends_reverse(&self, reference: &[u8]) -> MaybeAligned<(u32, usize, usize)> {
        crate::alignment::sw::sw_simd_score_ends_reverse::<T, N, S>(reference, self)
    }

    /// Computes the Smith-Waterman local alignment between the query profile
    /// and the passed reference.
    ///
    /// Returns [`Overflowed`] if the score overflowed and [`Unmapped`] if no
    /// portion of the query mapped to the reference.
    ///
    /// For more info, see: [`sw_simd_alignment`].
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     alignment::{StripedProfile, sw::sw_simd_score},
    /// #     data::matrices::WeightMatrix
    /// # };
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = StripedProfile::<u8, 32, 5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let alignment = profile.smith_waterman_alignment(reference).unwrap();
    /// // alignment contains ref_range, states (cigar), and score
    /// ```
    ///
    /// [`Overflowed`]: MaybeAligned::Overflowed
    /// [`Unmapped`]: MaybeAligned::Unmapped
    #[inline]
    #[must_use]
    pub fn smith_waterman_alignment(&self, reference: &[u8]) -> MaybeAligned<Alignment<u32>> {
        sw_simd_alignment::<T, N, S>(reference, self)
    }

    /// Similar to [`smith_waterman_score`] but includes the reference and query
    /// alignment ranges.
    ///
    /// These are standard 0-based, half-open ranges for slicing.
    ///
    /// [`smith_waterman_score`]: Self::smith_waterman_score
    #[inline]
    #[must_use]
    #[cfg(feature = "dev-3pass")]
    pub fn smith_waterman_score_ranges(
        &self, reference: &[u8],
    ) -> MaybeAligned<(u32, std::ops::Range<usize>, std::ops::Range<usize>)> {
        sw_simd_score_ranges::<T, N, S>(reference, self)
    }

    /// Similar to [`smith_waterman_alignment`] but computes the Smith-Waterman
    /// local alignment using a 3-pass algorithm.
    ///
    /// This approach was inspired by (7), albeit this implementation is for
    /// large reference / small query because banded is not used.
    ///
    /// See **[module citations](crate::alignment::sw#module-citations)**.
    ///
    /// [`smith_waterman_alignment`]: Self::smith_waterman_alignment
    #[inline]
    #[must_use]
    #[cfg(feature = "dev-3pass")]
    #[allow(clippy::missing_panics_doc)]
    pub fn smith_waterman_alignment_3pass(
        &self, reference: &[u8], original_query: &[u8], matrix: &WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> MaybeAligned<Alignment<u32>> {
        sw_simd_score_ranges::<T, N, S>(reference, self).and_then(|(score, ref_range, query_range)| {
            use crate::data::views::IndexAdjustable;

            if query_range.is_empty() {
                return MaybeAligned::Unmapped;
            } else if query_range.len() == ref_range.len()
            // TODO: this incantation can likely be done better
                && reference[ref_range.clone()]
                    .iter()
                    .zip(original_query[query_range.clone()].iter())
                    .map(|(&r, &q)| i32::from(matrix.get_weight(r, q)))
                    .sum::<i32>()
                    .try_into()
                    .unwrap_or(0)
                    == score
            {
                let states = super::AlignmentStates::new_no_gaps(query_range.clone(), original_query.len());
                return MaybeAligned::Some(Alignment {
                    score,
                    ref_range,
                    query_range,
                    states,
                    ref_len: reference.len(),
                    query_len: original_query.len(),
                });
            }

            // Validity: query_range is non-empty per check, and weights are
            // same as those in self
            let query_new =
                ScalarProfile::<S>::new_unchecked(&original_query[query_range.clone()], matrix, gap_open, gap_extend);
            let reference_new = &reference[ref_range.clone()];

            let mut band_width = ref_range.len().abs_diff(query_range.len()) + 1;

            let max_bandwidth = (query_range.len() - 1) / 2;
            let mut banded_alignment = None;

            while band_width <= max_bandwidth {
                if let MaybeAligned::Some(alignment) = sw_banded_alignment(reference_new, &query_new, band_width)
                    && alignment.score == score
                {
                    banded_alignment = Some(alignment);
                    break;
                }
                band_width *= 2;
            }

            let mut alignment = match banded_alignment {
                Some(alignment) => alignment,
                None => sw_scalar_alignment(reference_new, &query_new).unwrap(),
            };

            // The alignment found may not span the full bounding box, since
            // multiple alignments of different lengths could have the same
            // score
            let adjusted_query_range = alignment.query_range.add(query_range.start);
            let adjusted_ref_range = alignment.ref_range.add(ref_range.start);

            // Note: banded alignment may not cover the full reference/query ranges like SIMD
            alignment.states.prepend_soft_clip(adjusted_query_range.start);
            alignment.states.soft_clip(original_query.len() - adjusted_query_range.end);

            MaybeAligned::Some(Alignment {
                score,
                ref_range: adjusted_ref_range,
                query_range: adjusted_query_range,
                states: alignment.states,
                ref_len: reference.len(),
                query_len: original_query.len(),
            })
        })
    }
}

#[cfg(test)]
mod bench {
    use test::Bencher;
    extern crate test;
    use super::*;
    use crate::alignment::sw::test_data::{GAP_EXTEND, GAP_OPEN};
    use crate::data::{constants::mappings::DNA_PROFILE_MAP, matrices::WeightMatrix};
    pub(crate) static DATA: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/KJ907631.1.txt")); // H5 HA, complete CDS
    pub(crate) static MATRIX: WeightMatrix<u8, 5> =
        WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N')).to_biased_matrix();

    #[bench]
    fn build_profile_u8(b: &mut Bencher) {
        b.iter(|| std::hint::black_box(StripedProfile::<u8, 32, 5>::new(DATA, &MATRIX, GAP_OPEN, GAP_EXTEND)));
    }

    #[bench]
    fn build_profile_u8_rev_half_recreate(b: &mut Bencher) {
        b.iter(|| {
            let data_rev: Vec<u8> = DATA[..DATA.len() / 2].iter().copied().rev().collect();
            std::hint::black_box(StripedProfile::<u8, 32, 5>::new(&data_rev, &MATRIX, GAP_OPEN, GAP_EXTEND))
        });
    }

    #[cfg(feature = "dev-3pass")]
    #[bench]
    fn build_profile_u8_rev_half_mapping(b: &mut Bencher) {
        let prof = StripedProfile::<u8, 32, 5>::new(DATA, &MATRIX, GAP_OPEN, GAP_EXTEND).unwrap();

        b.iter(|| std::hint::black_box(prof.reverse_from_forward(DATA.len() / 2)));
    }
}
