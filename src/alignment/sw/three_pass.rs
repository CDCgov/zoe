use crate::{
    alignment::{
        Alignment, AlignmentStates, MaybeAligned, ScalarProfile, ScoreAndRanges, StripedProfile,
        sw::{sw_banded_align, sw_scalar_align, sw_simd_score_ranges},
    },
    data::{WeightMatrix, views::IndexAdjustable},
    math::{AlignableIntWidth, SimdAnyInt},
};
use std::simd::Simd;

/// Similar to [`sw_simd_align`] but computes the Smith-Waterman local
/// alignment using a 3-pass algorithm.
///
/// This approach was inspired by (7).
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
///
/// [`sw_simd_align`]: crate::alignment::sw::sw_simd_align
#[inline]
#[must_use]
pub fn sw_align_3pass<T, const N: usize, const S: usize>(
    reference: &[u8], query_profile: &StripedProfile<T, N, S>, query: &[u8], matrix: &WeightMatrix<i8, S>, gap_open: i8,
    gap_extend: i8,
) -> MaybeAligned<Alignment<u32>>
where
    T: AlignableIntWidth,
    Simd<T, N>: SimdAnyInt<T, N>, {
    sw_simd_score_ranges(reference, query_profile).and_then(|score_and_ranges| {
        let ScoreAndRanges {
            score,
            ref_range,
            query_range,
        } = score_and_ranges;

        if query_range.is_empty() {
            return MaybeAligned::Unmapped;
        } else if query_range.len() == ref_range.len()
            // TODO: this incantation can likely be done better
                && reference[ref_range.clone()]
                    .iter()
                    .zip(query[query_range.clone()].iter())
                    .map(|(&r, &q)| i32::from(matrix.get_weight(r, q)))
                    .sum::<i32>()
                    .try_into()
                    .unwrap_or(0)
                    == score
        {
            let states = AlignmentStates::new_no_gaps(query_range.clone(), query.len());
            return MaybeAligned::Some(Alignment {
                score,
                ref_range,
                query_range,
                states,
                ref_len: reference.len(),
                query_len: query.len(),
            });
        }

        // Validity: query_range is non-empty per check, and weights are
        // same as those in self
        let query_new = ScalarProfile::new_unchecked(&query[query_range.clone()], matrix, gap_open, gap_extend);
        let reference_new = &reference[ref_range.clone()];

        let mut band_width = ref_range.len().abs_diff(query_range.len()) + 1;

        let max_bandwidth = (query_range.len() - 1) / 2;
        let mut banded_alignment = None;

        while band_width <= max_bandwidth {
            if let MaybeAligned::Some(alignment) = sw_banded_align(reference_new, &query_new, band_width)
                && alignment.score == score
            {
                banded_alignment = Some(alignment);
                break;
            }
            band_width *= 2;
        }

        let mut alignment = match banded_alignment {
            Some(alignment) => alignment,
            None => sw_scalar_align(reference_new, &query_new).unwrap(),
        };

        // The alignment found may not span the full bounding box, since
        // multiple alignments of different lengths could have the same
        // score
        let adjusted_query_range = alignment.query_range.add(query_range.start);
        let adjusted_ref_range = alignment.ref_range.add(ref_range.start);

        // Note: banded alignment may not cover the full reference/query ranges
        // like SIMD
        alignment.states.prepend_soft_clip(adjusted_query_range.start);
        alignment.states.soft_clip(query.len() - adjusted_query_range.end);

        MaybeAligned::Some(Alignment {
            score,
            ref_range: adjusted_ref_range,
            query_range: adjusted_query_range,
            states: alignment.states,
            ref_len: reference.len(),
            query_len: query.len(),
        })
    })
}
