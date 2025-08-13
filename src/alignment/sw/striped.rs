use super::*;
use crate::{math::AlignableIntWidth, simd::SimdAnyInt};
use std::simd::{LaneCount, SupportedLaneCount};

/// Smith-Waterman algorithm (vectorized), yielding the optimal score.
///
/// Provides the locally optimal sequence alignment score (1) using affine gap
/// penalties (2). We adapt Farrar's striped SIMD implemention (5) for portable
/// SIMD.
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
///
/// ## Complexity
///
/// For query length $m$, reference length $n$, and $N$ SIMD lanes:
///
/// - Time: $O(mn)$, with the average case as $O(mn/N)$
/// - Space: $O(n)$
///
/// ## Limitations
///
/// - The SIMD algorithm may not perform as well when both the query and
///   reference are short. If the query and reference are both shorter than 25
///   bases, ensure that this algorithm is called with `N=16` or less.
/// - If the query and reference are both shorter than 10 bases, consider using
///   the scalar algorithm [`sw_scalar_score`].
/// - For general use the `multiversion` feature is recommended.
///
/// ## Example
///
/// ```
/// # use zoe::{alignment::{StripedProfile, sw::sw_simd_score}, data::WeightMatrix};
/// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
/// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
///
/// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
/// const GAP_OPEN: i8 = -3;
/// const GAP_EXTEND: i8 = -1;
///
/// let profile = StripedProfile::<u8, 32, 5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
/// let score = sw_simd_score(&reference, &profile).unwrap();
/// assert_eq!(score, 26);
/// ```
#[allow(non_snake_case)]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn sw_simd_score<T, const N: usize, const S: usize>(reference: &[u8], query: &StripedProfile<T, N, S>) -> Option<u32>
where
    T: AlignableIntWidth,
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdAnyInt<T, N>, {
    let num_vecs = query.number_vectors();
    let profile: &Vec<Simd<T, N>> = &query.profile;

    let min = T::MIN;
    let gap_opens = Simd::splat(query.gap_open);
    let gap_extends = Simd::splat(query.gap_extend);
    let minimums = Simd::splat(T::MIN);
    let biases = Simd::splat(query.bias);

    let mut load = vec![minimums; num_vecs];
    let mut store = vec![minimums; num_vecs];
    let mut e_scores = vec![minimums; num_vecs];
    let mut max_scores = minimums; // Minimum value for unsigned

    let map = query.mapping;

    for ref_index in reference.iter().copied().map(|r| map.to_index(r)) {
        let mut F = minimums;
        let mut H = store[num_vecs - 1].shift_elements_right::<1>(min);

        (load, store) = (store, load);

        // This statement helps with bounds checks.
        let scores_vec = &profile[(ref_index * num_vecs)..(ref_index * num_vecs + num_vecs)];

        for j in 0..num_vecs {
            let mut E = e_scores[j];

            H = H.saturating_add(scores_vec[j]);
            if !T::SIGNED {
                H = H.saturating_sub(biases);
            }

            H = H.simd_max(E).simd_max(F);
            max_scores = max_scores.simd_max(H);

            store[j] = H;

            H = H.saturating_sub(gap_opens);
            E = E.saturating_sub(gap_extends).simd_max(H);
            F = F.saturating_sub(gap_extends).simd_max(H);

            // Store E; Load H
            e_scores[j] = E;
            H = load[j];
        }

        let mut j = 0;
        H = store[j];
        F = F.shift_elements_right::<1>(min);

        // ¬∀x (F = (H - Go))
        //  ∃x (F > (H - Go)), given 1-sided
        let mut mask = F.simd_gt(H.saturating_sub(gap_opens));
        while mask.any() {
            H = H.simd_max(F);
            store[j] = H;

            F = F.saturating_sub(gap_extends);

            j += 1;
            if j >= num_vecs {
                j = 0;
                F = F.shift_elements_right::<1>(min);
            }

            H = store[j];

            mask = F.simd_gt(H.saturating_sub(gap_opens));
        }
    }

    let best = max_scores.reduce_max();

    if T::SIGNED {
        // Map best score to an unsigned range. Note that: MAX+1 = abs(MIN). If
        // we would have overflowed (saturated), return None, otherwise return
        // the best score. T is at most i32, so T::MAX + 1 fits in u32. Since
        // best < i32::MAX, the wrapping add will never wrap
        (best < T::MAX).then(|| (T::MAX.cast_as::<u32>() + 1).wrapping_add_signed(best.cast_as::<i32>()))
    } else {
        // If we would have overflowed, return None, otherwise return the best
        // score. We add one because we care if the value is equal to the MAX.
        best.checked_add(query.bias + T::ONE).map(|_| best.cast_as::<u32>())
    }
}

/// Smith-Waterman algorithm (vectorized), yielding the optimal alignment.
///
/// Provides the locally optimal sequence alignment (1) using affine gap
/// penalties (2). We adapt Farrar's striped SIMD implemention (5) for portable
/// SIMD. We take inspiration from previous implementations for the F-loop (6),
/// span calculation (7), and traceback/other optimizations (8).
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
///
/// ## Complexity
///
/// For query length $m$, reference length $n$, and $N$ SIMD lanes:
///
/// - Time: $O(mn)$, with the average case as $O(mn/N)$
/// - Space: $O(mn)$
///
/// ## Limitations
///
/// - This algorithm may not be suitable for large sequence pairs due to high
///   memory usage.
/// - For general use the `multiversion` feature is recommended.
///
/// ## Example
///
/// ```
/// # use zoe::{
/// #     alignment::{Alignment, AlignmentStates, StripedProfile, sw::sw_simd_alignment},
/// #     data::WeightMatrix,
/// # };
///
/// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
/// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
/// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
/// const GAP_OPEN: i8 = -3;
/// const GAP_EXTEND: i8 = -1;
/// let profile = StripedProfile::<u8, 8, 5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
/// let alignment = sw_simd_alignment(reference, &profile).unwrap();
///
/// let Alignment {
///     score,
///     ref_range,
///     query_range,
///     states,
///     ..
/// } = alignment;
/// assert_eq!(states, AlignmentStates::try_from("6M2D9M3S").unwrap());
/// assert_eq!(score, 26);
/// ```
#[allow(non_snake_case, clippy::too_many_lines)]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn sw_simd_alignment<T, const N: usize, const S: usize>(
    reference: &[u8], query: &StripedProfile<T, N, S>,
) -> MaybeAligned<u32>
where
    T: AlignableIntWidth,
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdAnyInt<T, N>, {
    if reference.is_empty() {
        return MaybeAligned::Unmapped;
    }

    let num_vecs = query.number_vectors();
    let profile: &Vec<Simd<T, N>> = &query.profile;

    let min = T::MIN;
    let gap_opens = Simd::splat(query.gap_open);
    let gap_extends = Simd::splat(query.gap_extend);
    let minimums = Simd::splat(T::MIN);
    let biases = Simd::splat(query.bias);

    // Any value strictly less than this threshold did not saturate
    let saturating_threshold = if T::SIGNED { T::MAX } else { T::MAX - query.bias };

    let mut load = vec![minimums; num_vecs];
    let mut store = vec![minimums; num_vecs];
    let mut e_scores = vec![minimums; num_vecs];
    let mut max_row = vec![minimums; num_vecs];

    let mut best = min;
    let mut r_end = reference.len() - 1;

    let mut backtrack = BacktrackMatrixStriped::make_uninit_data(reference.len() * num_vecs);
    let map = query.mapping;

    for (r, ref_index) in reference.iter().copied().map(|r| map.to_index(r)).enumerate() {
        let mut F = minimums;
        let mut H = store[num_vecs - 1].shift_elements_right::<1>(min);

        if r > 1 && r_end == r - 2 {
            (max_row, load) = (load, max_row);
        }
        (load, store) = (store, load);

        // This statement helps with bounds checks.
        let scores_vec = &profile[(ref_index * num_vecs)..(ref_index * num_vecs + num_vecs)];
        let backtrack_row = &mut backtrack[(r * num_vecs)..(r * num_vecs + num_vecs)];
        let mut max_scores = minimums;

        for v in 0..num_vecs {
            let mut E = e_scores[v];

            H = H.saturating_add(scores_vec[v]);
            if !T::SIGNED {
                H = H.saturating_sub(biases);
            }

            H = H.simd_max(E).simd_max(F);

            let mut flags = Simd::<u8, N>::simd_match();
            max_scores = max_scores.simd_max(H);

            flags.simd_up(E.simd_eq(H).cast());
            flags.simd_left(F.simd_eq(H).cast());

            let stopped = H.simd_eq(minimums).cast();

            store[v] = H;

            H = H.saturating_sub(gap_opens);
            E = E.saturating_sub(gap_extends).simd_max(H);
            F = F.saturating_sub(gap_extends).simd_max(H);

            flags.simd_up_extending(E.simd_gt(H).cast());
            flags.simd_left_extending(F.simd_gt(H).cast());
            flags.simd_stop(stopped);

            backtrack_row[v].write(flags);
            e_scores[v] = E;
            H = load[v];
        }

        'lazy_f: for _ in 0..N {
            F = F.shift_elements_right::<1>(min);

            for v in 0..num_vecs {
                H = store[v];

                if !F.simd_gt(H.saturating_sub(gap_opens)).any() {
                    break 'lazy_f;
                }

                H = H.simd_max(F);
                store[v] = H;
                let mut flags = unsafe { backtrack_row[v].assume_init() };

                let stopped = H.simd_eq(minimums);

                flags.simd_correct_and_set_left(F.simd_eq(H).cast());

                H = H.saturating_sub(gap_opens);
                F = F.saturating_sub(gap_extends);

                flags.simd_left_extending(F.simd_gt(H).cast());
                flags.simd_stop(stopped.cast());
                backtrack_row[v].write(flags);
            }
        }

        let row_best = max_scores.reduce_max();
        if row_best > best {
            if row_best >= saturating_threshold {
                return MaybeAligned::Overflowed;
            }
            best = row_best;
            r_end = r;
        }
    }

    if r_end == reference.len() - 1 {
        max_row = store;
    } else if r_end == reference.len() - 2 {
        max_row = load;
    }

    let mut c_end = query.seq_len - 1;

    for ci in 0..query.seq_len {
        let v = ci % num_vecs;
        let lane = (ci - v) / num_vecs;

        if max_row[v][lane] == best {
            c_end = ci;
            break;
        }
    }

    // SAFETY: we have initialized all members of the table in the main loop.
    //
    // Also, this should not re-allocate thanks to equivalent size and
    // alignment:
    // https://doc.rust-lang.org/nightly/src/alloc/vec/in_place_collect.rs.html
    let mut backtrack = BacktrackMatrixStriped::new(
        backtrack.into_iter().map(|uninit| unsafe { uninit.assume_init() }).collect(),
        num_vecs,
    );

    if T::SIGNED {
        // Map best score to an unsigned range. Note that: MAX+1 = abs(MIN). If
        // we would have overflowed (saturated), return None, otherwise return
        // the best score. T is at most i32, so T::MAX + 1 fits in u32. Since
        // best < i32::MAX, the wrapping add will never wrap
        (best < T::MAX).then(|| (T::MAX.cast_as::<u32>() + 1).wrapping_add_signed(best.cast_as::<i32>()))
    } else {
        // If we would have overflowed, return None, otherwise return the best
        // score. We add one because we care if the value is equal to the MAX.
        // TODO: Justify safety
        best.checked_add(query.bias + T::ONE).map(|_| best.cast_as::<u32>())
    }
    .map_or(MaybeAligned::Overflowed, |score| {
        if score == 0 {
            MaybeAligned::Unmapped
        } else {
            backtrack.to_alignment(score, r_end, c_end, reference.len(), query.seq_len)
        }
    })
}
