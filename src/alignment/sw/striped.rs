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
pub fn sw_simd_score<T, const N: usize, const S: usize>(
    reference: &[u8], query: &StripedProfile<T, N, S>,
) -> MaybeAligned<u32>
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
    let mut max_scores = minimums;

    for ref_index in reference.iter().copied().map(|r| query.mapping.to_index(r)) {
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
    score_to_maybe_aligned(best, query.bias, |score| score)
}

/// Similar to [`sw_simd_score`] but includes reference and query 0-based,
/// exclusive end indices.
///
/// Note: these coordinates are equivalent to the 1-based end positions.
#[inline]
#[must_use]
pub fn sw_simd_score_ends<T, const N: usize, const S: usize>(
    reference: &[u8], query: &StripedProfile<T, N, S>,
) -> MaybeAligned<(u32, usize, usize)>
where
    T: AlignableIntWidth,
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdAnyInt<T, N>, {
    sw_simd_score_ends_dir::<T, N, S, true>(reference, query)
}

/// Similar to [`sw_simd_score`], but includes the reference and query
/// inclusive start index (0-based).
///
/// Important: the query profile should also be in reverse orientation.
#[inline]
#[must_use]
#[cfg(feature = "dev-3pass")]
pub(crate) fn sw_simd_score_ends_reverse<T, const N: usize, const S: usize>(
    reference: &[u8], query: &StripedProfile<T, N, S>,
) -> MaybeAligned<(u32, usize, usize)>
where
    T: AlignableIntWidth,
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdAnyInt<T, N>, {
    sw_simd_score_ends_dir::<T, N, S, false>(reference, query)
}

/// Similar to [`sw_simd_score`] but also returns the reference and query end or
/// start coordinates.
///
/// In the forward direction, the end coordinate is equivalently:
/// - The 1-based position of the last aligning character
/// - The 0-based exclusive end for the alignment range
///
/// For reverse, the start coordinate is:
/// - The 0-based inclusive start for the alignment range
#[allow(non_snake_case)]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
fn sw_simd_score_ends_dir<T, const N: usize, const S: usize, const FORWARD: bool>(
    reference: &[u8], query: &StripedProfile<T, N, S>,
) -> MaybeAligned<(u32, usize, usize)>
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

    let saturating_threshold = if T::SIGNED { T::MAX } else { T::MAX - query.bias };

    let mut load = vec![minimums; num_vecs];
    let mut store = vec![minimums; num_vecs];
    let mut e_scores = vec![minimums; num_vecs];
    let mut max_row = vec![minimums; num_vecs];

    let mut best = min;
    let mut r_end = reference.len() - 1;

    let len = reference.len();
    for r in 0..len {
        let ref_index = query.mapping.to_index(reference[if FORWARD { r } else { len - 1 - r }]);
        let mut F = minimums;
        let mut H = store[num_vecs - 1].shift_elements_right::<1>(min);

        if r > 1 && r_end == r - 2 {
            (max_row, load) = (load, max_row);
        }
        (load, store) = (store, load);

        // This statement helps with bounds checks.
        let scores_vec = &profile[(ref_index * num_vecs)..(ref_index * num_vecs + num_vecs)];
        let mut max_scores = minimums;

        for v in 0..num_vecs {
            let mut E = e_scores[v];

            H = H.saturating_add(scores_vec[v]);
            if !T::SIGNED {
                H = H.saturating_sub(biases);
            }

            H = H.simd_max(E).simd_max(F);

            max_scores = max_scores.simd_max(H);

            store[v] = H;

            H = H.saturating_sub(gap_opens);
            E = E.saturating_sub(gap_extends).simd_max(H);
            F = F.saturating_sub(gap_extends).simd_max(H);

            e_scores[v] = E;
            H = load[v];
        }

        'lazy_f: for _ in 0..N {
            F = F.shift_elements_right::<1>(min);

            for store_v in &mut store {
                H = *store_v;

                if !F.simd_gt(H.saturating_sub(gap_opens)).any() {
                    break 'lazy_f;
                }

                H = H.simd_max(F);
                *store_v = H;

                F = F.saturating_sub(gap_extends);
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
        let lane = ci / num_vecs;

        if max_row[v][lane] == best {
            c_end = ci;
            break;
        }
    }

    if FORWARD {
        r_end += 1;
        c_end += 1;
    } else {
        r_end = reference.len() - 1 - r_end;
        c_end = query.seq_len - 1 - c_end;
    }

    score_to_maybe_aligned(best, query.bias, |score| (score, r_end, c_end))
}

/// Similar to [`sw_simd_score`] but also returns the reference and query
/// alignment ranges for 0-based slicing. The algorithm performs a truncated two
/// pass approach. This approach was inspired by (7).
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
#[cfg(feature = "dev-3pass")]
#[inline]
#[must_use]
pub fn sw_simd_score_ranges<T, const N: usize, const S: usize>(
    reference: &[u8], query: &StripedProfile<T, N, S>,
) -> MaybeAligned<(u32, Range<usize>, Range<usize>)>
where
    T: AlignableIntWidth,
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdAnyInt<T, N>, {
    sw_simd_score_ends::<T, N, S>(reference, query).and_then(|(score, ref_end, query_end)| {
        let Some(query_rev) = query.reverse_from_forward(query_end) else {
            return MaybeAligned::Unmapped;
        };
        sw_simd_score_ends_reverse::<T, N, S>(&reference[..ref_end], &query_rev).map(|(score2, ref_start, query_start)| {
            debug_assert_eq!(score, score2);
            (score, ref_start..ref_end, query_start..query_end)
        })
    })
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
) -> MaybeAligned<Alignment<u32>>
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

    for (r, ref_index) in reference.iter().copied().map(|r| query.mapping.to_index(r)).enumerate() {
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
        // Also equal to (ci - v) / num_vecs. Subtracting the division remainder
        // before dividing is equivalent to integer (truncating) division.
        let lane = ci / num_vecs;

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

    score_to_maybe_aligned(best, query.bias, |score| {
        backtrack.to_alignment(score, r_end, c_end, reference.len(), query.seq_len)
    })
}

/// For test purposes only
/// This approach was inspired by (7), albeit this implementation is for large
/// reference / small query because banded is not used.
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
#[must_use]
#[cfg(feature = "dev-3pass")]
pub fn sw_simd_alignment_3pass<T, const N: usize, const S: usize>(
    reference: &[u8], query: &StripedProfile<T, N, S>,
) -> MaybeAligned<Alignment<u32>>
where
    T: AlignableIntWidth,
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdAnyInt<T, N>, {
    sw_simd_score_ranges::<T, N, S>(reference, query).and_then(|(score, ref_range, query_range)| {
        let Some(query_new) = query.new_with_range(query_range.clone()) else {
            return MaybeAligned::Unmapped;
        };
        let reference_new = &reference[ref_range.clone()];
        sw_simd_alignment(reference_new, &query_new).map(|mut alignment| {
            debug_assert_eq!(alignment.score, score);
            debug_assert_eq!(alignment.ref_range, 0..reference_new.len());
            debug_assert_eq!(alignment.query_range, 0..query_new.seq_len);
            alignment.states.prepend_soft_clip(query_range.start);
            alignment.states.soft_clip(query.seq_len - query_range.end);
            Alignment {
                score,
                ref_range,
                query_range,
                states: alignment.states,
                ref_len: reference.len(),
                query_len: query.seq_len,
            }
        })
    })
}

/// Converts a the `best` score seen so far by a Striped Smith Waterman
/// algorithm to [`MaybeAligned`].
///
/// Given `best` and the corresponding `bias` used by the algorithm, convert
/// `best` to an unbiased `u32` score, then potentially return
/// [`MaybeAligned::Overflowed`] or [`MaybeAligned::Unmapped`]. If neither of
/// those are returned, call `f` to convert the `u32` score into the desired
/// output type `R`.
#[inline]
#[must_use]
fn score_to_maybe_aligned<T, F, R>(best: T, bias: T, f: F) -> MaybeAligned<R>
where
    T: AlignableIntWidth,
    F: FnOnce(u32) -> R, {
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
        best.checked_add(bias + T::ONE).map(|_| best.cast_as::<u32>())
    }
    .map_or(MaybeAligned::Overflowed, |score| {
        if score == 0 {
            MaybeAligned::Unmapped
        } else {
            MaybeAligned::Some(f(score))
        }
    })
}
