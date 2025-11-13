// TODO: revisit truncation issues
use super::*;
use std::ops::Add;

/// Smith-Waterman algorithm (non-vectorized), yielding the optimal score.
///
/// Provides the locally optimal sequence alignment score (1) using affine gap
/// penalties (2). Our implementation adapts the algorithm provided by [*Flouri
/// et al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html)
/// (4).
///
/// This implementation will not return [`MaybeAligned::Overflowed`].
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
///
/// ## Complexity
///
/// For query length $m$ and reference length $n$:
///
/// - Time: $O(mn)$
/// - Space: $O(m)$
///
/// ## Limitations
///
/// - Not suitable for production due to high runtimes.
/// - Useful for testing correctness.
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
/// let score = sw_scalar_score(&reference, &profile).unwrap();
/// assert_eq!(score, 27);
/// ```
#[must_use]
#[allow(clippy::cast_sign_loss)]
pub fn sw_scalar_score<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> MaybeAligned<u32> {
    // Definitions:
    // * Q[c] = c-th base in the query, zero-indexed
    // * R[r] = r-th base in the reference, zero-indexed
    // * W[r,c] = contribution to score of aligning Q[c] and R[r]
    // * H[r,c] = maximum score for aligning Q[0..c] to R[0..r]
    // * E[r,c] = maximum score for aligning Q[0..c] to R[0..r] such that the
    //   alignment ends with a gap consuming a base in the reference
    // * F[r,c] = maximum score for aligning Q[0..c] to R[0..r] such that the
    //   alignment ends with a gap consuming a base in the query

    // The maximum value of H[r,c] for all r and c, to be updated during the
    // algorithm. Initialized to the minimum of 0
    let mut best_score = 0;
    // A row H[r, ..] of H. Initialized for row -1 (no bases in the reference
    // consumed)
    let mut h_row = vec![0; query.seq.len()];
    // A row E[r, ..] of E. Initialized for row 0 (one base in the reference
    // consumed, via a gap)
    let mut e_row = vec![query.gap_open; query.seq.len()];

    // Iterate over R (index=r)
    for reference_base in reference.iter().copied() {
        // A value in F, initialized to F[0,0]. F[0,0] is reached by excluding
        // R[0] from the alignment, then consuming Q[0] via a gap
        let mut f = query.gap_open;
        // A value of H. Initialized to H[r-1,-1] (r bases consumed in R and no
        // bases consumed in Q)
        let mut h = 0;

        // Iterate over Q (index=c). Assumes that h contains H[r-1,c-1].
        for c in 0..query.seq.len() {
            // W[r,c]
            let match_score = i32::from(query.matrix.get_weight(reference_base, query.seq[c]));
            // Temporarily update to H[r-1,c-1] + W[r,c]
            h += match_score;

            // Extract E[r,c]
            let mut e = e_row[c];
            // Update H[r-1,c-1] to H[r,c] by computing the maximum of
            // H[r-1,c-1] + W[r,c], E[r,c], F[r,c], and 0
            h = h.max(e).max(f).max(0);

            // Update best_score with the new H[r,c]
            best_score = best_score.max(h);

            // Update E[r,c] to E[r,c+1]
            e = e.add(query.gap_extend).max(h + query.gap_open);
            // Update F[r,c] to F[r,c+1]
            f = f.add(query.gap_extend).max(h + query.gap_open);

            // Update H[r,c] in h_row, and set h to H[r-1,c] (becomes H[r-1,c-1]
            // next iteration)
            (h, h_row[c]) = (h_row[c], h);
            // Update E[r,c] in e_row
            e_row[c] = e;
        }
    }

    // score is i32 and is non-negative due to this algorithm, so this is a
    // valid cast
    let best_score = best_score as u32;
    if best_score > 0 {
        MaybeAligned::Some(best_score)
    } else {
        MaybeAligned::Unmapped
    }
}

/// Smith-Waterman algorithm (non-vectorized), yielding the optimal alignment.
///
/// Provides the locally optimal sequence alignment (1) using affine gap
/// penalties (2) with improvements and corrections by (3-4). Our implementation
/// adapts the corrected algorithm provided by [*Flouri et
/// al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html)(4).
///
/// This implementation will not return [`MaybeAligned::Overflowed`].
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
///
/// ## Complexity
///
/// For query length $m$ and reference length $n$:
///
/// - Time: $O(mn)$
/// - Space: $O(mn)$
///
/// ## Limitations
///
/// - Not suitable for large sequences due to high memory usage.
/// - Not suitable for production due to high runtimes.
/// - Useful for testing correctness.
///
/// ## Example
///
///  ```
/// # use zoe::{
/// #     alignment::{Alignment, ScalarProfile, sw::sw_scalar_alignment},
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
/// let alignment = sw_scalar_alignment(&reference, &profile).unwrap();
/// assert_eq!(alignment.ref_range.start, 3);
/// assert_eq!(alignment.states, Cigar::from_slice_unchecked("5M1D4M"));
/// assert_eq!(alignment.score, 27);
/// ```
#[must_use]
#[allow(clippy::cast_sign_loss)]
pub fn sw_scalar_alignment<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> MaybeAligned<Alignment<u32>> {
    // See dev comments in sw_scalar_score for more details

    // TODO: Potentially remove this, since it isn't necessary and isn't an
    // expected case
    if reference.is_empty() {
        return MaybeAligned::Unmapped;
    }

    let mut best_score = 0;
    // The value of r and c yielding the highest H[r,c] in best_score. The
    // initialization is irrelevant, since they will be overwritten when
    // best_score is updated, or they won't be accessed (since best_score==0
    // returns Unmapped)
    let (mut r_end, mut c_end) = (0, 0);
    let mut h_row = vec![0; query.seq.len()];
    let mut e_row = vec![query.gap_open; query.seq.len()];

    let mut backtrack = BacktrackMatrix::new(reference.len(), query.seq.len());

    for (r, reference_base) in reference.iter().copied().enumerate() {
        let mut f = query.gap_open;
        let mut h = 0;

        for c in 0..query.seq.len() {
            backtrack.move_to(r, c);

            let match_score = i32::from(query.matrix.get_weight(reference_base, query.seq[c]));
            h += match_score;

            let mut e = e_row[c];
            h = h.max(e).max(f).max(0);

            // Update best_score with the new H[r,c], along with r_end and c_end
            if h > best_score {
                best_score = h;
                r_end = r;
                c_end = c;
            }

            if e == h {
                backtrack.up();
            }

            if f == h {
                backtrack.left();
            }

            if h == 0 {
                backtrack.stop();
            }

            // H[r-1,c]
            let next_diag = h_row[c];
            // Update H[r,c] in h_row
            h_row[c] = h;

            // Temporarily store H[r,c]+go in h
            h += query.gap_open;
            // Update E[r,c] to E[r,c+1] by computing the maximum of E[r,c]+ge
            // and H[r,c]+go
            e = e.add(query.gap_extend).max(h);
            // Update F[r,c] to F[r,c+1] by computing the maximum of F[r,c]+ge
            // and F[r,c]+go
            f = f.add(query.gap_extend).max(h);

            // Check H[r,c]≠0
            if h != query.gap_open {
                // Mark the backtrack as up_extending only if it cannot be
                // reached another way. It cannot be reached via clipping since
                // H[r,c]≠0. This check ensures E[r,c+1] > H[r,c]+go, so that
                // E[r,c+1] = E[r,c]+ge.
                if e > h {
                    backtrack.up_extending();
                }

                // Mark the backtrack as left_extending only if it cannot be
                // reached another way. It cannot be reached via clipping since
                // H[r,c]≠0. This check ensures F[r,c+1] > H[r,c]+go, so that
                // F[r,c+1] = F[r,c]+ge.
                if f > h {
                    backtrack.left_extending();
                }
            }

            // Store H[r-1,c] in h (becomes H[r-1,c-1] next iteration)
            h = next_diag;
            e_row[c] = e;
        }
    }

    if best_score == 0 {
        MaybeAligned::Unmapped
    } else {
        // score is i32 and is non-negative due to this algorithm, so this is a
        // valid cast
        MaybeAligned::Some(backtrack.to_alignment(best_score as u32, r_end, c_end, reference.len(), query.seq.len()))
    }
}

/// Similar to [`sw_scalar_alignment`], but allows the user to pass a closure to
/// selectively alter certain scores in the DP table.
///
/// The closure `alter_score` accepts the row index, the column index, and the
/// old score as arguments and should yield the new score.
#[must_use]
#[cfg(feature = "alignment-diagnostics")]
#[allow(clippy::cast_sign_loss)]
pub fn sw_scalar_alignment_override<F, const S: usize>(
    reference: &[u8], query: &ScalarProfile<S>, mut alter_score: F,
) -> MaybeAligned<Alignment<u32>>
where
    F: FnMut(usize, usize, i32) -> i32, {
    if reference.is_empty() {
        return MaybeAligned::Unmapped;
    }

    let (mut best_score, mut r_end, mut c_end) = (0, 0, 0);
    let mut h_row = vec![0; query.seq.len()];
    let mut e_row = vec![query.gap_open; query.seq.len()];

    let mut backtrack = BacktrackMatrix::new(reference.len(), query.seq.len());

    for (r, reference_base) in reference.iter().copied().enumerate() {
        let mut f = query.gap_open;
        let mut h = 0;

        for c in 0..query.seq.len() {
            backtrack.move_to(r, c);

            // matching is the default direction
            let match_score = i32::from(query.matrix.get_weight(reference_base, query.seq[c]));
            h += match_score;

            let mut e = e_row[c];
            h = h.max(e).max(f).max(0);
            h = alter_score(r, c, h);

            if h > best_score {
                best_score = h;
                r_end = r;
                c_end = c;
            }

            if e == h {
                backtrack.up();
            }

            if f == h {
                backtrack.left();
            }

            if h == 0 {
                backtrack.stop();
            }

            let next_diag = h_row[c];
            h_row[c] = h;

            h += query.gap_open;
            e = e.add(query.gap_extend).max(h);
            f = f.add(query.gap_extend).max(h);

            if h != query.gap_open {
                if e > h {
                    backtrack.up_extending();
                }

                if f > h {
                    backtrack.left_extending();
                }
            }

            h = next_diag;
            e_row[c] = e;
        }
    }

    if best_score == 0 {
        MaybeAligned::Unmapped
    } else {
        // score is i32 and is non-negative due to this algorithm, so this is a
        // valid cast
        MaybeAligned::Some(backtrack.to_alignment(best_score as u32, r_end, c_end, reference.len(), query.seq.len()))
    }
}
