// TODO: revisit truncation issues
use super::*;
use std::ops::Add;

/// Smith-Waterman algorithm, yielding the optimal score.
///
/// Provides the locally optimal sequence alignment score (1) using affine gap
/// penalties (2). Our implementation adapts the algorithm provided by [*Flouri
/// et al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html)
/// (4).
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
///
/// ## Complexity
///
/// Time: $O(mn)$
///
/// Space: $O(m)$ where $m$ is the length of the query
///
/// ## Limitations
///
/// - Not suitable for production due to high runtimes.
/// - Useful for testing correctness.
///
/// ## Example
///
/// ```
/// # use zoe::{alignment::{ScalarProfile, sw::sw_scalar_score}, data::WeightMatrix};
/// let reference: &[u8] = b"GGCCACAGGATTGAG";
/// let query: &[u8] = b"CTCAGATTG";
///
/// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
/// const GAP_OPEN: i8 = -3;
/// const GAP_EXTEND: i8 = -1;
///
/// let profile = ScalarProfile::<5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
/// let score = sw_scalar_score(&reference, &profile);
/// assert_eq!(score, 27);
/// ```
#[must_use]
#[allow(clippy::cast_sign_loss)]
pub fn sw_scalar_score<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> u64 {
    let mut best_score = 0;
    let mut h_row = vec![0; query.seq.len()];
    let mut e_row: Vec<i32> = vec![query.gap_open; query.seq.len()];

    for reference_base in reference.iter().copied() {
        let mut f = query.gap_open;
        let mut h = 0;

        for c in 0..query.seq.len() {
            let match_score = i32::from(query.matrix.get_weight(reference_base, query.seq[c]));
            h += match_score;

            let mut e = e_row[c];
            h = h.max(e).max(f).max(0);

            best_score = best_score.max(h);

            e = e.add(query.gap_extend).max(h + query.gap_open);
            f = f.add(query.gap_extend).max(h + query.gap_open);

            (h, h_row[c]) = (h_row[c], h);
            e_row[c] = e;
        }
    }

    // Score must be non-negative
    best_score as u64
}

/// Perform a local Smith-Waterman alignment.
///
/// Provides the locally optimal sequence alignment (1) using affine gap
/// penalties (2) with improvements and corrections by (3-4). Our implementation
/// adapts the corrected algorithm provided by [*Flouri et
/// al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html)(4).
///
/// See **[module citations](crate::alignment::sw#module-citations)**.
///
/// ## Complexity
///
/// Time: $O(mn)$
///
/// Space: $O(mn)$
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
/// # use zoe::{alignment::{Alignment, ScalarProfile, sw::sw_scalar_alignment}, data::{WeightMatrix, cigar::Cigar}};
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
pub fn sw_scalar_alignment<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> MaybeAligned<i32> {
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
            // matching is the default direction
            backtrack.move_to(r, c);

            let match_score = i32::from(query.matrix.get_weight(reference_base, query.seq[c]));
            h += match_score;

            let mut e = e_row[c];
            h = h.max(e).max(f).max(0);

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
        return MaybeAligned::Unmapped;
    }

    backtrack.to_alignment(best_score, r_end, c_end, reference.len(), query.seq.len())
}
