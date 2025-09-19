use super::*;
use crate::alignment::BackTrackable;
use std::ops::Add;

/// Banded Smith-Waterman algorithm, yielding a local alignment.
///
/// This implementation restricts the search space to a diagonal band and
/// performs traceback to generate a local alignment with CIGAR string.
/// Note that this is a heuristic approach that may not find the globally
/// optimal alignment due to the band constraint.
///
/// ## Arguments
///
/// * `reference` - The reference sequence
/// * `query` - The query profile
/// * `band_width` - Number of positions away from diagonal to search
///
/// ## Example
///
/// ```
/// # use zoe::{
/// #     alignment::{Alignment, ScalarProfile, sw::sw_banded_alignment},
/// #     data::{matrices::WeightMatrix, cigar::Cigar}
/// # };
/// let reference: &[u8] = b"GGCCACAGGATTGAG";
/// let query: &[u8] = b"CTCAGATTG";
///
/// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
/// const GAP_OPEN: i8 = -3;
/// const GAP_EXTEND: i8 = -1;
/// const BAND_WIDTH: usize = 5;
///
/// let profile = ScalarProfile::<5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
/// let alignment = sw_banded_alignment(reference, &profile, BAND_WIDTH).unwrap();
/// ```
#[must_use]
#[allow(clippy::cast_sign_loss)]
pub fn sw_banded_alignment<const S: usize>(
    reference: &[u8], query: &ScalarProfile<S>, band_width: usize,
) -> MaybeAligned<Alignment<u32>> {
    if reference.is_empty() {
        return MaybeAligned::Unmapped;
    }

    let (mut best_score, mut r_end, mut c_end) = (0, 0, 0);
    let q_len = query.seq.len();

    let mut h_row = vec![0i32; q_len];
    let mut e_row = vec![query.gap_open; q_len];
    let mut backtrack = BandedBacktrackMatrix::new(reference.len(), band_width);

    let mut h_store = 0;

    for (r, reference_base) in reference.iter().copied().enumerate() {
        let mut f = query.gap_open;
        let mut h = h_store;

        let start_col = r.saturating_sub(band_width);
        let end_col = (r + band_width + 1).min(q_len);

        for c in start_col..end_col {
            backtrack.move_to(r, c);

            // matching is the default direction
            let match_score = i32::from(query.matrix.get_weight(reference_base, query.seq[c]));
            h += match_score;

            let mut e = e_row[c];
            h = h.max(e).max(f).max(0);

            // TODO: THIS IS INEFFICIENT
            if c + band_width == r {
                h_store = h;
            }

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
        MaybeAligned::Some(backtrack.to_alignment(best_score as u32, r_end, c_end, reference.len(), q_len))
    }
}
