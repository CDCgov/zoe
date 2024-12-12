#![allow(
    clippy::cast_possible_wrap,
    clippy::cast_possible_truncation,
    clippy::needless_range_loop,
    unused_imports
)]
// TO-DO: revisit truncation issues

use super::*;
use crate::data::{
    byte_types::ByteMappings,
    mappings::{to_dna_index, to_dna_profile_index},
    types::{cigar::Cigar, Uint},
};
use std::{
    iter::zip,
    ops::{AddAssign, Shl},
    simd::{prelude::*, LaneCount, SupportedLaneCount},
};

/// Smith-Waterman algorithm, yielding the **ENDING** reference and query position and optimal score.
///
/// Locally optimal sequence alignment (1) for affine gaps (2) with improvements
/// and corrections by (3-4). Our implementation adapts the corrected algorithm
/// provided by [*Flouri et al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html) (4).
///
/// We use the affine gap formula of `W(k) = u(k-1) + v`. In order to use `W(k) =
/// uk + v`, simply pass gap open as the gap open plus gap extension.
///
/// ### Complexity
///
/// Time: *O(mn)*
///
/// Space: *O(mn)*
///
/// ### Citations
///
/// 1. Smith, Temple F. & Waterman, Michael S. (1981). "Identification of Common
///    Molecular Subsequences" (PDF). Journal of Molecular Biology. 147 (1):
///    195–197.
///
/// 2. Osamu Gotoh (1982). "An improved algorithm for matching biological
///    sequences". Journal of Molecular Biology. 162 (3): 705–708.
///
/// 3. Stephen F. Altschul & Bruce W. Erickson (1986). "Optimal sequence
///    alignment using affine gap costs". Bulletin of Mathematical Biology. 48
///    (5–6): 603–616.
///
/// 4. Tomáš Flouri, Kassian Kobert, Torbjørn Rognes, Alexandros
///    Stamatakis(2015). "Are all global alignment algorithms and implementations
///    correct?" bioRxiv 031500. doi: <https://doi.org/10.1101/031500>
///
#[must_use]
pub fn sw_score(
    reference: &[u8], query: &[u8], gap_open: i32, gap_extend: i32, substitution_matrix: &[[i32; 256]; 256],
) -> (usize, usize, i32) {
    let mut current = vec![ScoreCell::default(); query.len() + 1];
    let mut best_score = BestScore::new();

    for c in 1..=query.len() {
        // first row
        current[c].endgap_up = gap_open + (c as i32 - 1) * gap_extend;
    }

    for r in 1..=reference.len() {
        // first column
        let mut endgap_left = gap_open + (r as i32 - 1) * gap_extend;
        let mut diag = 0;

        for c in 1..=query.len() {
            // matching is the default direction

            let match_score = substitution_matrix[reference[r - 1] as usize][query[c - 1] as usize];
            let mut score = diag + match_score;
            let mut endgap_up = current[c].endgap_up;

            score = score.max(endgap_up).max(endgap_left).max(0);

            diag = current[c].matching;
            current[c].matching = score;
            best_score.add_score(r, c, score);

            score += gap_open;
            endgap_up += gap_extend;
            endgap_left += gap_extend;

            endgap_left = endgap_left.max(score);
            current[c].endgap_up = endgap_up.max(score);
        }
    }

    best_score.get_best_score()
}

/// An alternative Smith-Waterman implementation for benching and testing.
#[must_use]
#[allow(dead_code)]
pub(crate) fn sw_score_alt(
    reference: &[u8], query: &[u8], gap_open: i32, gap_extend: i32, matching: i32, mismatching: i32,
) -> i32 {
    struct SWCells {
        matching:  i32,
        endgap_up: i32,
        q_idx:     u8,
    }

    let mut best_score = 0;
    let reference: Vec<_> = reference.iter().copied().map(to_dna_index).collect();
    let mut current: Vec<_> = query
        .iter()
        .copied()
        .enumerate()
        .map(|(c,q)|
            // first row
            SWCells {
                matching:  0,
                endgap_up: gap_open + (c as i32) * gap_extend,
                q_idx:     to_dna_index(q),
            })
        .collect();

    for (r, r_idx) in reference.into_iter().enumerate() {
        // first column
        let mut endgap_left = gap_open + (r as i32) * gap_extend;
        let mut diag = 0;

        for curr in &mut current {
            // matching is the default direction

            let match_score = if r_idx == curr.q_idx { matching } else { mismatching };
            let mut score = diag + match_score;
            let mut endgap_up = curr.endgap_up;

            score = score.max(endgap_up).max(endgap_left).max(0);

            diag = curr.matching;
            curr.matching = score;
            best_score = best_score.max(score);

            score += gap_open;
            endgap_up += gap_extend;
            endgap_left += gap_extend;

            endgap_left = endgap_left.max(score);
            curr.endgap_up = endgap_up.max(score);
        }
    }

    best_score
}

/// Smith-Waterman alignment, yielding the reference starting position, cigar, and optimal score.
///
/// Locally optimal sequence alignment (1) for affine gaps (2) with improvements
/// and corrections by (3-4). Our implementation adapts the corrected algorithm
/// provided by [*Flouri et al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html) (4).
///
/// We use the affine gap formula of `W(k) = u(k-1) + v`. In order to use `W(k) =
/// uk + v`, simply pass gap open as the gap open plus gap extension.
///
/// ### Complexity
///
/// Time: *O(mn)*
///
/// Space: *O(mn)*
///
/// ### Citations
///
/// 1. Smith, Temple F. & Waterman, Michael S. (1981). "Identification of Common
///    Molecular Subsequences" (PDF). Journal of Molecular Biology. 147 (1):
///    195–197.
///
/// 2. Osamu Gotoh (1982). "An improved algorithm for matching biological
///    sequences". Journal of Molecular Biology. 162 (3): 705–708.
///
/// 3. Stephen F. Altschul & Bruce W. Erickson (1986). "Optimal sequence
///    alignment using affine gap costs". Bulletin of Mathematical Biology. 48
///    (5–6): 603–616.
///
/// 4. Tomáš Flouri, Kassian Kobert, Torbjørn Rognes, Alexandros
///    Stamatakis(2015). "Are all global alignment algorithms and implementations
///    correct?" bioRxiv 031500. doi: <https://doi.org/10.1101/031500>
///
#[must_use]
pub fn sw_alignment(
    reference: &[u8], query: &[u8], gap_open: i32, gap_extend: i32, substitution_matrix: &[[i32; 256]; 256],
) -> (usize, Cigar, i32) {
    let mut current = vec![ScoreCell::default(); query.len() + 1];
    let mut backtrack = BacktrackMatrix::new(reference.len() + 1, query.len() + 1);
    let mut best_score = BestScore::new();

    for c in 1..=query.len() {
        // first row
        current[c].endgap_up = gap_open + (c as i32 - 1) * gap_extend;
    }

    backtrack.stop();

    for r in 1..=reference.len() {
        // first column
        let mut endgap_left = gap_open + (r as i32 - 1) * gap_extend;
        let mut diag = 0;

        for c in 1..=query.len() {
            // matching is the default direction

            backtrack.move_to(r, c);
            let match_score = substitution_matrix[reference[r - 1] as usize][query[c - 1] as usize];
            let mut score = diag + match_score;
            let mut endgap_up = current[c].endgap_up;

            score = score.max(endgap_up).max(endgap_left);

            if endgap_up == score {
                backtrack.up();
            }

            if endgap_left == score {
                backtrack.left();
            }

            if score < 1 {
                score = 0;
                backtrack.stop();
            }

            diag = current[c].matching;
            current[c].matching = score;
            best_score.add_score(r, c, score);

            score += gap_open;
            endgap_up += gap_extend;
            endgap_left += gap_extend;

            if score != gap_open {
                if endgap_up > score {
                    backtrack.up_extending();
                }

                if endgap_left > score {
                    backtrack.left_extending();
                }
            }

            endgap_left = endgap_left.max(score);
            current[c].endgap_up = endgap_up.max(score);
        }
    }

    let mut states = AlignmentStates::new();
    let mut op = 0;
    let (mut r, mut c, best_score) = best_score.get_best_score();
    // soft clip 3'
    backtrack.move_to(r, c);
    states.soft_clip(query.len() - c);

    while !backtrack.is_stop() && r > 0 && c > 0 {
        if op == b'D' && backtrack.is_up_extending() {
            op = b'D';
            r -= 1;
        } else if op == b'I' && backtrack.is_left_extending() {
            op = b'I';
            c -= 1;
        } else if backtrack.is_up() {
            op = b'D';

            r -= 1;
        } else if backtrack.is_left() {
            op = b'I';
            c -= 1;
        } else {
            op = b'M';
            r -= 1;
            c -= 1;
        }
        states.add_state(op);
        backtrack.move_to(r, c);
    }

    // soft clip 5'
    states.soft_clip(c);
    let cigar = states.reverse().to_cigar();

    (r + 1, cigar, best_score)
}

/*
#[allow(dead_code, unused_variables, unused_mut, non_snake_case)]
fn sw_simd_score<T, const N: usize>(
    reference: &[u8], query: &QueryProfileStripedDNA<T, N>, gap_open: T, gap_extend: T,
    substitution_matrix: &[[i32; 256]; 256],
) -> (usize, usize, T)
where
    T: Uint,
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdUint<Scalar = T> + Shl<u8, Output = Simd<T, N>> + AddAssign<Simd<T, N>>, {
    let gap_open = Simd::from_array([gap_open; N]);
    let gap_extend = Simd::from_array([gap_extend; N]);
    let zeroes = Simd::from_array([T::zero(); N]);
    let mut max_scores = zeroes;

    let num_vecs = query.profile.len();
    let profile: &Vec<Simd<T, N>> = &query.profile;
    let mut load = vec![zeroes; num_vecs];
    let mut store = vec![zeroes; num_vecs];
    let mut e_scores = vec![zeroes; num_vecs];

    for ref_index in reference.iter().copied().map(ByteMappings::to_dna_index) {
        let mut F = zeroes;
        let mut H = store[num_vecs - 1] << 1;
        (load, store) = (store, load);

        let scores = &profile[(ref_index * num_vecs)..(ref_index * num_vecs + num_vecs)];
        for (score, (e_ref, (storing, loading))) in
            zip(scores, zip(e_scores.iter_mut(), zip(store.iter_mut(), load.iter_mut())))
        {
            H += *score;
            max_scores = max_scores.max(H);

            let mut E = *e_ref;
            H = H.max(E).max(F);
            *storing = H;

            H = H.saturating_sub(gap_open);

            E = E.saturating_sub(gap_extend);
            E = E.max(H);

            F = F.saturating_sub(gap_extend);
            F = F.max(H);

            // Store E; Load H
            *e_ref = E;
            H = *loading;
        }
    }

    let best = max_scores.saturating_sub(Simd::splat(query.bias)).reduce_max();
    (0, 0, best)
}
*/

#[cfg(test)]
pub(crate) mod test_data {
    // H5 HA: <https://www.ncbi.nlm.nih.gov/nuccore/NC_007362.1>
    pub(crate) static REFERENCE: &[u8] = b"GCAGGGGTATAATCTGTCAAAATGGAGAAAATAGTGCTTCTTCTTGCAATAGTCAGTCTTGTCAAAAGTGATCAGATTTGCATTGGTTACCATGCAAACAACTCGACAGAGCAGGTTGACACAATAATGGAAAAGAACGTTACTGTTACACATGCCCAAGACATACTGGAAAAGACACACAATGGGAAGCTCTGCGATCTAAATGGAGTGAAGCCTCTCATTTTGAGAGATTGTAGTGTAGCTGGATGGCTCCTCGGAAACCCTATGTGTGACGAATTCATCAATGTGCCGGAATGGTCTTACATAGTGGAGAAGGCCAGTCCAGCCAATGACCTCTGTTACCCAGGGGATTTCAACGACTATGAAGAACTGAAACACCTATTGAGCAGAACAAACCATTTTGAGAAAATTCAGATCATCCCCAAAAGTTCTTGGTCCAATCATGATGCCTCATCAGGGGTGAGCTCAGCATGTCCATACCATGGGAGGTCCTCCTTTTTCAGAAATGTGGTATGGCTTATCAAAAAGAACAGTGCATACCCAACAATAAAGAGGAGCTACAATAATACCAACCAAGAAGATCTTTTAGTACTGTGGGGGATTCACCATCCTAATGATGCGGCAGAGCAGACAAAGCTCTATCAAAACCCAACCACTTACATTTCCGTTGGAACATCAACACTGAACCAGAGATTGGTTCCAGAAATAGCTACTAGACCCAAAGTAAACGGGCAAAGTGGAAGAATGGAGTTCTTCTGGACAATTTTAAAGCCGAATGATGCCATCAATTTCGAGAGTAATGGAAATTTCATTGCTCCAGAATATGCATACAAAATTGTCAAGAAAGGGGACTCAGCAATTATGAAAAGTGAATTGGAATATGGTAACTGCAACACCAAGTGTCAAACTCCAATGGGGGCGATAAACTCTAGTATGCCATTCCACAACATACACCCCCTCACCATCGGGGAATGCCCCAAATATGTGAAATCAAACAGATTAGTCCTTGCGACTGGACTCAGAAATACCCCTCAGAGAGAGAGAAGAAGAAAAAAGAGAGGACTATTTGGAGCTATAGCAGGTTTTATAGAGGGAGGATGGCAGGGAATGGTAGATGGTTGGTATGGGTACCACCATAGCAATGAGCAGGGGAGTGGATACGCTGCAGACAAAGAATCCACTCAAAAGGCAATAGATGGAGTCACCAATAAGGTCAACTCGATCATTGACAAAATGAACACTCAGTTTGAGGCCGTTGGAAGGGAATTTAATAACTTGGAAAGGAGGATAGAGAATTTAAACAAGCAGATGGAAGACGGATTCCTAGATGTCTGGACTTATAATGCTGAACTTCTGGTTCTCATGGAAAATGAGAGAACTCTAGACTTTCATGACTCAAATGTCAAGAACCTTTATGACAAGGTCCGACTACAGCTTAGGGATAATGCAAAGGAGCTGGGTAATGGTTGTTTCGAGTTCTATCACAAATGTGATAATGAATGTATGGAAAGTGTAAAAAACGGAACGTATGACTACCCGCAGTATTCAGAAGAAGCAAGACTAAACAGAGAGGAAATAAGTGGAGTAAAATTGGAATCAATGGGAACTTACCAAATACTGTCAATTTATTCAACAGTGGCGAGTTCCCTAGCACTGGCAATCATGGTAGCTGGTCTATCTTTATGGATGTGCTCCAATGGATCGTTACAATGCAGAATTTGCATTTAAATTTGTGAGTTCAGATTGTAGTTAAAAACACC";

    /// H1 HA1: <https://www.ncbi.nlm.nih.gov/nuccore/NC_026433.1>
    pub(crate) static QUERY: &[u8] = b"gacacattatgtataggttatcatgcgaacaattcaacagacactgtagacacagtactagaaaagaatgtaacagtaacacactctgttaaccttctagaagacaagcataacgggaaactatgcaaactaagaggggtagccccattgcatttgggtaaatgtaacattgctggctggatcctgggaaatccagagtgtgaatcactctccacagcaagctcatggtcctacattgtggaaacacctagttcagacaatggaacgtgttacccaggagatttcatcgattatgaggagctaagagagcaattgagctcagtgtcatcatttgaaaggtttgagatattccccaagacaagttcatggcccaatcatgactcgaacaaaggtgtaacggcagcatgtcctcatgctggagcaaaaagcttctacaaaaatttaatatggctagttaaaaaaggaaattcatacccaaagctcagcaaatcctacattaatgataaagggaaagaagtcctcgtgctatggggcattcaccatccatctactagtgctgaccaacaaagtctctatcagaatgcagatgcatatgtttttgtggggtcatcaagatacagcaagaagttcaagccggaaatagcaataagacccaaagtgagggatcaagaagggagaatgaactattactggacactagtagagccgggagacaaaataacattcgaagcaactggaaatctagtggtaccgagatatgcattcgcaatggaaagaaatgctggatctggtattatcatttcagatacaccagtccacgattgcaatacaacttgtcaaacacccaagggtgctataaacaccagcctcccatttcagaatatacatccgatcacaattggaaaatgtccaaaatatgtaaaaagcacaaaattgagactggccacaggattgaggaatatcccgtctattcaatctaga";

    pub(crate) const GAP_OPEN: i32 = -10;
    pub(crate) const GAP_EXTEND: i32 = -1;

    pub(crate) const SM: [[i32; 256]; 256] = {
        const MISMATCH: i32 = -5;
        const MATCH: i32 = 2;
        let mut bytes = [[0i32; 256]; 256];
        let bases = b"AGCTatcg";

        let mut i = 0;
        while i < bases.len() {
            let mut j = 0;
            while j < bases.len() {
                bytes[bases[i] as usize][bases[j] as usize] =
                    if bases[i].to_ascii_uppercase() == bases[j].to_ascii_uppercase() {
                        MATCH
                    } else {
                        MISMATCH
                    };
                j += 1;
            }
            i += 1;
        }
        bytes
    };
}

#[cfg(test)]
mod test {
    use super::{test_data::*, *};

    #[test]
    fn sw() {
        let (r, cigar, score) = sw_alignment(REFERENCE, QUERY, GAP_OPEN, GAP_EXTEND, &SM);
        assert_eq!((337, 37), (r, score));

        let matches = cigar.match_length();

        let (r, _, score) = sw_score(REFERENCE, QUERY, GAP_OPEN, GAP_EXTEND, &SM);
        assert_eq!((337, 37), (r - matches + 1, score));

        let score = sw_score_alt(REFERENCE, QUERY, GAP_OPEN, GAP_EXTEND, 2, -5);
        assert_eq!(37, score);
    }
}

#[cfg(test)]
mod bench {
    extern crate test;
    use super::{test_data::*, *};
    use test::Bencher;

    #[bench]
    fn sw_alignment_scalar(b: &mut Bencher) {
        b.iter(|| sw_alignment(REFERENCE, QUERY, GAP_OPEN, GAP_EXTEND, &SM));
    }

    #[bench]
    fn sw_score_scalar(b: &mut Bencher) {
        b.iter(|| sw_score(REFERENCE, QUERY, GAP_OPEN, GAP_EXTEND, &SM));
    }

    #[bench]
    fn sw_score_scalar_alt(b: &mut Bencher) {
        b.iter(|| sw_score_alt(REFERENCE, QUERY, GAP_OPEN, GAP_EXTEND, 2, -5));
    }
}
