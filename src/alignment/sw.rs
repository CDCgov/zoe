#![allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation, clippy::needless_range_loop)]
// TO-DO: revisit truncation issues

use super::*;
use crate::{
    data::types::{cigar::Cigar, Uint},
    simd::SimdExt,
};
use std::{
    ops::{AddAssign, Shl},
    simd::{prelude::*, LaneCount, SimdElement, SupportedLaneCount},
};

/// Smith-Waterman algorithm, yielding the optimal score.
///
/// Provides the locally optimal sequence alignment (1) for affine gaps (2). Our
/// implementation adapts the algorithm provided by [*Flouri et
/// al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html)
/// (3).
///
/// We use the affine gap formula of $W(k) = u(k-1) + v$. In order to use $W(k) =
/// uk + v$, simply pass gap open as the gap open plus gap extension.
///
/// ### Complexity
///
/// Time: $O(mn)$
///
/// Space: $O(m)$ where $m$ is the length of the query
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
/// 3. Tomáš Flouri, Kassian Kobert, Torbjørn Rognes, Alexandros
///    Stamatakis(2015). "Are all global alignment algorithms and implementations
///    correct?" bioRxiv 031500. doi: <https://doi.org/10.1101/031500>
///
#[must_use]
pub fn sw_scalar_score<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> i32 {
    struct SWCells {
        matching:   i32,
        endgap_up:  i32,
        query_base: u8,
    }

    let mut best_score = 0;
    let mut current: Vec<_> = query
        .query
        .iter()
        .copied()
        .enumerate()
        .map(|(c, q)|
            // first row
            SWCells {
                matching:  0,
                endgap_up: query.gap_open + (c as i32) * query.gap_extend,
                query_base:     q
            })
        .collect();

    for (r, reference_base) in reference.iter().copied().enumerate() {
        // first column
        let mut endgap_left = query.gap_open + (r as i32) * query.gap_extend;
        let mut diag = 0;

        for curr in &mut current {
            // matching is the default direction

            let match_score = i32::from(query.matrix.get_weight(reference_base, curr.query_base));
            let mut score = diag + match_score;
            let mut endgap_up = curr.endgap_up;

            score = score.max(endgap_up).max(endgap_left).max(0);

            diag = curr.matching;
            curr.matching = score;
            best_score = best_score.max(score);

            score += query.gap_open;
            endgap_up += query.gap_extend;
            endgap_left += query.gap_extend;

            endgap_left = endgap_left.max(score);
            curr.endgap_up = endgap_up.max(score);
        }
    }

    best_score
}

/// Smith-Waterman alignment, yielding the reference starting position, cigar,
/// and optimal score.
///
/// Provides the locally optimal sequence alignment (1) for affine gaps (2) with
/// improvements and corrections by (3-4). Our implementation adapts the
/// corrected algorithm provided by [*Flouri et
/// al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html)
/// (4).
///
/// We use the affine gap formula of $W(k) = u(k-1) + v$. In order to use $W(k) =
/// uk + v$, simply pass gap open as the gap open plus gap extension.
///
/// ### Complexity
///
/// Time: $O(mn)$
///
/// Space: $O(mn)$
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
pub fn sw_scalar_alignment<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> (usize, Cigar, i32) {
    let mut current = vec![ScoreCell::default(); query.query.len() + 1];
    let mut backtrack = BacktrackMatrix::new(reference.len() + 1, query.query.len() + 1);
    let mut best_score = BestScore::new();

    for c in 1..=query.query.len() {
        // first row
        current[c].endgap_up = query.gap_open + (c as i32 - 1) * query.gap_extend;
    }

    backtrack.stop();

    for r in 1..=reference.len() {
        // first column
        let mut endgap_left = query.gap_open + (r as i32 - 1) * query.gap_extend;
        let mut diag = 0;

        for c in 1..=query.query.len() {
            // matching is the default direction

            backtrack.move_to(r, c);
            let match_score = i32::from(query.matrix.get_weight(reference[r - 1], query.query[c - 1]));
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

            score += query.gap_open;
            endgap_up += query.gap_extend;
            endgap_left += query.gap_extend;

            if score != query.gap_open {
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
    states.soft_clip(query.query.len() - c);

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

/// Smith-Waterman algorithm (vectorized), yielding the optimal score.
///
/// Provides the locally optimal sequence alignment (1) for affine gaps (2). We
/// adapt Farrar's striped SIMD implemention (3) for portable SIMD.
///
/// We use the affine gap formula of $W(k) = u(k-1) + v$. In order to use $W(k) =
/// uk + v$, simply pass gap open as the gap open plus gap extension.
///
/// ### Complexity
///
/// Time: $O(mn)$, with the average case as $mn/N$ for $N$ SIMD lanes
///
/// Space: $O(n)$
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
/// 3. Farrar, Michael (2006). "Striped Smith-Waterman speeds database searches
///    six times over other SIMD implementations". Bioinformatics, 23(2),
///    156-161. doi: <https://doi.org/10.1093/bioinformatics/btl582>
///
#[allow(non_snake_case)]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn sw_simd_score<T, const N: usize, const S: usize>(reference: &[u8], query: &StripedProfile<T, N, S>) -> Option<T>
where
    T: Uint + From<u8> + Default + PartialEq + std::ops::Add<Output = T> + std::fmt::Debug,
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdUint<Scalar = T>
        + Shl<T, Output = Simd<T, N>>
        + AddAssign<Simd<T, N>>
        + SimdOrd
        + SimdPartialEq<Mask = Mask<<T as SimdElement>::Mask, N>>, {
    let num_vecs = query.number_vectors();
    let profile: &Vec<Simd<T, N>> = &query.profile;

    let zero = T::zero();
    let gap_opens = Simd::splat(query.gap_open);
    let gap_extends = Simd::splat(query.gap_extend);
    let zeroes = Simd::splat(T::zero());
    let biases = Simd::splat(query.bias);

    let mut load = vec![zeroes; num_vecs];
    let mut store = vec![zeroes; num_vecs];
    let mut e_scores = vec![zeroes; num_vecs];
    let mut max_scores = zeroes; // Minimum value for unsigned

    let map = query.mapping;

    for ref_index in reference.iter().copied().map(|r| map.to_index(r)) {
        let mut F = zeroes;
        let mut H = store[num_vecs - 1].shift_elements_right_z::<1>(zero);

        (load, store) = (store, load);

        // This statement helps with bounds checks.
        let scores_vec = &profile[(ref_index * num_vecs)..(ref_index * num_vecs + num_vecs)];

        for j in 0..num_vecs {
            let mut E = e_scores[j];

            H = H.saturating_add(scores_vec[j]).saturating_sub(biases);

            max_scores = max_scores.simd_max(H);

            H = H.simd_max(E).simd_max(F);
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
        F = F.shift_elements_right_z::<1>(zero);

        // (F - (H - Go)) --> 0
        let mut mask = F.saturating_sub(H.saturating_sub(gap_opens)).simd_eq(zeroes);
        while !mask.all() {
            // Update & save H
            H = H.simd_max(F);
            store[j] = H;

            // Update E in case H is greater
            // TODO: Remove this when we determine an API to guarantee it isn't necessary
            e_scores[j] = e_scores[j].simd_max(H.saturating_sub(gap_opens));

            F = F.saturating_sub(gap_extends);

            j += 1;
            if j >= num_vecs {
                j = 0;
                F = F.shift_elements_right_z::<1>(zero);
            }

            // New J here
            H = store[j];

            mask = F.saturating_sub(H.saturating_sub(gap_opens)).simd_eq(zeroes);
        }
    }

    let best = max_scores.reduce_max();

    // If we would have overflowed, return none, otherwise return the best score
    // We add one because we care if the value is equal to the MAX.
    best.checked_addition(query.bias + T::one()).map(|_| best)
}

#[cfg(test)]
pub(crate) mod test_data {
    use super::ScalarProfile;
    use crate::data::{BiasedWeightMatrix, SimpleWeightMatrix};
    use std::sync::LazyLock;

    pub(crate) static REFERENCE: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/NC_007362.1.txt")); // H5 HA
    pub(crate) static QUERY: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/NC_026433.1.txt")); // H1 HA1

    pub(crate) static WEIGHTS: SimpleWeightMatrix<5> = SimpleWeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
    pub(crate) static BIASED_WEIGHTS: BiasedWeightMatrix<5> = WEIGHTS.into_biased_matrix();
    pub(crate) static SCALAR_PROFILE: LazyLock<ScalarProfile<5>> =
        LazyLock::new(|| ScalarProfile::new(QUERY, WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap());

    pub(crate) const GAP_OPEN: u8 = 10;
    pub(crate) const GAP_EXTEND: u8 = 1;
}

#[cfg(test)]
mod test {
    use super::{test_data::*, *};
    use crate::alignment::profile::StripedProfile;
    use crate::data::constants::matrices::SimpleWeightMatrix;
    use crate::data::mappings::DNA_PROFILE_MAP;

    #[test]
    fn sw() {
        let weights = SimpleWeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
        let profile = ScalarProfile::new(QUERY, weights, GAP_OPEN, GAP_EXTEND).unwrap();
        let (r, _, score) = sw_scalar_alignment(REFERENCE, &profile);
        assert_eq!((337, 37), (r, score));

        let score = sw_scalar_score(REFERENCE, &profile);
        assert_eq!(37, score);

        let v: Vec<_> = std::iter::repeat(b'A').take(100).collect();
        let profile = ScalarProfile::new(&v, weights, GAP_OPEN, GAP_EXTEND).expect("Sequence is non-empty");
        let score = sw_scalar_score(&v, &profile);
        assert_eq!(200, score);
    }

    #[test]
    fn sw_simd() {
        let matrix = SimpleWeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N')).into_biased_matrix();
        let profile = matrix
            .to_striped_profile::<u8, 16>(REFERENCE, GAP_OPEN, GAP_EXTEND)
            .expect("Sequence is non-empty");
        let score = profile.smith_waterman_score(QUERY);
        assert_eq!(Some(37), score);
    }

    #[test]
    fn sw_simd_poly_a() {
        let v: Vec<_> = std::iter::repeat(b'A').take(100).collect();
        let matrix = SimpleWeightMatrix::new_dna_matrix(2, -5, Some(b'N')).into_biased_matrix();
        let profile = StripedProfile::<u16, 16, 5>::new(&v, &matrix, GAP_OPEN, GAP_EXTEND).expect("Sequence is non-empty");
        let score = profile.smith_waterman_score(&v);
        assert_eq!(Some(200), score);
    }

    #[test]
    fn sw_simd_single() {
        let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
        let matrix = SimpleWeightMatrix::<5>::new_dna_matrix(2, -5, Some(b'N')).into_biased_matrix();
        let profile = matrix
            .to_striped_profile::<u16, 16>(v, GAP_OPEN, GAP_EXTEND)
            .expect("Sequence is non-empty");
        let score = profile.smith_waterman_score(v);
        assert_eq!(Some(3372), score);
    }

    #[test]
    fn sw_simd_regression() {
        let query = b"AGA";
        let reference = b"AA";
        let matrix = SimpleWeightMatrix::new(&DNA_PROFILE_MAP, 10, -10, Some(b'N')).into_biased_matrix();
        let profile = matrix
            .to_striped_profile::<u16, 4>(query, 5, 5)
            .expect("Sequence is non-empty");
        let score = profile.smith_waterman_score(reference);
        assert_eq!(Some(15), score);
    }

    #[test]
    fn sw_simd_overflow_check() {
        let query = b"AAAA";
        let reference = b"AAAA";

        let matrix = SimpleWeightMatrix::new(&DNA_PROFILE_MAP, 127, 0, Some(b'N')).into_biased_matrix();
        let profile = matrix
            .to_striped_profile::<u8, 8>(query, GAP_OPEN, GAP_EXTEND)
            .expect("Sequence is non-empty");
        let score = profile.smith_waterman_score(reference);
        assert!(score.is_none());
    }

    #[test]
    fn sw_simd_profile_set() {
        let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
        let matrix = SimpleWeightMatrix::<5>::new_dna_matrix(2, -5, Some(b'N')).into_biased_matrix();
        let profile = LocalProfile::<16, 5>::new(&v, &matrix, GAP_OPEN, GAP_EXTEND).expect("Sequence is non-empty");
        let score = profile.smith_waterman_score_from_u8(&v);
        assert_eq!(Some(3372), score);
    }
}

#[cfg(test)]
mod bench {
    extern crate test;
    use super::{test_data::*, *};
    use crate::alignment::profile::StripedProfile;
    use test::Bencher;

    #[bench]
    fn sw_alignment_scalar(b: &mut Bencher) {
        let query_profile = &*SCALAR_PROFILE;
        b.iter(|| sw_scalar_alignment(REFERENCE, query_profile));
    }

    #[bench]
    fn sw_score_scalar(b: &mut Bencher) {
        let query_profile = &*SCALAR_PROFILE;
        b.iter(|| sw_scalar_score(REFERENCE, query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_n16u8(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u8, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_n16u16(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u16, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_n32u8(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u8, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_n32u16(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u16, 32, 5>(QUERY, &query_profile));
    }
}
