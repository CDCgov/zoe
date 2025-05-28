#![allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation, clippy::needless_range_loop)]
// TODO: revisit truncation issues

use super::*;
use crate::{data::cigar::Ciglet, math::AnyInt, simd::SimdAnyInt};
use std::{
    cmp::Ordering::{Equal, Greater, Less},
    ops::Add,
    simd::{LaneCount, SimdElement, SupportedLaneCount, prelude::*},
};

/// Compute the Smith-Waterman score for the alignment given by `cigar` between
/// `query` and the reference slice corresponding to the
/// [`Alignment::ref_range`].
///
/// ## Errors
///
/// - A `cigar` must be a valid CIGAR string
/// - All of `query`, `ref_in_alignment`, and `cigar` must be fully consumed
/// - The final score should be nonnegative
pub fn sw_score_from_path<const S: usize>(
    ciglets: impl IntoIterator<Item = Ciglet>, ref_in_alignment: &[u8], query: &ScalarProfile<S>,
) -> Result<u64, ScoringError> {
    let mut score = 0;
    let mut r: usize = 0;
    let mut q = 0;

    for Ciglet { inc, op } in ciglets {
        match op {
            b'M' | b'=' | b'X' => {
                for _ in 0..inc {
                    let Some(reference_base) = ref_in_alignment.get(r).copied() else {
                        return Err(ScoringError::ReferenceEnded);
                    };
                    let Some(query_base) = query.seq.get(q).copied() else {
                        return Err(ScoringError::QueryEnded);
                    };
                    score += i32::from(query.matrix.get_weight(reference_base, query_base));
                    q += 1;
                    r += 1;
                }
            }
            b'I' => {
                score += query.gap_open + query.gap_extend * (inc - 1) as i32;
                q += inc;
            }
            b'D' => {
                score += query.gap_open + query.gap_extend * (inc - 1) as i32;
                r += inc;
            }
            b'S' => q += inc,
            b'N' => r += inc,
            b'H' | b'P' => {}
            op => return Err(ScoringError::InvalidCigarOp(op)),
        }
    }

    match q.cmp(&query.seq.len()) {
        Less => return Err(ScoringError::FullQueryNotUsed),
        Greater => return Err(ScoringError::QueryEnded),
        Equal => {}
    }

    match r.cmp(&ref_in_alignment.len()) {
        Less => return Err(ScoringError::FullReferenceNotUsed),
        Greater => return Err(ScoringError::ReferenceEnded),
        Equal => {}
    }

    match u64::try_from(score) {
        Ok(score) => Ok(score),
        Err(_) => Err(ScoringError::NegativeScore(score)),
    }
}

/// Smith-Waterman algorithm, yielding the optimal score.
///
/// Provides the locally optimal sequence alignment (1) using affine gap
/// penalties (2). Our implementation adapts the algorithm provided by [*Flouri
/// et al.*](https://cme.h-its.org/exelixis/web/software/alignment/correct.html)
/// (3).
///
/// We use the affine gap formula, $W(k) = u(k-1) + v$, where $k$ is the gap
/// length, $u$ is the gap extend penalty, $v$ is the gap open penalty, and
/// $W(k)$ is the total penalty for the gap. In order to use $W(k) = uk + v$,
/// simply pass gap open as the gap open plus gap extension.
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
///
/// ## Complexity
///
/// Time: $O(mn)$
///
/// Space: $O(m)$ where $m$ is the length of the query
///
/// ## Citations
///
/// 1. Smith, Temple F. & Waterman, Michael S. (1981). "Identification of Common
///    Molecular Subsequences" (PDF). Journal of Molecular Biology. 147 (1):
///    195–197.
///
/// 2. Osamu Gotoh (1982). "An improved algorithm for matching biological
///    sequences". Journal of Molecular Biology. 162 (3): 705–708.
///
/// 3. Tomáš Flouri, Kassian Kobert, Torbjørn Rognes, Alexandros
///    Stamatakis(2015). "Are all global alignment algorithms and
///    implementations correct?" bioRxiv 031500. doi:
///    <https://doi.org/10.1101/031500>
///
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
/// We use the affine gap formula, $W(k) = u(k-1) + v$, where $k$ is the gap
/// length, $u$ is the gap extend penalty, $v$ is the gap open penalty, and
/// $W(k)$ is the total penalty for the gap. In order to use $W(k) = uk + v$,
/// simply pass gap open as the gap open plus gap extension.
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
/// let alignment = sw_scalar_alignment(&reference, &profile);
/// assert_eq!(alignment.ref_range.start, 3);
/// assert_eq!(alignment.states, Cigar::from_slice_unchecked("5M1D4M"));
/// assert_eq!(alignment.score, 27);
/// ```
///
/// ## Complexity
///
/// Time: $O(mn)$
///
/// Space: $O(mn)$
///
/// ## Citations
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
///    Stamatakis(2015). "Are all global alignment algorithms and
///    implementations correct?" bioRxiv 031500. doi:
///    <https://doi.org/10.1101/031500>
///
#[must_use]
pub fn sw_scalar_alignment<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> Alignment<i32> {
    let (mut best_score, mut r_end, mut c_end) = (0, 0, 0);
    let mut h_row = vec![0; query.seq.len()];
    let mut e_row = vec![query.gap_open; query.seq.len()];

    let mut backtrack = BacktrackMatrix::new(reference.len(), query.seq.len());

    //backtrack.stop();

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

            if h >= best_score {
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

    let mut states = AlignmentStates::new();
    let mut op = 0;

    backtrack.move_to(r_end, c_end); // 0-based move to max

    // 1-based as though we had a padded matrix
    r_end += 1;
    c_end += 1;

    let (mut r, mut c) = (r_end, c_end);

    // soft clip 3'
    states.soft_clip(query.seq.len() - c);

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
        backtrack.move_to(r.saturating_sub(1), c.saturating_sub(1)); // 0-based next position
    }

    // soft clip 5'
    states.soft_clip(c);
    states.make_reverse();

    // r and c are decremented and becomes 0-based
    Alignment {
        score: best_score,
        ref_range: r..r_end,
        query_range: c..c_end,
        states,
        ref_len: reference.len(),
        query_len: query.seq.len(),
    }
}

/// Smith-Waterman algorithm (vectorized), yielding the optimal score.
///
/// Provides the locally optimal sequence alignment (1) using affine gap
/// penalties (2). We adapt Farrar's striped SIMD implemention (3) for portable
/// SIMD.
///
/// We use the affine gap formula, $W(k) = u(k-1) + v$, where $k$ is the gap
/// length, $u$ is the gap extend penalty, $v$ is the gap open penalty, and
/// $W(k)$ is the total penalty for the gap. In order to use $W(k) = uk + v$,
/// simply pass gap open as the gap open plus gap extension.
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
///
/// ## Complexity
///
/// Time: $O(mn)$, with the average case as $mn/N$ for $N$ SIMD lanes
///
/// Space: $O(n)$
///
/// ## Limitations
///
/// The SIMD algorithm used here may not perform as well when both the query and
/// reference are short. If the query and reference are both shorter than 25
/// bases, ensure that this algorithm is called with `N=16` or less. If the
/// query and reference are both shorter than 10 bases, consider using the
/// scalar algorithm [`sw_scalar_score`].
///
/// ## Citations
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
pub fn sw_simd_score<T, const N: usize, const S: usize>(reference: &[u8], query: &StripedProfile<T, N, S>) -> Option<u64>
where
    T: AnyInt + SimdElement,
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
            // Update & save H
            H = H.simd_max(F);
            store[j] = H;

            F = F.saturating_sub(gap_extends);

            j += 1;
            if j >= num_vecs {
                j = 0;
                F = F.shift_elements_right::<1>(min);
            }

            // New J here
            H = store[j];

            mask = F.simd_gt(H.saturating_sub(gap_opens));
        }
    }

    let best = max_scores.reduce_max();

    if T::SIGNED {
        // Map best score to an unsigned range. Note that: MAX+1 = abs(MIN). If
        // we would have overflowed (saturated), return None, otherwise return
        // the best score.
        (best < T::MAX).then(|| (T::MAX.cast_as::<u64>() + 1).wrapping_add_signed(best.cast_as::<i64>()))
    } else {
        // If we would have overflowed, return None, otherwise return the best
        // score. We add one because we care if the value is equal to the MAX.
        best.checked_add(query.bias + T::ONE).map(|_| best.cast_as::<u64>())
    }
}

#[cfg(test)]
pub(crate) mod test_data {
    use super::ScalarProfile;
    use crate::data::WeightMatrix;
    use std::sync::LazyLock;

    pub(crate) static REFERENCE: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/NC_007362.1.txt")); // H5 HA
    pub(crate) static QUERY: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/NC_026433.1.txt")); // H1 HA1

    pub(crate) static WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
    pub(crate) static BIASED_WEIGHTS: WeightMatrix<u8, 5> = WEIGHTS.into_biased_matrix();
    pub(crate) static SCALAR_PROFILE: LazyLock<ScalarProfile<5>> =
        LazyLock::new(|| ScalarProfile::new(QUERY, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap());

    pub(crate) const GAP_OPEN: i8 = -10;
    pub(crate) const GAP_EXTEND: i8 = -1;
}

#[cfg(test)]
mod test {
    use super::{test_data::*, *};
    use crate::alignment::profile::StripedProfile;
    use crate::data::constants::matrices::WeightMatrix;
    use crate::data::mappings::DNA_PROFILE_MAP;

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn sw_verify_score_from_path() {
        let reference = b"ATTCCTTTTGCCGGG";
        let weights: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(3, -1, Some(b'N'));
        let profile = ScalarProfile::new(b"ATTGCGCCCGG", &weights, -4, -1).unwrap();

        let Alignment {
            score,
            ref_range,
            states,
            ..
        } = sw_scalar_alignment(reference, &profile);

        assert_eq!(Ok(score as u64), sw_score_from_path(&states, &reference[ref_range], &profile));
    }

    #[test]
    fn sw() {
        let weights: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
        let profile = ScalarProfile::new(QUERY, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
        let Alignment { score, ref_range, .. } = sw_scalar_alignment(REFERENCE, &profile);
        assert_eq!((336, 37), (ref_range.start, score));

        let score = sw_scalar_score(REFERENCE, &profile);
        assert_eq!(37, score);

        let v: Vec<_> = std::iter::repeat_n(b'A', 100).collect();
        let profile = ScalarProfile::new(&v, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = sw_scalar_score(&v, &profile);
        assert_eq!(200, score);
    }

    #[test]
    fn sw_t_u_check() {
        let profile = ScalarProfile::new(b"ACGTUNacgtun", &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = sw_scalar_score(b"ACGTTNACGTTN", &profile);
        assert_eq!(score, 20);

        let profile = StripedProfile::<u16, 16, 5>::new(b"ACGTUNacgtun", &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = sw_simd_score(b"ACGTTNACGTTN", &profile);
        assert_eq!(score, Some(20));

        let profile = StripedProfile::<i16, 16, 5>::new(b"ACGTUNacgtun", &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = sw_simd_score(b"ACGTTNACGTTN", &profile);
        assert_eq!(score, Some(20));
    }

    #[test]
    fn sw_simd() {
        let matrix_i = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
        let matrix_u = matrix_i.into_biased_matrix();

        let profile = StripedProfile::<u8, 16, 5>::new(REFERENCE, &matrix_u, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = profile.smith_waterman_score(QUERY);
        assert_eq!(Some(37), score);

        let profile = StripedProfile::<i16, 16, 5>::new(REFERENCE, &matrix_i, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = profile.smith_waterman_score(QUERY);
        assert_eq!(Some(37), score);
    }

    #[test]
    fn sw_simd_poly_a() {
        let v: Vec<_> = std::iter::repeat_n(b'A', 100).collect();
        let matrix = WeightMatrix::new_dna_matrix(2, -5, Some(b'N')).into_biased_matrix();
        let profile = StripedProfile::<u16, 16, 5>::new(&v, &matrix, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = profile.smith_waterman_score(&v);
        assert_eq!(Some(200), score);
    }

    #[test]
    fn sw_simd_single() {
        let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
        let matrix = WeightMatrix::<i8, 5>::new_dna_matrix(2, -5, Some(b'N')).into_biased_matrix();
        let profile = StripedProfile::<u16, 16, 5>::new(v, &matrix, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = profile.smith_waterman_score(v);
        assert_eq!(Some(3372), score);
    }

    #[test]
    fn sw_simd_regression() {
        let query = b"AGA";
        let reference = b"AA";
        let matrix = WeightMatrix::new(&DNA_PROFILE_MAP, 10, -10, Some(b'N')).into_biased_matrix();
        let profile = StripedProfile::<u16, 4, 5>::new(query, &matrix, -5, -5).unwrap();
        let score: Option<u64> = profile.smith_waterman_score(reference);
        assert_eq!(Some(15), score);
    }

    #[test]
    fn sw_simd_overflow_check() {
        let query = b"AAAA";
        let reference = b"AAAA";

        let matrix = WeightMatrix::new(&DNA_PROFILE_MAP, 127, 0, Some(b'N')).into_biased_matrix();
        let profile = StripedProfile::<u8, 8, 5>::new(query, &matrix, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = profile.smith_waterman_score(reference);
        assert!(score.is_none());
    }

    #[test]
    fn sw_simd_profile_set() {
        let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
        let matrix_i = WeightMatrix::<i8, 5>::new_dna_matrix(2, -5, Some(b'N'));

        let profile = LocalProfiles::<16, 8, 4, 2, 5>::new(&v, &matrix_i, GAP_OPEN, GAP_EXTEND).unwrap();
        let score = profile.smith_waterman_score_from_i8(&v);
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
    fn sw_simd_no_profile_16n08u(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u8, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_16n16u(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u16, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_32n08u(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u8, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_32n16u(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u16, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_16n08i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i8, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_16n16i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i16, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_32n08i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i8, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_no_profile_32n16i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i16, 32, 5>(QUERY, &query_profile));
    }
}
