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
        query_range: _,
        states,
        ..
    } = sw_scalar_alignment(reference, &profile).unwrap();

    assert_eq!(Ok(score as u64), sw_score_from_path(&states, &reference[ref_range], &profile));
}

#[test]
fn sw() {
    let weights: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
    let profile = ScalarProfile::new(QUERY, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let Alignment { score, ref_range, .. } = sw_scalar_alignment(REFERENCE, &profile).unwrap();
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
fn sw_simd_aln() {
    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(REFERENCE, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<i16, 8, 5>::new(REFERENCE, &weights, GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(QUERY);

    let aln_scalar = profile_scalar.smith_waterman_alignment(QUERY).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(QUERY).unwrap();
    assert_eq!(score, Some(aln_scalar.score));
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln2() {
    let reference = b"AAACTA";
    let query = b"TTTAG";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<i8, 2, 5>::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln3() {
    let reference = b"AAAAAAAAAA";
    let query = b"AAAAAATAAA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<i8, 4, 5>::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln4() {
    let reference = b"TAAAA";
    let query = b"CCCCA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<i8, 4, 5>::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln5() {
    let reference = b"TCCCC";
    let query = b"CCCCC";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<i8, 4, 5>::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln6() {
    let reference = b"GCTTTTC";
    let query = b"CCCCT";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<i8, 4, 5>::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));
    assert_eq!(aln_scalar, aln_simd);
}
#[test]
fn sw_simd_aln7() {
    let reference = b"TTGTTTTTTTTTGTT";
    let query = b"TTTTTGTTTTCTTTTTTGTTTA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<i8, 16, 5>::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));

    let recaculated_score = sw_score_from_path(&aln_simd.states, &reference[aln_simd.ref_range.clone()], &profile_scalar);
    assert_eq!(Ok(score.unwrap()), recaculated_score, "for {aln_simd:?}",);
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln8() {
    let reference = b"TTTTTGTTTGGGAAAAATTCTT";
    let query = b"TTGTTTTGGGGAAAAA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<u8, 8, 5>::new(query, &weights.into_biased_matrix(), GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));

    let recaculated_score = sw_score_from_path(&aln_simd.states, &reference[aln_simd.ref_range.clone()], &profile_scalar);
    assert_eq!(
        Ok(score.unwrap()),
        recaculated_score,
        "{aln_simd:?} (v {cigar})",
        cigar = aln_scalar.states
    );
    assert_eq!(aln_scalar, aln_simd);
}

// F > H works for both, albeit it doesn't match the order of the scalar implementation.
// F==H seemed to yield this monstrosity. Methinks to fix it you'd need to detect stops again.
// E.g., F==H && !stopped. This is extra cycles though.

#[test]
fn sw_simd_aln9() {
    let reference = b"TTTTTGTTTTCTTGGT";
    let query = b"TTTTTTTCTTGTTTTTG";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<u8, 16, 5>::new(query, &weights.into_biased_matrix(), GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));

    let recaculated_score = sw_score_from_path(&aln_simd.states, &reference[aln_simd.ref_range.clone()], &profile_scalar);
    assert_eq!(
        Ok(score.unwrap()),
        recaculated_score,
        "{aln_simd:?} (v {cigar})",
        cigar = aln_scalar.states
    );
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln10() {
    let reference = b"TTTTTTTTTTTTAAAATTTGTAAACGTTTTGTTA";
    let query = b"TTTTTTTTACTATTTTTAAATTTATGTTTTGTTA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<u8, 8, 5>::new(query, &weights.into_biased_matrix(), GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));

    assert_eq!(aln_scalar.ref_range.end, aln_simd.ref_range.end);
    assert_eq!(aln_scalar.query_range.end, aln_simd.query_range.end);
    let recaculated_score = sw_score_from_path(&aln_simd.states, &reference[aln_simd.ref_range.clone()], &profile_scalar);
    assert_eq!(
        Ok(score.unwrap()),
        recaculated_score,
        "{aln_simd:?} (v {cigar})",
        cigar = aln_scalar.states
    );
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln11() {
    let reference = b"TTTATTTTTTTTTTTTTTCCCCCCCTTTTTTTTTTTTTTTTTCCCCCCTTT";
    let query = b"TTTTTTTTTTTTTTTTTTTCCTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCTTTA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<u8, 8, 5>::new(query, &weights.into_biased_matrix(), GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));

    assert_eq!(aln_scalar.ref_range.end, aln_simd.ref_range.end);
    assert_eq!(aln_scalar.query_range.end, aln_simd.query_range.end);
    let recaculated_score = sw_score_from_path(&aln_simd.states, &reference[aln_simd.ref_range.clone()], &profile_scalar);
    assert_eq!(
        Ok(score.unwrap()),
        recaculated_score,
        "{aln_simd:?} (v {cigar})",
        cigar = aln_scalar.states
    );
    assert_eq!(aln_scalar, aln_simd);
}

#[test]
fn sw_simd_aln12() {
    let reference = b"TTTTTTTTTTTTTTTCCCCCTTTTTTTTTTCCCCCCCCCTT";
    let query = b"TTTTTTTTTTTTTTTCCTTTTTTTTTTTTTTTTTTTCCCCCCCCCTA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<u8, 8, 5>::new(query, &weights.into_biased_matrix(), GAP_OPEN, GAP_EXTEND).unwrap();

    let score = profile_simd.smith_waterman_score(reference);

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap().as_u64();
    let aln_simd = profile_simd.smith_waterman_alignment(reference).unwrap();
    assert_eq!(score, Some(aln_scalar.score));

    assert_eq!(aln_scalar.ref_range.end, aln_simd.ref_range.end);
    assert_eq!(aln_scalar.query_range.end, aln_simd.query_range.end);
    let recaculated_score = sw_score_from_path(&aln_simd.states, &reference[aln_simd.ref_range.clone()], &profile_scalar);
    assert_eq!(
        Ok(score.unwrap()),
        recaculated_score,
        "{aln_simd:?} (v {cigar})",
        cigar = aln_scalar.states
    );
    assert_eq!(aln_scalar, aln_simd);
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

    let profile = LocalProfiles::new_with_w128(&v, &matrix_i, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = profile.smith_waterman_score_from_i8(&v);
    assert_eq!(Some(3372), score);
}
