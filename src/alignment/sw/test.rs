use super::{test_data::*, *};
use crate::alignment::profile::StripedProfile;
use crate::data::mappings::DNA_PROFILE_MAP;
use crate::data::matrices::WeightMatrix;

macro_rules! test_sw_simd_alignment {
    ($profile_seq:expr, $other_seq:expr, $int_type:ty, $uint_type:ty, $lanes:expr) => {{
        let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
        let profile_scalar = ScalarProfile::new($profile_seq, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
        let profile_simd =
            StripedProfile::<$int_type, $lanes, 5>::new($profile_seq, &weights, GAP_OPEN, GAP_EXTEND).unwrap();

        let score = profile_simd.smith_waterman_score($other_seq);
        let aln_scalar = profile_scalar.smith_waterman_alignment($other_seq).unwrap();
        let aln_simd = profile_simd.smith_waterman_alignment($other_seq).unwrap();

        assert_eq!(score, MaybeAligned::Some(aln_scalar.score));
        assert_eq!(aln_scalar, aln_simd);

        let profile_simd =
            StripedProfile::<$uint_type, $lanes, 5>::new($profile_seq, &weights.to_biased_matrix(), GAP_OPEN, GAP_EXTEND)
                .unwrap();
        let aln_simd = profile_simd.smith_waterman_alignment($other_seq).unwrap();
        assert_eq!(aln_scalar, aln_simd, "SIMD ALIGNMENT");

        let (score, ref_end, query_end) = profile_simd.smith_waterman_score_ends($other_seq).unwrap();
        assert_eq!(score, aln_scalar.score);
        assert_eq!(aln_scalar.ref_range.end, ref_end, "REFERENCE END");
        assert_eq!(aln_scalar.query_range.end, query_end, "QUERY END");
    }};
}

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

    assert_eq!(Ok(score), sw_score_from_path(&states, &reference[ref_range], &profile));
}

#[test]
fn sw() {
    let weights: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
    let profile = ScalarProfile::new(QUERY, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let Alignment { score, ref_range, .. } = sw_scalar_alignment(REFERENCE, &profile).unwrap();
    assert_eq!((336, 37), (ref_range.start, score));

    let score = sw_scalar_score(REFERENCE, &profile).unwrap();
    assert_eq!(37, score);

    let v: Vec<_> = std::iter::repeat_n(b'A', 100).collect();
    let profile = ScalarProfile::new(&v, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = sw_scalar_score(&v, &profile).unwrap();
    assert_eq!(200, score);
}

#[test]
fn sw_t_u_check() {
    let profile = ScalarProfile::new(b"ACGTUNacgtun", &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = sw_scalar_score(b"ACGTTNACGTTN", &profile);
    assert_eq!(score, MaybeAligned::Some(20));

    let profile = StripedProfile::<u16, 16, 5>::new(b"ACGTUNacgtun", &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = sw_simd_score(b"ACGTTNACGTTN", &profile);
    assert_eq!(score, MaybeAligned::Some(20));

    let profile = StripedProfile::<i16, 16, 5>::new(b"ACGTUNacgtun", &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = sw_simd_score(b"ACGTTNACGTTN", &profile);
    assert_eq!(score, MaybeAligned::Some(20));
}

#[test]
fn sw_simd() {
    let matrix_i = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let matrix_u = matrix_i.to_biased_matrix();

    let profile = StripedProfile::<u8, 16, 5>::new(REFERENCE, &matrix_u, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = profile.smith_waterman_score(QUERY);
    assert_eq!(MaybeAligned::Some(37), score);

    let profile = StripedProfile::<i16, 16, 5>::new(REFERENCE, &matrix_i, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = profile.smith_waterman_score(QUERY);
    assert_eq!(MaybeAligned::Some(37), score);
}

#[test]
fn sw_simd_aln() {
    test_sw_simd_alignment!(REFERENCE, QUERY, i16, u16, 8);
}

#[test]
fn sw_simd_aln2() {
    let reference = b"AAACTA";
    let query = b"TTTAG";
    test_sw_simd_alignment!(query, reference, i8, u8, 2);
}

#[test]
fn sw_simd_aln3() {
    let reference = b"AAAAAAAAAA";
    let query = b"AAAAAATAAA";
    test_sw_simd_alignment!(query, reference, i8, u8, 4);
}

#[test]
fn sw_simd_aln4() {
    let reference = b"TAAAA";
    let query = b"CCCCA";
    test_sw_simd_alignment!(query, reference, i8, u8, 4);
}

#[test]
fn sw_simd_aln5() {
    let reference = b"TCCCC";
    let query = b"CCCCC";
    test_sw_simd_alignment!(query, reference, i8, u8, 4);
}

#[test]
fn sw_simd_aln6() {
    let reference = b"GCTTTTC";
    let query = b"CCCCT";
    test_sw_simd_alignment!(query, reference, i8, u8, 4);
}
#[test]
fn sw_simd_aln7() {
    let reference = b"TTGTTTTTTTTTGTT";
    let query = b"TTTTTGTTTTCTTTTTTGTTTA";
    test_sw_simd_alignment!(query, reference, i8, u8, 16);
}

#[test]
fn sw_simd_aln8() {
    let reference = b"TTTTTGTTTGGGAAAAATTCTT";
    let query = b"TTGTTTTGGGGAAAAA";
    test_sw_simd_alignment!(query, reference, i8, u8, 8);
}

#[test]
fn sw_simd_aln9() {
    let reference = b"TTTTTGTTTTCTTGGT";
    let query = b"TTTTTTTCTTGTTTTTG";
    test_sw_simd_alignment!(query, reference, i8, u8, 16);
}

#[test]
fn sw_simd_aln10() {
    let reference = b"TTTTTTTTTTTTAAAATTTGTAAACGTTTTGTTA";
    let query = b"TTTTTTTTACTATTTTTAAATTTATGTTTTGTTA";
    test_sw_simd_alignment!(query, reference, i8, u8, 8);
}

#[test]
fn sw_simd_aln11() {
    let reference = b"TTTATTTTTTTTTTTTTTCCCCCCCTTTTTTTTTTTTTTTTTCCCCCCTTT";
    let query = b"TTTTTTTTTTTTTTTTTTTCCTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCTTTA";
    test_sw_simd_alignment!(query, reference, i8, u8, 8);
}

#[test]
fn sw_simd_aln12() {
    let reference = b"TTTTTTTTTTTTTTTCCCCCTTTTTTTTTTCCCCCCCCCTT";
    let query = b"TTTTTTTTTTTTTTTCCTTTTTTTTTTTTTTTTTTTCCCCCCCCCTA";
    test_sw_simd_alignment!(query, reference, i8, u8, 8);
}

#[cfg(feature = "dev-3pass")]
#[test]
fn sw_simd_locations() {
    let reference = b"TTTTTTCCTTTTTTTTCCCCCTTTTT";
    let query = b"GGGGGGGCCCCCAAAA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<u8, 8, 5>::new(query, &weights.to_biased_matrix(), GAP_OPEN, GAP_EXTEND).unwrap();

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap();
    let (score, ref_end, query_end) = profile_simd.smith_waterman_score_ends(reference).unwrap();

    assert_eq!(score, aln_scalar.score);
    assert_eq!(aln_scalar.ref_range.end, ref_end, "REFERENCE END");
    assert_eq!(aln_scalar.query_range.end, query_end, "QUERY END");

    let query_rev: Vec<u8> = query[..query_end].iter().copied().rev().collect();

    let profile_simd_rev =
        StripedProfile::<u8, 8, 5>::new(&query_rev, &weights.to_biased_matrix(), GAP_OPEN, GAP_EXTEND).unwrap();

    let (score2, ref_start, query_start) = profile_simd_rev
        .smith_waterman_score_ends_reverse(&reference[..ref_end])
        .unwrap();
    assert_eq!(score, score2);
    assert_eq!(ref_start..ref_end, aln_scalar.ref_range, "REFERENCE");
    assert_eq!(query_start..query_end, aln_scalar.query_range, "QUERY");

    let profile_simd_rev2 = profile_simd.reverse_from_forward(query_end).unwrap();
    assert_eq!(profile_simd_rev, profile_simd_rev2);

    let (_, ref_range, query_range) = profile_simd.smith_waterman_score_ranges(reference).unwrap();
    assert_eq!(query_range, aln_scalar.query_range, "3pass QUERY RANGE");
    assert_eq!(ref_range, aln_scalar.ref_range, "3pass REFERENCE RANGE");
}

#[cfg(feature = "dev-3pass")]
#[test]
fn sw_simd_ranges() {
    let reference = b"TAAAA";
    let query = b"CCCCA";

    let weights = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));
    let profile_scalar = ScalarProfile::new(query, &weights, GAP_OPEN, GAP_EXTEND).unwrap();
    let profile_simd = StripedProfile::<u8, 8, 5>::new(query, &weights.to_biased_matrix(), GAP_OPEN, GAP_EXTEND).unwrap();

    let aln_scalar = profile_scalar.smith_waterman_alignment(reference).unwrap();
    let (score, ref_range, query_range) = profile_simd.smith_waterman_score_ranges(reference).unwrap();

    assert_eq!(score, aln_scalar.score);
    assert_eq!(ref_range, aln_scalar.ref_range, "REFERENCE");
    assert_eq!(query_range, aln_scalar.query_range, "QUERY");
}

#[test]
fn sw_simd_poly_a() {
    let v: Vec<_> = std::iter::repeat_n(b'A', 100).collect();
    let matrix = WeightMatrix::new_dna_matrix(2, -5, Some(b'N')).to_biased_matrix();
    let profile = StripedProfile::<u16, 16, 5>::new(&v, &matrix, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = profile.smith_waterman_score(&v);
    assert_eq!(MaybeAligned::Some(200), score);
}

#[test]
fn sw_simd_single() {
    let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
    let matrix = WeightMatrix::<i8, 5>::new_dna_matrix(2, -5, Some(b'N')).to_biased_matrix();
    let profile = StripedProfile::<u16, 16, 5>::new(v, &matrix, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = profile.smith_waterman_score(v);
    assert_eq!(MaybeAligned::Some(3372), score);
}

#[test]
fn sw_simd_regression() {
    let query = b"AGA";
    let reference = b"AA";
    let matrix = WeightMatrix::new(&DNA_PROFILE_MAP, 10, -10, Some(b'N')).to_biased_matrix();
    let profile = StripedProfile::<u16, 4, 5>::new(query, &matrix, -5, -5).unwrap();
    let score = profile.smith_waterman_score(reference);
    assert_eq!(MaybeAligned::Some(15), score);
}

#[test]
fn sw_simd_overflow_check() {
    let query = b"AAAA";
    let reference = b"AAAA";

    let matrix = WeightMatrix::new(&DNA_PROFILE_MAP, 127, 0, Some(b'N')).to_biased_matrix();
    let profile = StripedProfile::<u8, 8, 5>::new(query, &matrix, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = profile.smith_waterman_score(reference);
    assert_eq!(score, MaybeAligned::Overflowed);
}

#[test]
fn sw_simd_profile_set() {
    let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
    let matrix_i = WeightMatrix::<i8, 5>::new_dna_matrix(2, -5, Some(b'N'));

    let profile = LocalProfiles::new_with_w128(&v, &matrix_i, GAP_OPEN, GAP_EXTEND).unwrap();
    let score = profile.smith_waterman_score_from_i8(&v);
    assert_eq!(MaybeAligned::Some(3372), score);
}
