use crate::assert_fp_eq;

use super::*;

pub(crate) static LONG_READ_1: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/KJ907623.1.txt")); // H5 HA, complete CDS
pub(crate) static LONG_READ_2: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/KJ907631.1.txt")); // H5 HA, complete CDS

static AMBIGUOUS_READ_1: &[u8] = b"ACGTX";
static AMBIGUOUS_READ_2: &[u8] = b"ACGAT";
static SHORT_READ_1: &[u8] = b"ATCGATCGATCG";
static SHORT_READ_2: &[u8] = b"ATCGATCGATCA";
static MID_READ_GAPS_1: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/KP326324.1.txt")); // H5 HA, complete CDS, aligned
static MID_READ_GAPS_2: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/KC684446.1.txt")); // H5 HA, complete CDS, aligned

#[test]
fn test_hamming_from_matrix() {
    let sub_matrix = dna_substitution_matrix(SHORT_READ_1, SHORT_READ_2);
    let hamming = hamming_dist_from_sub_matrix(&sub_matrix);
    assert_eq!(hamming, 1);

    let sub_matrix = dna_substitution_matrix(MID_READ_GAPS_1, MID_READ_GAPS_2);
    let hamming = hamming_dist_from_sub_matrix(&sub_matrix);
    assert_eq!(hamming, 9);

    let sub_matrix = dna_substitution_matrix(LONG_READ_1, LONG_READ_2);
    let hamming = hamming_dist_from_sub_matrix(&sub_matrix);
    assert_eq!(hamming, 148);
}

#[test]
fn test_jc() {
    let data = [
        (SHORT_READ_1, SHORT_READ_2, Some(0.088_337_276_7)),
        (MID_READ_GAPS_1, MID_READ_GAPS_2, Some(0.040_188_184_6)),
        (LONG_READ_1, LONG_READ_2, Some(0.092_309_042_4)),
        (LONG_READ_1, LONG_READ_1, Some(0.0)),
        (LONG_READ_1, b"", None),
        (&LONG_READ_1[0..20], LONG_READ_2, Some(0.051_744_653_615_213_576)),
        (
            AMBIGUOUS_READ_1,
            AMBIGUOUS_READ_2,
            dna_substitution_matrix(b"ACGT", b"ACGA").jc69_distance(),
        ),
    ];

    for (a, b, expected) in data {
        let sm = dna_substitution_matrix(a, b);
        assert_fp_eq!(sm.jc69_distance(), expected);
        assert_fp_eq!(jukes_cantor_69(a, b), expected);
    }
}

#[test]
fn compare_jc_regression() {
    let data: Vec<(&[u8], &[u8])> = vec![
        (b"gug", b"ggg"),
        (b"gggcaaaccacaaaaaaaaaaaaaaaaaaaaaa", b"aaaaaaaaaaaaaaaaaaaaaataatttttttt"),
    ];
    for (a, b) in data.into_iter().skip(1) {
        let sm = dna_substitution_matrix(a, b);
        assert_eq!(sm.jc69_distance(), jukes_cantor_69(a, b));
    }
}

#[test]
fn test_k80() {
    let data = [
        (SHORT_READ_1, SHORT_READ_2, Some(0.091_160_778_3)),
        (MID_READ_GAPS_1, MID_READ_GAPS_2, Some(0.040_576_991_6)),
        (LONG_READ_1, LONG_READ_2, Some(0.093_274_373_1)),
        (LONG_READ_1, LONG_READ_1, Some(0.0)),
        (b"", b"", None),
        (LONG_READ_1, b"", None),
        (AMBIGUOUS_READ_1, AMBIGUOUS_READ_2, kimura_80(b"ACGT", b"ACGA")),
    ];

    for (a, b, expected) in data {
        assert_fp_eq!(kimura_80(a, b), expected);
    }
}

#[test]
fn test_k81() {
    let data = [
        (SHORT_READ_1, SHORT_READ_2, Some(0.091_160_778_396_977_3)),
        (MID_READ_GAPS_1, MID_READ_GAPS_2, Some(0.040_582_502_030_320_3)),
        (LONG_READ_1, LONG_READ_2, Some(0.093_297_615_589_788_4)),
        (LONG_READ_1, LONG_READ_1, Some(0.0)),
        (LONG_READ_1, b"", None),
        (b"", b"", None),
        (AMBIGUOUS_READ_1, AMBIGUOUS_READ_2, kimura_81(b"ACGT", b"ACGA")),
    ];
    for (a, b, expected) in data {
        assert_fp_eq!(kimura_81(a, b), expected);
    }
}

#[test]
fn test_f81() {
    let data = [
        (SHORT_READ_1, SHORT_READ_2, Some(0.088_362_461_866_048_9)),
        (MID_READ_GAPS_1, MID_READ_GAPS_2, Some(0.040_201_884_674_593_0)),
        (LONG_READ_1, LONG_READ_2, Some(0.092_428_092_165_986_7)),
        (LONG_READ_1, LONG_READ_1, Some(0.0)),
        (LONG_READ_1, b"", None),
        (b"", b"", None),
        (AMBIGUOUS_READ_1, AMBIGUOUS_READ_2, felsenstein_81(b"ACGT", b"ACGA")),
    ];

    for (a, b, expected) in data {
        assert_fp_eq!(felsenstein_81(a, b), expected);
    }
}

#[test]
fn test_tn93() {
    let data = [
        (SHORT_READ_1, SHORT_READ_2, Some(0.102_047_809_6)),
        (MID_READ_GAPS_1, MID_READ_GAPS_2, Some(0.040_633_116_7)),
        (LONG_READ_1, LONG_READ_2, Some(0.093_467_734_2)),
        (LONG_READ_1, LONG_READ_1, Some(0.0)),
        (LONG_READ_1, b"", None),
        (b"", b"", None),
        (AMBIGUOUS_READ_1, AMBIGUOUS_READ_2, tamura_nei_93(b"ACGT", b"ACGA")),
    ];

    for (a, b, expected) in data {
        assert_fp_eq!(tamura_nei_93(a, b), expected);
    }
}
