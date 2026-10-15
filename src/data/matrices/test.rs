use super::*;
use crate::data::{ByteIndexMap, DNA_PROFILE_MAP};

static AA_MATS: [&WeightMatrix<'static, i8, 25>; 22] = [
    &BLOSUM_30,
    &BLOSUM_35,
    &BLOSUM_40,
    &BLOSUM_45,
    &BLOSUM_50,
    &BLOSUM_55,
    &BLOSUM_60,
    &BLOSUM_62,
    &BLOSUM_65,
    &BLOSUM_70,
    &BLOSUM_75,
    &BLOSUM_80,
    &BLOSUM_85,
    &BLOSUM_90,
    &BLOSUM_95,
    &BLOSUM_100,
    &PAM_30,
    &PAM_40,
    &PAM_70,
    &PAM_120,
    &PAM_200,
    &PAM_250,
];

#[test]
fn create_simple() {
    static RESIDUE_MAP: ByteIndexMap<2> = ByteIndexMap::new(*b"AB", b'A');

    let result1 = WeightMatrix {
        weights: [[1, 0], [0, 1]],
        mapping: &RESIDUE_MAP,
        bias:    0,
    };

    let result2 = WeightMatrix::new(&RESIDUE_MAP, 1, 0, None);

    let result3 = WeightMatrix::new_custom(&RESIDUE_MAP, [[1, 0], [0, 1]]);

    assert_eq!(result1, result2);
    assert_eq!(result1, result3);
}

#[test]
#[should_panic(expected = "An invalid byte was specified for the ignoring field.")]
fn test_invalid_char() {
    let _ = WeightMatrix::new_dna_matrix(1, 0, Some(b'U'));
}

#[allow(non_snake_case)]
#[test]
fn create_IRMA_matrix() {
    let result1 = WeightMatrix {
        weights: [
            [2, -5, -5, -5, 0],
            [-5, 2, -5, -5, 0],
            [-5, -5, 2, -5, 0],
            [-5, -5, -5, 2, 0],
            [0, 0, 0, 0, 0],
        ],
        mapping: &DNA_PROFILE_MAP,
        bias:    0,
    };

    let result2 = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));

    let result3 = WeightMatrix::new_custom(
        &DNA_PROFILE_MAP,
        [
            [2, -5, -5, -5, 0],
            [-5, 2, -5, -5, 0],
            [-5, -5, 2, -5, 0],
            [-5, -5, -5, 2, 0],
            [0, 0, 0, 0, 0],
        ],
    );

    assert_eq!(result1, result2);
    assert_eq!(result1, result3);
}

#[test]
fn test_blosum_matrices() {
    assert!(aa_mat_from_name("BLOSUM30").is_some());
    assert!(aa_mat_from_name("blosum 62").is_some());
    assert!(aa_mat_from_name("BLOSUM_80").is_some());
    assert!(aa_mat_from_name("blosum-100").is_some());
    assert!(aa_mat_from_name("BLOSUM 45").is_some());

    assert!(aa_mat_from_name("BLOSUM29").is_none());
    assert!(aa_mat_from_name("BLOSUM101").is_none());
    assert!(aa_mat_from_name("BLOSUM").is_none());
}

#[test]
fn test_pam_matrices() {
    assert!(aa_mat_from_name("PAM30").is_some());
    assert!(aa_mat_from_name("pam-70").is_some());
    assert!(aa_mat_from_name("PAM_120").is_some());
    assert!(aa_mat_from_name("PAM 200").is_some());

    assert!(aa_mat_from_name("PAM39").is_none());
    assert!(aa_mat_from_name("PAM2500").is_none());
    assert!(aa_mat_from_name("PAM").is_none());
}

#[test]
fn test_invalid_names() {
    assert!(aa_mat_from_name("INVALID30").is_none());
    assert!(aa_mat_from_name("UNKNOWN").is_none());
    assert!(aa_mat_from_name("PAM_").is_none());
    assert!(aa_mat_from_name("").is_none());
}

#[test]
fn test_subset_unbiased() {
    const OTHER_MAP: ByteIndexMap<21> = ByteIndexMap::new_ignoring_case(*b"DIEVRSAWGMCYTXNPFLHKQ", b'X');

    for mat in AA_MATS {
        let other_mat = mat.get_subset(&OTHER_MAP);
        for aa1 in u8::MIN..=u8::MAX {
            for aa2 in u8::MIN..=u8::MAX {
                let other_aa1 = if OTHER_MAP.byte_keys().contains(&aa1.to_ascii_uppercase()) {
                    aa1
                } else {
                    b'X'
                };
                let other_aa2 = if OTHER_MAP.byte_keys().contains(&aa2.to_ascii_uppercase()) {
                    aa2
                } else {
                    b'X'
                };

                let score = mat.get_weight(other_aa1, other_aa2);
                let other_score = other_mat.get_weight(aa1, aa2);
                assert_eq!(score, other_score, "{aa1} vs {aa2}");
            }
        }
    }
}

#[test]
#[allow(clippy::cast_possible_wrap)]
fn test_subset_biased() {
    const OTHER_MAP: ByteIndexMap<4> = ByteIndexMap::new_ignoring_case(*b"AWGX", b'X');

    for mat in AA_MATS.map(WeightMatrix::to_biased_matrix) {
        let other_mat = mat.get_subset(&OTHER_MAP);
        for aa1 in u8::MIN..=u8::MAX {
            for aa2 in u8::MIN..=u8::MAX {
                let other_aa1 = if OTHER_MAP.byte_keys().contains(&aa1.to_ascii_uppercase()) {
                    aa1
                } else {
                    b'X'
                };
                let other_aa2 = if OTHER_MAP.byte_keys().contains(&aa2.to_ascii_uppercase()) {
                    aa2
                } else {
                    b'X'
                };

                let score = i8::try_from(mat.get_weight(other_aa1, other_aa2)).unwrap() - i8::try_from(mat.bias).unwrap();
                let other_score =
                    i8::try_from(other_mat.get_weight(aa1, aa2)).unwrap() - i8::try_from(other_mat.bias).unwrap();
                assert_eq!(score, other_score, "{aa1} vs {aa2}");
            }
        }
    }
}

#[test]
#[allow(clippy::cast_possible_wrap)]
fn test_subset_correct_bias() {
    for byte_keys in generate_subsets::<u8, 14, 4>(b"ADFHIKMQTVY*BZ") {
        let mapping = ByteIndexMap::new_ignoring_case([byte_keys[0], byte_keys[1], byte_keys[2], byte_keys[3], b'X'], b'X');
        for mat in AA_MATS.map(WeightMatrix::to_biased_matrix) {
            let other_mat = mat.get_subset(&mapping);
            let min_weight = *other_mat.weights.map(|row| *row.iter().min().unwrap()).iter().min().unwrap();
            assert_eq!(min_weight, 0);
        }
    }
}

fn generate_subsets<T: Copy + Default, const S: usize, const L: usize>(arr: &[T; S]) -> Vec<[T; L]> {
    fn generate_subsets_inner<T: Copy, const S: usize, const L: usize>(
        arr: &[T; S], start: usize, current: [T; L], current_len: usize, out: &mut Vec<[T; L]>,
    ) {
        if current_len == L {
            out.push(current);
            return;
        }

        for i in start..S {
            let mut new_current = current;
            new_current[current_len] = arr[i];
            generate_subsets_inner::<T, S, L>(arr, i + 1, new_current, current_len + 1, out);
        }
    }

    let mut out = Vec::new();
    generate_subsets_inner::<T, S, L>(arr, 0, [T::default(); L], 0, &mut out);
    out
}
