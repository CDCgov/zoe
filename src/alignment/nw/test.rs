use crate::{
    alignment::{
        Alignment, AlignmentStates, ScalarProfile,
        nw::{nw_scalar_align, nw_scalar_score},
    },
    data::WeightMatrix,
};

#[test]
fn nw_empty_ref() {
    let query = b"ACGT";
    let reference = b"";

    let weights: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
    let gap_open = -10;
    let gap_extend = -1;

    let profile = ScalarProfile::new(query, &weights, gap_open, gap_extend).unwrap();

    let score = nw_scalar_score(reference, &profile);
    assert_eq!(score, -13);

    let Alignment {
        score,
        ref_range,
        query_range,
        states,
        ref_len,
        query_len,
    } = nw_scalar_align(reference, &profile);

    assert_eq!(score, -13);
    assert_eq!(ref_range, 0..reference.len());
    assert_eq!(query_range, 0..query.len());
    assert_eq!(states, AlignmentStates::try_from(b"4I").unwrap());
    assert_eq!(ref_len, reference.len());
    assert_eq!(query_len, query.len());
}

#[test]
fn nw_empty_query() {
    let query = b"";
    let reference = b"ACGT";

    let weights: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
    let gap_open = -10;
    let gap_extend = -1;

    let profile = ScalarProfile::new(query, &weights, gap_open, gap_extend).unwrap();

    let score = nw_scalar_score(reference, &profile);
    assert_eq!(score, -13);

    let Alignment {
        score,
        ref_range,
        query_range,
        states,
        ref_len,
        query_len,
    } = nw_scalar_align(reference, &profile);

    assert_eq!(score, -13);
    assert_eq!(ref_range, 0..reference.len());
    assert_eq!(query_range, 0..query.len());
    assert_eq!(states, AlignmentStates::try_from(b"4D").unwrap());
    assert_eq!(ref_len, reference.len());
    assert_eq!(query_len, query.len());
}

#[test]
fn nw_empty_both() {
    let query = b"";
    let reference = b"";

    let weights: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
    let gap_open = -10;
    let gap_extend = -1;

    let profile = ScalarProfile::new(query, &weights, gap_open, gap_extend).unwrap();

    let score = nw_scalar_score(reference, &profile);
    assert_eq!(score, 0);

    let Alignment {
        score,
        ref_range,
        query_range,
        states,
        ref_len,
        query_len,
    } = nw_scalar_align(reference, &profile);

    assert_eq!(score, 0);
    assert_eq!(ref_range, 0..reference.len());
    assert_eq!(query_range, 0..query.len());
    assert_eq!(states, AlignmentStates::new());
    assert_eq!(ref_len, reference.len());
    assert_eq!(query_len, query.len());
}
