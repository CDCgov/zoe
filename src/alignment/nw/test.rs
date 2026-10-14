use crate::alignment::{
    Alignment, AlignmentStates, ScalarProfile, ScoreParams,
    nw::{nw_scalar_align, nw_scalar_score},
};

#[test]
fn nw_empty_ref() {
    let query = b"ACGT";
    let reference = b"";

    let params = ScoreParams::new_dna_params(2, -5, -10, -1, Some(b'N')).unwrap();

    let profile = ScalarProfile::new(query, &params);

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

    let params = ScoreParams::new_dna_params(2, -5, -10, -1, Some(b'N')).unwrap();

    let profile = ScalarProfile::new(query, &params);

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

    let params = ScoreParams::new_dna_params(2, -5, -10, -1, Some(b'N')).unwrap();

    let profile = ScalarProfile::new(query, &params);

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
