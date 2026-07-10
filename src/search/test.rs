use crate::search::fuzzy_substring_match;

#[test]
fn fuzzy_matching_no_diff() {
    let haystack = b"ATAGATTTCTGCTACTCTCCTCATAAGCAGTCCGGTGTATCGAAAGTACA";
    let needle = b"TCCTCATAAG";
    assert_eq!(fuzzy_substring_match(haystack, needle, 0), Some(17));
}

#[test]
fn fuzzy_matching_one_diff() {
    let haystack = b"ATAGATTTCTGCTACTCTCCTCATAAGCAGTCCGGTGTATCGAAAGTACA";
    let needle = b"TCCTTATAAG";
    assert_eq!(fuzzy_substring_match(haystack, needle, 1), Some(17));
}

#[test]
fn fuzzy_matching_two_diff() {
    let haystack = b"ATAGATTTCTGCTACTCTCCTCATAAGCAGTCCGGTGTATCGAAAGTACA";
    let needle = b"TCATTATAAG";
    assert_eq!(fuzzy_substring_match(haystack, needle, 2), Some(17));
}

#[test]
fn fuzzy_matching_long_needle() {
    let haystack = b"ACGT";
    let needle = b"GACTAG";
    assert_eq!(fuzzy_substring_match(haystack, needle, 1), None);
}

#[test]
fn fuzzy_matching_no_haystack() {
    let haystack = b"";
    let needle = b"ACGT";
    assert_eq!(fuzzy_substring_match(haystack, needle, 1), None);
}

#[test]
fn fuzzy_matching_no_needle() {
    let haystack = b"GACTAG";
    let needle = b"";
    assert_eq!(fuzzy_substring_match(haystack, needle, 1), None);
}

#[test]
fn fuzzy_matching_many_diffs() {
    let haystack = b"NNNNNN";
    let needle = b"AC";
    assert_eq!(fuzzy_substring_match(haystack, needle, 2), Some(0));
}

#[test]
fn fuzzy_matching_u8_max() {
    let haystack = vec![b'A'; 256];
    let needle = vec![b'B'; 256];
    assert_eq!(fuzzy_substring_match(&haystack, &needle, u8::MAX), None);
}

#[test]
fn fuzzy_matching_same_len() {
    let haystack = b"TCCTCATAAG";
    let needle = b"TCCTCTTAAG";
    assert_eq!(fuzzy_substring_match(haystack, needle, 1), Some(0));
}

#[test]
fn next_aa_range_search_short_tail() {
    use crate::{prelude::Nucleotides, search::ToRangeSearch};

    let seq: Nucleotides = b"AA".as_slice().into();
    assert_eq!(seq.search_in(1..2).find_next_aa_in_frame(b'M'), None);
}
