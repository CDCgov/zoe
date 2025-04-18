use super::{
    ByteIndexMap, DNA_PROFILE_MAP, IS_DNA_ACGTN_NO_GAPS, IS_DNA_ACGTN_NO_GAPS_UC, IS_DNA_IUPAC_NO_GAPS,
    IS_DNA_IUPAC_NO_GAPS_UC, TO_DNA_IUPAC_NO_GAPS_UC, TO_DNA_IUPAC_WITH_GAPS_UC,
};
use crate::{
    composition::ByteIndexCounts,
    data::{
        alphas::{DNA_ACGTN_NO_GAPS, DNA_ACGTN_NO_GAPS_UC, DNA_IUPAC_NO_GAPS, DNA_IUPAC_NO_GAPS_UC, DNA_IUPAC_WITH_GAPS},
        mappings::{
            ANY_TO_DNA_ACGTN_NO_GAPS_UC, ANY_TO_DNA_IUPAC_WITH_GAPS, ANY_TO_DNA_IUPAC_WITH_GAPS_UC, DNA_UNAMBIG_PROFILE_MAP,
            IS_DNA_IUPAC_WITH_GAPS, IUPAC_TO_DNA_ACGTN_WITH_GAPS, IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC, StdGeneticCode,
        },
    },
    prelude::Nucleotides,
};

#[test]
fn gc_self_test() {
    #[rustfmt::skip]
        let gc = [
            (b"TAA", b'*'), (b"TAG", b'*'), (b"TAR", b'*'), (b"TGA", b'*'), (b"TRA", b'*'), (b"GCA", b'A'), (b"GCB", b'A'), (b"GCC", b'A'), (b"GCD", b'A'), (b"GCG", b'A'), (b"GCH", b'A'),
            (b"GCK", b'A'), (b"GCM", b'A'), (b"GCN", b'A'), (b"GCR", b'A'), (b"GCS", b'A'), (b"GCT", b'A'), (b"GCV", b'A'), (b"GCW", b'A'), (b"GCY", b'A'), (b"TGC", b'C'), (b"TGT", b'C'),
            (b"TGY", b'C'), (b"GAC", b'D'), (b"GAT", b'D'), (b"GAY", b'D'), (b"GAA", b'E'), (b"GAG", b'E'), (b"GAR", b'E'), (b"TTC", b'F'), (b"TTT", b'F'), (b"TTY", b'F'), (b"GGA", b'G'),
            (b"GGB", b'G'), (b"GGC", b'G'), (b"GGD", b'G'), (b"GGG", b'G'), (b"GGH", b'G'), (b"GGK", b'G'), (b"GGM", b'G'), (b"GGN", b'G'), (b"GGR", b'G'), (b"GGS", b'G'), (b"GGT", b'G'),
            (b"GGV", b'G'), (b"GGW", b'G'), (b"GGY", b'G'), (b"CAC", b'H'), (b"CAT", b'H'), (b"CAY", b'H'), (b"ATA", b'I'), (b"ATC", b'I'), (b"ATH", b'I'), (b"ATM", b'I'), (b"ATT", b'I'),
            (b"ATW", b'I'), (b"ATY", b'I'), (b"AAA", b'K'), (b"AAG", b'K'), (b"AAR", b'K'), (b"CTA", b'L'), (b"CTB", b'L'), (b"CTC", b'L'), (b"CTD", b'L'), (b"CTG", b'L'), (b"CTH", b'L'),
            (b"CTK", b'L'), (b"CTM", b'L'), (b"CTN", b'L'), (b"CTR", b'L'), (b"CTS", b'L'), (b"CTT", b'L'), (b"CTV", b'L'), (b"CTW", b'L'), (b"CTY", b'L'), (b"TTA", b'L'), (b"TTG", b'L'),
            (b"TTR", b'L'), (b"YTA", b'L'), (b"YTG", b'L'), (b"YTR", b'L'), (b"ATG", b'M'), (b"AAC", b'N'), (b"AAT", b'N'), (b"AAY", b'N'), (b"CCA", b'P'), (b"CCB", b'P'), (b"CCC", b'P'),
            (b"CCD", b'P'), (b"CCG", b'P'), (b"CCH", b'P'), (b"CCK", b'P'), (b"CCM", b'P'), (b"CCN", b'P'), (b"CCR", b'P'), (b"CCS", b'P'), (b"CCT", b'P'), (b"CCV", b'P'), (b"CCW", b'P'),
            (b"CCY", b'P'), (b"CAA", b'Q'), (b"CAG", b'Q'), (b"CAR", b'Q'), (b"AGA", b'R'), (b"AGG", b'R'), (b"AGR", b'R'), (b"CGA", b'R'), (b"CGB", b'R'), (b"CGC", b'R'), (b"CGD", b'R'),
            (b"CGG", b'R'), (b"CGH", b'R'), (b"CGK", b'R'), (b"CGM", b'R'), (b"CGN", b'R'), (b"CGR", b'R'), (b"CGS", b'R'), (b"CGT", b'R'), (b"CGV", b'R'), (b"CGW", b'R'), (b"CGY", b'R'),
            (b"MGA", b'R'), (b"MGG", b'R'), (b"MGR", b'R'), (b"AGC", b'S'), (b"AGT", b'S'), (b"AGY", b'S'), (b"TCA", b'S'), (b"TCB", b'S'), (b"TCC", b'S'), (b"TCD", b'S'), (b"TCG", b'S'),
            (b"TCH", b'S'), (b"TCK", b'S'), (b"TCM", b'S'), (b"TCN", b'S'), (b"TCR", b'S'), (b"TCS", b'S'), (b"TCT", b'S'), (b"TCV", b'S'), (b"TCW", b'S'), (b"TCY", b'S'), (b"ACA", b'T'),
            (b"ACB", b'T'), (b"ACC", b'T'), (b"ACD", b'T'), (b"ACG", b'T'), (b"ACH", b'T'), (b"ACK", b'T'), (b"ACM", b'T'), (b"ACN", b'T'), (b"ACR", b'T'), (b"ACS", b'T'), (b"ACT", b'T'),
            (b"ACV", b'T'), (b"ACW", b'T'), (b"ACY", b'T'), (b"GTA", b'V'), (b"GTB", b'V'), (b"GTC", b'V'), (b"GTD", b'V'), (b"GTG", b'V'), (b"GTH", b'V'), (b"GTK", b'V'), (b"GTM", b'V'),
            (b"GTN", b'V'), (b"GTR", b'V'), (b"GTS", b'V'), (b"GTT", b'V'), (b"GTV", b'V'), (b"GTW", b'V'), (b"GTY", b'V'), (b"TGG", b'W'), (b"TAC", b'Y'), (b"TAT", b'Y'), (b"TAY", b'Y'),
            (b"...", b'.'), (b"---", b'-'), (b"NNN", b'X'),

            (b"UAA", b'*'), (b"UAG", b'*'), (b"UAR", b'*'), (b"UGA", b'*'), (b"URA", b'*'), (b"GCU", b'A'), (b"UGC", b'C'), (b"UGU", b'C'), (b"UGY", b'C'), (b"GAU", b'D'), (b"UUC", b'F'),
            (b"UUU", b'F'), (b"UUY", b'F'), (b"GGU", b'G'), (b"CAU", b'H'), (b"AUA", b'I'), (b"AUC", b'I'), (b"AUH", b'I'), (b"AUM", b'I'), (b"AUU", b'I'), (b"AUW", b'I'), (b"AUY", b'I'),
            (b"CUA", b'L'), (b"CUB", b'L'), (b"CUC", b'L'), (b"CUD", b'L'), (b"CUG", b'L'), (b"CUH", b'L'), (b"CUK", b'L'), (b"CUM", b'L'), (b"CUN", b'L'), (b"CUR", b'L'), (b"CUS", b'L'),
            (b"CUU", b'L'), (b"CUV", b'L'), (b"CUW", b'L'), (b"CUY", b'L'), (b"UUA", b'L'), (b"UUG", b'L'), (b"UUR", b'L'), (b"YUA", b'L'), (b"YUG", b'L'), (b"YUR", b'L'), (b"AUG", b'M'),
            (b"AAU", b'N'), (b"CCU", b'P'), (b"CGU", b'R'), (b"AGU", b'S'), (b"UCA", b'S'), (b"UCB", b'S'), (b"UCC", b'S'), (b"UCD", b'S'), (b"UCG", b'S'), (b"UCH", b'S'), (b"UCK", b'S'),
            (b"UCM", b'S'), (b"UCN", b'S'), (b"UCR", b'S'), (b"UCS", b'S'), (b"UCU", b'S'), (b"UCV", b'S'), (b"UCW", b'S'), (b"UCY", b'S'), (b"ACU", b'T'), (b"GUA", b'V'), (b"GUB", b'V'),
            (b"GUC", b'V'), (b"GUD", b'V'), (b"GUG", b'V'), (b"GUH", b'V'), (b"GUK", b'V'), (b"GUM", b'V'), (b"GUN", b'V'), (b"GUR", b'V'), (b"GUS", b'V'), (b"GUU", b'V'), (b"GUV", b'V'),
            (b"GUW", b'V'), (b"GUY", b'V'), (b"UGG", b'W'), (b"UAC", b'Y'), (b"UAU", b'Y'), (b"UAY", b'Y')
        ];

    for (codon, aa) in gc {
        assert_eq!(aa, StdGeneticCode::get(codon).unwrap());
        assert_eq!(aa, StdGeneticCode::translate_codon(codon));
        assert_eq!(aa == b'*', StdGeneticCode::is_stop_codon(codon));
    }
}

#[test]
fn test_dna_map() {
    for i in 0..=255 {
        match i {
            b'A' | b'a' => assert!(DNA_PROFILE_MAP.to_index(i) == 0),
            b'C' | b'c' => assert!(DNA_PROFILE_MAP.to_index(i) == 1),
            b'G' | b'g' => assert!(DNA_PROFILE_MAP.to_index(i) == 2),
            b'T' | b't' | b'U' | b'u' => assert!(DNA_PROFILE_MAP.to_index(i) == 3),
            _ => assert!(DNA_PROFILE_MAP.to_index(i) == 4),
        }
    }
}

#[test]
fn test_dna_map_ignores_case() {
    const MAP1: ByteIndexMap<5> = ByteIndexMap::new_ignoring_case(*b"acgtn", b'N').add_synonym_ignoring_case(b'u', b'T');
    const MAP2: ByteIndexMap<5> = ByteIndexMap::new_ignoring_case(*b"AcGtN", b'n').add_synonym_ignoring_case(b'U', b't');

    assert_eq!(DNA_PROFILE_MAP, MAP1);
    assert_eq!(DNA_PROFILE_MAP, MAP2);
}

#[test]
fn test_def_a_map() {
    const DEF_A_MAP: ByteIndexMap<4> = ByteIndexMap::new_ignoring_case(*b"ACGT", b'A').add_synonym_ignoring_case(b'U', b'T');
    for i in 0..=255 {
        match i {
            b'C' | b'c' => assert!(DEF_A_MAP.to_index(i) == 1),
            b'G' | b'g' => assert!(DEF_A_MAP.to_index(i) == 2),
            b'T' | b't' | b'U' | b'u' => assert!(DEF_A_MAP.to_index(i) == 3),
            _ => assert!(DEF_A_MAP.to_index(i) == 0),
        }
    }
}

#[test]
fn test_case_sensitive() {
    const DNA_MAP: ByteIndexMap<10> = ByteIndexMap::new(*b"ACGTNacgtn", b'N')
        .add_synonym(b'U', b'T')
        .add_synonym(b'u', b't');

    for i in 0..=255 {
        match i {
            b'A' => assert!(DNA_MAP.to_index(i) == 0),
            b'C' => assert!(DNA_MAP.to_index(i) == 1),
            b'G' => assert!(DNA_MAP.to_index(i) == 2),
            b'T' | b'U' => assert!(DNA_MAP.to_index(i) == 3),
            b'a' => assert!(DNA_MAP.to_index(i) == 5),
            b'c' => assert!(DNA_MAP.to_index(i) == 6),
            b'g' => assert!(DNA_MAP.to_index(i) == 7),
            b't' | b'u' => assert!(DNA_MAP.to_index(i) == 8),
            b'n' => assert!(DNA_MAP.to_index(i) == 9),
            _ => assert!(DNA_MAP.to_index(i) == 4),
        }
    }
}

#[test]
#[should_panic = "assertion failed: array_types::is_unique(&byte_keys)"]
fn test_duplicate() {
    let _ = ByteIndexMap::new(*b"ACGTNA", b'N');
}

#[test]
#[should_panic = "assertion failed: array_types::is_unique(&byte_keys)"]
fn test_duplicate_nocase() {
    let _ = ByteIndexMap::new_ignoring_case(*b"ACGTNA", b'N');
}

#[test]
#[should_panic = "The catch_all must be present in the byte_keys."]
fn test_missing_catch_all() {
    let _ = ByteIndexMap::new(*b"ACGT", b'N');
}

#[test]
#[should_panic = "The catch_all must be present in the byte_keys."]
fn test_missing_catch_all_nocase() {
    let _ = ByteIndexMap::new_ignoring_case(*b"ACGT", b'N');
}

const VALIDATOR_PAIRS: [(&[u8], [bool; 256]); 5] = [
    (DNA_IUPAC_WITH_GAPS, IS_DNA_IUPAC_WITH_GAPS),
    (DNA_IUPAC_NO_GAPS, IS_DNA_IUPAC_NO_GAPS),
    (DNA_IUPAC_NO_GAPS_UC, IS_DNA_IUPAC_NO_GAPS_UC),
    (DNA_ACGTN_NO_GAPS, IS_DNA_ACGTN_NO_GAPS),
    (DNA_ACGTN_NO_GAPS_UC, IS_DNA_ACGTN_NO_GAPS_UC),
];

#[test]
fn test_alphabet_validators() {
    for (alphabet, validator) in VALIDATOR_PAIRS {
        for c in u8::MIN..=u8::MAX {
            assert_eq!(alphabet.contains(&c), validator[c as usize]);
        }
    }
}

const TO_DNA_PAIRS: [(&[u8], [u8; 256]); 2] = [
    (DNA_IUPAC_NO_GAPS, TO_DNA_IUPAC_NO_GAPS_UC),
    (DNA_IUPAC_WITH_GAPS, TO_DNA_IUPAC_WITH_GAPS_UC),
];

#[test]
fn test_to_dna() {
    for (alphabet, converter) in TO_DNA_PAIRS {
        for c in u8::MIN..=u8::MAX {
            let in_alphabet = alphabet.contains(&c);
            if in_alphabet {
                if matches!(c, b'u' | b'U') {
                    assert_eq!(converter[c as usize], b'T');
                } else if c.is_ascii_lowercase() {
                    assert_eq!(converter[c as usize], c.to_ascii_uppercase());
                } else {
                    assert_eq!(converter[c as usize], c);
                }
            } else {
                assert_eq!(converter[c as usize], 0);
            }
        }
    }
}

#[test]
fn test_iupac_to_dna_acgtn_uc() {
    for c in u8::MIN..=u8::MAX {
        if matches!(c, b'u' | b'U') {
            assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC[c as usize], b'T');
        } else if DNA_ACGTN_NO_GAPS.contains(&c) {
            assert_eq!(
                IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC[c as usize],
                c.to_ascii_uppercase()
            );
        } else if DNA_IUPAC_NO_GAPS.contains(&c) {
            assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC[c as usize], b'N');
        } else {
            assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC[c as usize], c);
        }
    }
}

#[test]
fn test_iupac_to_dna_acgtn() {
    for c in u8::MIN..=u8::MAX {
        if c == b'u' {
            assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS[c as usize], b't');
        } else if c == b'U' {
            assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS[c as usize], b'T');
        } else if DNA_ACGTN_NO_GAPS.contains(&c) {
            assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS[c as usize], c);
        } else if DNA_IUPAC_NO_GAPS.contains(&c) {
            if c.is_ascii_lowercase() {
                assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS[c as usize], b'n');
            } else {
                assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS[c as usize], b'N');
            }
        } else {
            assert_eq!(IUPAC_TO_DNA_ACGTN_WITH_GAPS[c as usize], c);
        }
    }
}

#[test]
fn test_any_to_dna_acgtn_uc() {
    for c in u8::MIN..=u8::MAX {
        if matches!(c, b'u' | b'U') {
            assert_eq!(ANY_TO_DNA_ACGTN_NO_GAPS_UC[c as usize], b'T');
        } else if DNA_ACGTN_NO_GAPS.contains(&c) {
            assert_eq!(
                ANY_TO_DNA_ACGTN_NO_GAPS_UC[c as usize],
                c.to_ascii_uppercase()
            );
        } else {
            assert_eq!(ANY_TO_DNA_ACGTN_NO_GAPS_UC[c as usize], b'N');
        }
    }
}

#[test]
fn test_any_to_dna_iupac_with_gaps() {
    for c in u8::MIN..=u8::MAX {
        if c == b'u' {
            assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS[c as usize], b't');
        } else if c == b'U' {
            assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS[c as usize], b'T');
        } else if DNA_IUPAC_WITH_GAPS.contains(&c) {
            assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS[c as usize], c);
        } else {
            assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS[c as usize], b'N');
        }
    }
}

#[test]
fn test_any_to_dna_iupac_with_gaps_uc() {
    for c in u8::MIN..=u8::MAX {
        if matches!(c, b'u' | b'U') {
            assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS_UC[c as usize], b'T');
        } else if DNA_IUPAC_WITH_GAPS.contains(&c) {
            assert_eq!(
                ANY_TO_DNA_IUPAC_WITH_GAPS_UC[c as usize],
                c.to_ascii_uppercase()
            );
        } else {
            assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS_UC[c as usize], b'N');
        }
    }
}

#[test]
fn test_dna_profile_byte_index_counts() {
    let seq = Nucleotides::from(b"TGCTCCACGNGCAGGCGCCTGN");
    let mut counts = ByteIndexCounts::<usize, 5>::new(&DNA_PROFILE_MAP);
    counts.tally_from_seq(seq);

    assert_eq!(counts.into_inner(), [2, 8, 7, 3, 2]);
    assert_eq!(counts.get_count(b'A'), 2);
    assert_eq!(counts.get_count(b'C'), 8);
    assert_eq!(counts.get_count(b'G'), 7);
    assert_eq!(counts.get_count(b'T'), 3);
    assert_eq!(counts.get_count(b'N'), 2);

    counts += b'A';
    counts += b'C';

    assert_eq!(counts.into_inner(), [3, 9, 7, 3, 2]);
    assert_eq!(counts.get_count(b'A'), 3);
    assert_eq!(counts.get_count(b'C'), 9);
    assert_eq!(counts.get_count(b'G'), 7);
    assert_eq!(counts.get_count(b'T'), 3);
    assert_eq!(counts.get_count(b'N'), 2);
}

#[test]
fn test_unambig_dna_profile_byte_index_counts() {
    let seq = Nucleotides::from(b"TGCTCCACGNGCAGGCGCCTGN");
    let mut counts = ByteIndexCounts::<usize, 4>::new(&DNA_UNAMBIG_PROFILE_MAP);
    counts.tally_from_seq(seq);

    assert_eq!(counts.into_inner(), [4, 8, 7, 3]);
    assert_eq!(counts.get_count(b'A'), 4);
    assert_eq!(counts.get_count(b'C'), 8);
    assert_eq!(counts.get_count(b'G'), 7);
    assert_eq!(counts.get_count(b'T'), 3);

    counts += b'A';
    counts += b'C';

    assert_eq!(counts.into_inner(), [5, 9, 7, 3]);
    assert_eq!(counts.get_count(b'A'), 5);
    assert_eq!(counts.get_count(b'C'), 9);
    assert_eq!(counts.get_count(b'G'), 7);
    assert_eq!(counts.get_count(b'T'), 3);
}
