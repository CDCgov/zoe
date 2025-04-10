use super::*;
use crate::data::{alphas::DNA_IUPAC_NO_GAPS, types::nucleotides::make_reverse_complement};
use std::sync::LazyLock;

const N: usize = 1200;
const SEED: u64 = 42;

static SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(DNA_IUPAC_NO_GAPS, N, SEED));

#[test]
fn test_translate() {
    let s = Nucleotides(b"ATGTCAGATcccagagaaTGAgg".to_vec());
    assert_eq!(s.translate().as_bytes(), b"MSDPRE*~");
}

#[test]
fn test_translate_collect() {
    use super::super::amino_acids::AminoAcids;
    let s = Nucleotides(b"ATGTCAGATcccagagaaTGAgg".to_vec());
    assert_eq!(s.to_aa_iter().collect::<AminoAcids>().0, b"MSDPRE*~".to_vec());
}

#[test]
fn test_translate_collect_with_partial_codon() {
    use super::super::amino_acids::AminoAcids;
    let s = Nucleotides(b"ATGTCAGATcccagagaaTGAgg".to_vec());
    assert_eq!(s.to_aa_iter_with(b'X').collect::<AminoAcids>().0, b"MSDPRE*X".to_vec());
}

#[test]
fn retain_and_recode_iupac_with_gaps_uc() {
    let mut s: Nucleotides = b"U gotta get my gat back--ok?!".into();
    s.retain_and_recode_dna(RefineDNAStrat::IupacWithGapsUc);
    assert_eq!(s.to_string(), "TGTTAGTMYGATBACK--K");
}

#[test]
fn sanitize_dna() {
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::fn_params_excessive_bools)]
    fn check_is_valid_dna(
        seq: &Vec<u8>, iupac_no_gaps: bool, iupac_no_gaps_uc: bool, iupac_with_gaps: bool, iupac_with_gaps_uc: bool,
        acgtn_no_gaps: bool, acgtn_no_gaps_uc: bool, acgtn_std_gaps_uc: bool,
    ) {
        assert_eq!(seq.is_valid_dna(IsValidDNA::IupacNoGaps), iupac_no_gaps);
        assert_eq!(seq.is_valid_dna(IsValidDNA::IupacNoGapsUc), iupac_no_gaps_uc);
        assert_eq!(seq.is_valid_dna(IsValidDNA::IupacWithGaps), iupac_with_gaps);
        assert_eq!(seq.is_valid_dna(IsValidDNA::IupacWithGapsUc), iupac_with_gaps_uc);
        assert_eq!(seq.is_valid_dna(IsValidDNA::AcgtnNoGaps), acgtn_no_gaps);
        assert_eq!(seq.is_valid_dna(IsValidDNA::AcgtnNoGapsUc), acgtn_no_gaps_uc);
        assert_eq!(seq.is_valid_dna(IsValidDNA::AcgtnStdGapsUc), acgtn_std_gaps_uc);
    }

    for base in u8::MIN..=u8::MAX {
        let seq = vec![base];

        assert_eq!(seq.is_acgtn_uc(), seq.is_valid_dna(IsValidDNA::AcgtnNoGapsUc));

        match base {
            b'A' | b'C' | b'G' | b'T' | b'N' => {
                check_is_valid_dna(&seq, true, true, true, true, true, true, true);
            }
            b'a' | b'c' | b'g' | b't' | b'n' => {
                check_is_valid_dna(&seq, true, false, true, false, true, false, false);
            }
            b'B' | b'D' | b'H' | b'K' | b'M' | b'R' | b'S' | b'U' | b'V' | b'W' | b'Y' => {
                check_is_valid_dna(&seq, true, true, true, true, false, false, false);
            }
            b'b' | b'd' | b'h' | b'k' | b'm' | b'r' | b's' | b'u' | b'v' | b'w' | b'y' => {
                check_is_valid_dna(&seq, true, false, true, false, false, false, false);
            }
            b'-' => {
                check_is_valid_dna(&seq, false, false, true, true, false, false, true);
            }
            b'.' => {
                check_is_valid_dna(&seq, false, false, true, true, false, false, false);
            }
            _ => {
                check_is_valid_dna(&seq, false, false, false, false, false, false, false);
            }
        }
    }
}

#[test]
fn check_retain_dna() {
    for base in u8::MIN..=u8::MAX {
        for strategy in [
            IsValidDNA::IupacNoGaps,
            IsValidDNA::IupacNoGapsUc,
            IsValidDNA::IupacWithGaps,
            IsValidDNA::IupacWithGapsUc,
            IsValidDNA::AcgtnNoGaps,
            IsValidDNA::AcgtnNoGapsUc,
            IsValidDNA::AcgtnStdGapsUc,
        ] {
            let mut seq = Nucleotides::from(&[base]);
            if seq.is_valid_dna(strategy) {
                seq.retain_dna(strategy);
                assert_eq!(seq, Nucleotides::from(&[base]));
            } else {
                seq.retain_dna(strategy);
                assert!(seq.is_empty());
            }
        }
    }
}

#[test]
fn check_refine_dna() {
    fn get_refined(seq: u8, strategy: RefineDNAStrat) -> Vec<u8> {
        let mut seq = Nucleotides::from(&[seq]);
        seq.retain_and_recode_dna(strategy);
        seq.into_vec()
    }

    for base in u8::MIN..=u8::MAX {
        for (refine_strat, valid_strat) in [
            (RefineDNAStrat::IupacNoGapsUc, IsValidDNA::IupacNoGaps),
            (RefineDNAStrat::IupacWithGapsUc, IsValidDNA::IupacWithGaps),
            (RefineDNAStrat::AcgtnNoGapsUc, IsValidDNA::AcgtnNoGaps),
        ] {
            let seq = get_refined(base, refine_strat);
            if base.eq_ignore_ascii_case(&b'U') {
                assert_eq!(seq, b"T", "{} vs T in {valid_strat:?}", seq[0] as char);
            } else if valid_strat.is_valid(base) {
                let b = base.to_ascii_uppercase();
                assert_eq!(seq, &[b], "{} vs {} in {valid_strat:?}", seq[0] as char, b as char);
            } else {
                assert!(seq.is_empty(), "ø vs {} in {valid_strat:?}", base as char,);
            }
        }

        let (valid_strat, refine_strat) = (IsValidDNA::IupacNoGaps, RefineDNAStrat::IupacCorrectGapsUc);
        let seq = get_refined(base, refine_strat);
        if base.eq_ignore_ascii_case(&b'U') {
            assert_eq!(seq, b"T");
        } else if valid_strat.is_valid(base) {
            assert_eq!(seq, &[base.to_ascii_uppercase()]);
        } else if base == b':' || base == b'~' || base == b'-' {
            assert_eq!(seq, b"-");
        } else if base == b'.' {
            assert_eq!(seq, b".");
        } else {
            assert!(seq.is_empty());
        }

        let (valid_strat, refine_strat) = (IsValidDNA::AcgtnNoGaps, RefineDNAStrat::AcgtnWithGapsUc);
        let seq = get_refined(base, refine_strat);
        if base.eq_ignore_ascii_case(&b'U') {
            assert_eq!(seq, b"T", "{} vs T in {valid_strat:?}", seq[0] as char);
        } else if valid_strat.is_valid(base) {
            assert_eq!(seq, &[base.to_ascii_uppercase()]);
        } else if base == b'-' || base == b'.' {
            assert_eq!(seq, &[base]);
        } else {
            assert!(seq.is_empty());
        }

        let (valid_strat, refine_strat) = (IsValidDNA::AcgtnNoGaps, RefineDNAStrat::AcgtnStdGapsUc);
        let seq = get_refined(base, refine_strat);
        if base.eq_ignore_ascii_case(&b'U') {
            assert_eq!(seq, b"T", "{} vs T in {valid_strat:?}", seq[0] as char);
        } else if valid_strat.is_valid(base) {
            assert_eq!(seq, &[base.to_ascii_uppercase()]);
        } else if base == b'-' || base == b'.' || base == b':' || base == b'~' {
            assert_eq!(seq, b"-");
        } else {
            assert!(seq.is_empty());
        }
    }
}

#[test]
fn check_recode_dna() {
    fn get_recoded(seq: u8, strategy: RecodeDNAStrat) -> Vec<u8> {
        let mut seq = Nucleotides::from(&[seq]);
        seq.recode_dna(strategy);
        seq.into_vec()
    }

    for char in u8::MIN..=u8::MAX {
        let seq = vec![char];
        let recode = get_recoded(char, RecodeDNAStrat::IupacToAcgtnWithGaps);
        if seq.is_valid_dna(IsValidDNA::AcgtnNoGaps) || seq == b"-" || seq == b"." {
            assert_eq!(recode, seq);
        } else if seq == b"U" {
            assert_eq!(recode, b"T");
        } else if seq == b"u" {
            assert_eq!(recode, b"t");
        } else if seq.is_valid_dna(IsValidDNA::IupacWithGaps) {
            if seq[0].is_ascii_lowercase() {
                assert_eq!(&recode, b"n");
            } else {
                assert_eq!(&recode, b"N", "{}", seq[0]);
            }
        } else {
            assert_eq!(recode, seq);
        }

        let recode = get_recoded(char, RecodeDNAStrat::IupacToAcgtnWithGapsUpper);
        if seq.is_valid_dna(IsValidDNA::AcgtnNoGaps) || seq == b"-" || seq == b"." {
            assert_eq!(recode, seq.to_ascii_uppercase());
        } else if seq == b"U" || seq == b"u" {
            assert_eq!(recode, b"T");
        } else if seq.is_valid_dna(IsValidDNA::IupacWithGaps) {
            assert_eq!(&recode, b"N");
        } else {
            assert_eq!(recode, seq);
        }

        let recode = get_recoded(char, RecodeDNAStrat::AnyToAcgtnNoGapsUpper);
        if seq.is_valid_dna(IsValidDNA::AcgtnNoGaps) {
            assert_eq!(recode, seq.to_ascii_uppercase());
        } else if seq == b"U" || seq == b"u" {
            assert_eq!(recode, b"T");
        } else {
            assert_eq!(&recode, b"N");
        }

        let recode = get_recoded(char, RecodeDNAStrat::AnyToAcgtnWithGapsUpper);
        if seq.is_valid_dna(IsValidDNA::AcgtnNoGaps) || seq == b"-" || seq == b"." {
            assert_eq!(recode, seq.to_ascii_uppercase());
        } else if seq == b"U" || seq == b"u" {
            assert_eq!(recode, b"T");
        } else {
            assert_eq!(&recode, b"N");
        }

        let recode = get_recoded(char, RecodeDNAStrat::AnyToIupacWithGaps);
        if seq == b"U" {
            assert_eq!(recode, b"T");
        } else if seq == b"u" {
            assert_eq!(recode, b"t");
        } else if seq.is_valid_dna(IsValidDNA::IupacWithGaps) {
            assert_eq!(recode, seq);
        } else {
            assert_eq!(&recode, b"N");
        }

        let recode = get_recoded(char, RecodeDNAStrat::AnyToIupacWithGapsUpper);
        if seq == b"U" || seq == b"u" {
            assert_eq!(recode, b"T");
        } else if seq.is_valid_dna(IsValidDNA::IupacWithGaps) {
            assert_eq!(recode, seq.to_ascii_uppercase());
        } else {
            assert_eq!(&recode, b"N");
        }

        let recode = get_recoded(char, RecodeDNAStrat::AnyToIupacCorrectGapsUpper);
        if seq == b"U" || seq == b"u" {
            assert_eq!(recode, b"T");
        } else if seq.is_valid_dna(IsValidDNA::IupacWithGaps) {
            assert_eq!(recode, seq.to_ascii_uppercase());
        } else if &seq == b"~" || &seq == b":" {
            assert_eq!(&recode, b"-");
        } else {
            assert_eq!(&recode, b"N");
        }
    }
}

#[test]
fn simd_reverse_complement() {
    assert_eq!(
        String::from_utf8_lossy(&reverse_complement(&SEQ)),
        String::from_utf8_lossy(&reverse_complement_simd::<32>(&SEQ))
    );
}

#[test]
fn test_make_reverse_complement() {
    let mut seq2 = SEQ.clone();
    make_reverse_complement(&mut seq2);

    assert_eq!(
        String::from_utf8_lossy(&reverse_complement(&SEQ)),
        String::from_utf8_lossy(&seq2)
    );
}

#[test]
fn test_simd_make_reverse_complement() {
    let mut seq2 = SEQ.clone();
    make_reverse_complement_simd::<32>(&mut seq2);

    assert_eq!(
        String::from_utf8_lossy(&reverse_complement(&SEQ)),
        String::from_utf8_lossy(&seq2)
    );
}
