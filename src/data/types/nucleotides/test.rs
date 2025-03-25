use super::*;
use crate::data::{alphas::DNA_IUPAC, types::nucleotides::make_reverse_complement};
use std::sync::LazyLock;

const N: usize = 1200;
const SEED: u64 = 42;

static SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(DNA_IUPAC, N, SEED));

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
fn validate_iupac_with_gaps_uc() {
    let mut s: Nucleotides = b"U gotta get my gat back--ok?!".into();
    s.retain_iupac_with_gaps_uc();

    assert_eq!(s.to_string(), "TGTTAGTMYGATBACK--K");
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
