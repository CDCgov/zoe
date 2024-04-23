use super::*;
use crate::data::alphas::NUCLEIC_IUPAC_UNALIGNED;
use std::sync::LazyLock;

const N: usize = 1200;
const SEED: u64 = 42;

static SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(NUCLEIC_IUPAC_UNALIGNED, N, SEED));

#[test]
fn test_translate() {
    let s = Nucleotides(b"ATGTCAGATcccagagaaTGAgg".to_vec());
    assert_eq!(s.translate().as_bytes(), b"MSDPRE*~");
}

#[test]
fn test_translate_collect() {
    use super::super::amino_acids::AminoAcids;
    let s = Nucleotides(b"ATGTCAGATcccagagaaTGAgg".to_vec());
    assert_eq!(s.into_aa_iter().collect::<AminoAcids>().0, b"MSDPRE*~".to_vec());
}

#[test]
fn validate_nucleotides() {
    let mut s: Nucleotides = b"U gotta get my gat back--ok?!".into();
    s.retain_dna_uc();

    assert_eq!(s.to_string(), "TGTTAGTMYGATBACK--K");
}

#[test]
fn simd_reverse_complement() {
    assert_eq!(
        String::from_utf8_lossy(&reverse_complement(&SEQ)),
        String::from_utf8_lossy(&reverse_complement_simd::<32>(&SEQ))
    );
}
