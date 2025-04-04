use test::Bencher;
extern crate test;
use super::*;
use crate::data::{
    alphas::{DNA_ACGTN_NO_GAPS_UC, ENGLISH},
    mappings::TO_DNA_IUPAC_NO_GAPS_UC,
};
use std::sync::LazyLock;

const N: usize = 1200;
const SEED: u64 = 42;

static SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(ENGLISH, N, SEED));
static READ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(DNA_ACGTN_NO_GAPS_UC, 150, SEED));

#[bench]
fn translate_sequence_long(b: &mut Bencher) {
    b.iter(|| translate_sequence(&SEQ));
}

#[bench]
fn validate_retain_iupac_uc(b: &mut Bencher) {
    b.iter(|| {
        SEQ.clone().retain_mut(|b| {
            *b = TO_DNA_IUPAC_NO_GAPS_UC[*b as usize];
            *b > 0
        });
    });
}

#[bench]
fn validate_filtermap_iupac_uc(b: &mut Bencher) {
    b.iter(|| {
        let _: Vec<u8> = SEQ
            .clone()
            .iter_mut()
            .filter_map(|b| {
                *b = TO_DNA_IUPAC_NO_GAPS_UC[*b as usize];
                if *b > 0 { Some(*b) } else { None }
            })
            .collect();
    });
}

#[bench]
fn revcomp_scalar(b: &mut Bencher) {
    b.iter(|| reverse_complement(&SEQ));
}

#[bench]
fn revcomp_simd32(b: &mut Bencher) {
    b.iter(|| reverse_complement_simd::<32>(&SEQ));
}

#[bench]
fn is_acgtn_uc_read_scalar(b: &mut Bencher) {
    let s = READ.as_slice();
    b.iter(|| s.is_valid_dna(IsValidDNA::AcgtnNoGapsUc));
}

#[bench]
fn is_acgtn_uc_read_simd(b: &mut Bencher) {
    let s = READ.as_slice();
    b.iter(|| s.is_acgtn_uc());
}
