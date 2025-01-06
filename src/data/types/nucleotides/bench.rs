use test::Bencher;
extern crate test;
use super::*;
use crate::data::alphas::ENGLISH;
use std::sync::LazyLock;

const N: usize = 1200;
const SEED: u64 = 42;

static SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(ENGLISH, N, SEED));

#[bench]
fn translate_sequence_long(b: &mut Bencher) {
    b.iter(|| translate_sequence(&SEQ));
}

#[bench]
fn validate_retain_unaligned_base_uc(b: &mut Bencher) {
    b.iter(|| {
        SEQ.clone().retain_mut(|b| {
            *b = TO_UNALIGNED_DNA_UC[*b as usize];
            *b > 0
        });
    });
}

#[bench]
fn validate_filtermap_unaligned_base_uc(b: &mut Bencher) {
    b.iter(|| {
        let _: Vec<u8> = SEQ
            .clone()
            .iter_mut()
            .filter_map(|b| {
                *b = TO_UNALIGNED_DNA_UC[*b as usize];
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
