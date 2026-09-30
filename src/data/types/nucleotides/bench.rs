use test::Bencher;
extern crate test;
use super::*;
use crate::data::{alphas::DNA_ACGTN_UC, mappings::RETAIN_DNA_IUPAC_NO_GAPS_UC};
use std::sync::LazyLock;

/// Upper and lowercase English alphabet.
const ENGLISH: &[u8; 52] = b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

const LEN: usize = 1200;
const SEED: u64 = 42;

static SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(ENGLISH, LEN, SEED));
static READ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(DNA_ACGTN_UC, 150, SEED));

#[bench]
fn translate_sequence_long(b: &mut Bencher) {
    b.iter(|| translate_sequence(&SEQ));
}

#[bench]
fn validate_retain_iupac_uc(b: &mut Bencher) {
    b.iter(|| {
        SEQ.clone().retain_mut(|b| {
            *b = RETAIN_DNA_IUPAC_NO_GAPS_UC[*b];
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
                *b = RETAIN_DNA_IUPAC_NO_GAPS_UC[*b];
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

mod read_recode {
    use crate::DEFAULT_SIMD_LANES;

    use super::*;
    use crate::data::validation::recode::Recode;
    use crate::simd::SimdByteFunctions;
    use std::simd::prelude::*;

    #[bench]
    fn baseline(b: &mut Bencher) {
        let v: Nucleotides = SEQ.to_vec().into();
        b.iter(|| v.clone());
    }

    #[bench]
    fn scalar(b: &mut Bencher) {
        let v: Nucleotides = SEQ.to_vec().into();
        b.iter(|| v.clone().recode_dna(RecodeDNAStrat::AnyToAcgtnNoGapsUpper));
    }

    #[bench]
    fn as_chunks_simd(b: &mut Bencher) {
        let v: Nucleotides = SEQ.to_vec().into();
        b.iter(|| v.clone().recode_dna_reads());
    }

    #[bench]
    fn as_simd(b: &mut Bencher) {
        fn make_acgtn_uc_simd(s: &mut [u8]) {
            const A: Simd<u8, { DEFAULT_SIMD_LANES }> = Simd::splat(b'A');
            const G: Simd<u8, { DEFAULT_SIMD_LANES }> = Simd::splat(b'G');
            const C: Simd<u8, { DEFAULT_SIMD_LANES }> = Simd::splat(b'C');
            const T: Simd<u8, { DEFAULT_SIMD_LANES }> = Simd::splat(b'T');
            const N: Simd<u8, { DEFAULT_SIMD_LANES }> = Simd::splat(b'N');

            let (mut left, mid, mut right) = s.as_simd_mut::<{ DEFAULT_SIMD_LANES }>();

            left.recode(RecodeDNAStrat::AnyToAcgtnNoGapsUpper.mapping());

            for v in mid {
                v.make_ascii_uppercase();
                v.if_value_then_replace(b'U', b'T');
                let valid = v.simd_eq(A) | v.simd_eq(G) | v.simd_eq(C) | v.simd_eq(T);
                *v = valid.select(*v, N);
            }
            right.recode(RecodeDNAStrat::AnyToAcgtnNoGapsUpper.mapping());
        }

        let v = SEQ.to_vec();
        b.iter(|| make_acgtn_uc_simd(&mut v.clone()));
    }
}
