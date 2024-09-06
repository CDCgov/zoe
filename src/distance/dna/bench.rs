use test::Bencher;
extern crate test;
use super::{
    test::{LONG_READ_1, LONG_READ_2},
    *,
};

#[bench]
fn jukes_cantor_benchmark(b: &mut Bencher) {
    b.iter(|| jukes_cantor_69(LONG_READ_1, LONG_READ_2));
}

#[bench]
fn k80_benchmark(b: &mut Bencher) {
    b.iter(|| kimura_80(LONG_READ_1, LONG_READ_2));
}

#[bench]
fn k81_benchmark(b: &mut Bencher) {
    b.iter(|| kimura_81(LONG_READ_1, LONG_READ_2));
}

#[bench]
fn f81_benchmark(b: &mut Bencher) {
    b.iter(|| felsenstein_81(LONG_READ_1, LONG_READ_2));
}

#[bench]
fn tn93_benchmark(b: &mut Bencher) {
    b.iter(|| tamura_nei_93(LONG_READ_1, LONG_READ_2));
}
