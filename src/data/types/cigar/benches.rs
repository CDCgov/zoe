use super::*;
use test::Bencher;
extern crate test;

#[bench]
fn expand_cigar(b: &mut Bencher) {
    let cigar = Cigar::from_slice_unchecked("3S10M2I2D3M4H4P");

    b.iter(|| cigar.expand_cigar());
}

#[bench]
fn condense_cigar(b: &mut Bencher) {
    let cigar = Cigar::from_slice_unchecked("4S10M2I2D3M4H4P").expand_cigar();
    b.iter(|| cigar.clone().condense_to_cigar());
}

#[bench]
fn match_length(b: &mut Bencher) {
    let cigar = Cigar::from_slice_unchecked("4S10M2I2D3M4H4P");
    b.iter(|| cigar.match_length());
}
