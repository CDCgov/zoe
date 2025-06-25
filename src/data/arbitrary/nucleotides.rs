use super::{VecAsciiGraphic, VecAsciiGraphicBashSafe, impl_deref};
use crate::prelude::Nucleotides;
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a> Arbitrary<'a> for Nucleotides {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Nucleotides(Vec::<u8>::arbitrary(u)?))
    }
}

/// A wrapper around [`Nucleotides`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates graphic ASCII in the range `!`..=`~`.
#[derive(Debug, Clone)]
pub struct NucleotidesAsciiGraphic(pub Nucleotides);

impl_deref! {NucleotidesAsciiGraphic, Nucleotides}

impl<'a> Arbitrary<'a> for NucleotidesAsciiGraphic {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(NucleotidesAsciiGraphic(Nucleotides(VecAsciiGraphic::arbitrary(u)?.0)))
    }
}

/// A wrapper around [`Nucleotides`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates graphic ASCII in the range `!`..=`~` which also is not a
/// special character in bash.
///
/// The allowed symbols are
/// `%+,-./0123456789:@ABCDEFGHIJKLMNOPQRSTUVWXYZ^_abcdefghijklmnopqrstuvwxyz`.
#[derive(Debug, Clone)]
pub struct NucleotidesAsciiGraphicBashSafe(pub Nucleotides);

impl_deref! {NucleotidesAsciiGraphicBashSafe, Nucleotides}

impl<'a> Arbitrary<'a> for NucleotidesAsciiGraphicBashSafe {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(NucleotidesAsciiGraphicBashSafe(Nucleotides(
            VecAsciiGraphicBashSafe::arbitrary(u)?.0,
        )))
    }
}

/// A wrapper around [`Nucleotides`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates bases in `ACGTUNacgtun`.
#[derive(Debug, Clone)]
pub struct NucleotidesAcgtun(pub Nucleotides);

impl_deref! {NucleotidesAcgtun, Nucleotides}

impl<'a> Arbitrary<'a> for NucleotidesAcgtun {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        const ALPHA: &[u8] = b"ACGTUNacgtun";
        Ok(NucleotidesAcgtun(
            u.arbitrary_iter::<u8>()?
                .flatten()
                .map(|b| ALPHA[b as usize % ALPHA.len()])
                .collect(),
        ))
    }
}

/// A wrapper around [`Nucleotides`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates bases in `ACGTUacgtu`.
#[derive(Debug, Clone)]
pub struct NucleotidesAcgtu(pub Nucleotides);

impl_deref! {NucleotidesAcgtu, Nucleotides}

impl<'a> Arbitrary<'a> for NucleotidesAcgtu {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        const ALPHA: &[u8] = b"ACGTUacgtu";
        Ok(NucleotidesAcgtu(
            u.arbitrary_iter::<u8>()?
                .flatten()
                .map(|b| ALPHA[b as usize % ALPHA.len()])
                .collect(),
        ))
    }
}

/// A wrapper around [`Nucleotides`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates bases in `ACGT`.
#[derive(Debug, Clone)]
pub struct NucleotidesAcgtUpper(pub Nucleotides);

impl_deref! {NucleotidesAcgtUpper, Nucleotides}

impl<'a> Arbitrary<'a> for NucleotidesAcgtUpper {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        const ALPHA: &[u8] = b"ACGT";
        Ok(NucleotidesAcgtUpper(
            u.arbitrary_iter::<u8>()?
                .flatten()
                .map(|b| ALPHA[b as usize % ALPHA.len()])
                .collect(),
        ))
    }
}
