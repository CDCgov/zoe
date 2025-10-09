//! A module providing implementations of
//! [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
//! for many of *Zoe*'s types, as well as wrapper types to provide
//! [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
//! implementations with stronger assumptions.
//!
//! <div class="warning note">
//!
//! **Note**
//!
//! You must enable the *fuzzing* feature in your `Cargo.toml` to use these
//! functions.
//!
//! </div>

use super::{fasta::FastaSeq, fastq::FastQ};
use crate::prelude::*;
use arbitrary::{Arbitrary, Result, Unstructured};
use std::fmt::Display;

mod alignment;
mod kmer;
mod math;
mod sam;
mod string;
mod types;
mod vec;

pub use alignment::*;
pub use math::*;
pub use string::*;
pub use types::*;
pub use vec::*;

macro_rules! impl_deref {
    ($wrapper:ty, $inner:ty $(, $($generics:tt)*)?) => {
        impl$($($generics)*)? ::std::ops::Deref for $wrapper {
            type Target = $inner;

            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl$($($generics)*)? ::std::ops::DerefMut for $wrapper {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
            }
        }
    };
}

macro_rules! impl_from {
    ($wrapper:tt, $($ty:ty),*) => {
        $(
            impl ::std::convert::From<$wrapper<$ty>> for $ty {
                #[inline]
                fn from(val: $wrapper<$ty>) -> $ty {
                    val.0
                }
            }
        )*
    };
}

/// A wrapper around u8 such that the byte is graphic ASCII in the range
/// `!`..=`~`.
pub struct GraphicAsciiByte(u8);

impl_deref! {GraphicAsciiByte, u8}

impl<'a> Arbitrary<'a> for GraphicAsciiByte {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(GraphicAsciiByte(u.int_in_range(b'!'..=b'~')?))
    }

    fn size_hint(depth: usize) -> (usize, Option<usize>) {
        let _ = depth;
        (1, Some(1))
    }
}

/// A wrapper around u8 such that the byte is ASCII in the range
/// `0`..=`127`.
pub struct AsciiByte(u8);

impl_deref! {AsciiByte, u8}

impl<'a> Arbitrary<'a> for AsciiByte {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(AsciiByte(u.int_in_range(0..=127)?))
    }

    fn size_hint(depth: usize) -> (usize, Option<usize>) {
        let _ = depth;
        (1, Some(1))
    }
}

impl<'a> Arbitrary<'a> for QualityScores {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(QualityScores(VecAsciiGraphic::arbitrary(u)?.0))
    }
}

impl<'a> Arbitrary<'a> for FastQ {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let header = String::arbitrary(u)?;
        let sequence = Nucleotides::arbitrary(u)?;
        let quality = QualityScores::arbitrary(u)?;

        Ok(FastQ {
            header,
            sequence,
            quality,
        })
    }
}

/// A wrapper around [`FastQ`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates valid FASTQ records. This means:
///
/// - The header can only contain graphic ASCII or spaces
/// - The sequence must be graphic ASCII in the range `!`..=`~`
/// - The sequence and quality scores must be the same length
#[derive(Debug)]
pub struct FastQValid(pub FastQ);

impl_deref! {FastQValid, FastQ}

impl<'a> Arbitrary<'a> for FastQValid {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let header = StringAsciiGraphicOrSpace::arbitrary(u)?.0;
        let mut sequence = NucleotidesAsciiGraphic::arbitrary(u)?.0;
        let mut quality = QualityScores::arbitrary(u)?;
        let min_len = sequence.len().min(quality.len());
        sequence.shorten_to(min_len);
        quality.shorten_to(min_len);

        Ok(FastQValid(FastQ {
            header,
            sequence,
            quality,
        }))
    }
}

impl Display for FastQValid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl<'a> Arbitrary<'a> for FastaSeq {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let name = String::arbitrary(u)?.replace('>', " ");
        let mut sequence = Vec::<u8>::arbitrary(u)?;
        sequence.replace_all_bytes(b'>', b' ');

        Ok(FastaSeq { name, sequence })
    }
}

pub(crate) use impl_deref;
pub(crate) use impl_from;
