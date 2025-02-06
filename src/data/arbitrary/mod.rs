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
//! You must enable the *arbitrary* feature in your `Cargo.toml` to use these functions.
//!
//! </div>

use super::{fasta::FastaSeq, fastq::FastQ};
use crate::prelude::*;
use arbitrary::{Arbitrary, Result, Unstructured};
use std::{
    fmt::Display,
    ops::{Deref, DerefMut},
};

mod nucleotides;
mod string;
mod vec;

pub use nucleotides::*;
pub use string::*;
pub use vec::*;

macro_rules! impl_deref {
    ($wrapper:ty, $inner:ty) => {
        impl Deref for $wrapper {
            type Target = $inner;

            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl DerefMut for $wrapper {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
            }
        }
    };
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
/// * The header must begin with `@`
/// * The header can only contain graphic ASCII or spaces
/// * The sequence must be graphic ASCII in the range `!`..=`~`
/// * The sequence and quality scores must be the same length
#[derive(Debug)]
pub struct FastQValid(pub FastQ);

impl_deref! {FastQValid, FastQ}

impl<'a> Arbitrary<'a> for FastQValid {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let header = "@".to_string() + &StringAsciiGraphicOrSpace::arbitrary(u)?;
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
