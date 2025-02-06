use super::{VecAsciiGraphic, VecAsciiGraphicBashSafe, VecAsciiGraphicOrSpace, impl_deref};
use arbitrary::{Arbitrary, Result, Unstructured};
use std::ops::{Deref, DerefMut};

/// A wrapper around [`String`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates graphic ASCII in the range `!`..=`~`.
#[derive(Debug)]
pub struct StringAsciiGraphic(pub String);

impl_deref! {StringAsciiGraphic, String}

impl<'a> Arbitrary<'a> for StringAsciiGraphic {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(StringAsciiGraphic(
            String::from_utf8_lossy(&VecAsciiGraphic::arbitrary(u)?).to_string(),
        ))
    }
}

/// A wrapper around [`String`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates graphic ASCII in the range `!`..=`~` or spaces.
#[derive(Debug)]
pub struct StringAsciiGraphicOrSpace(pub String);

impl_deref! {StringAsciiGraphicOrSpace, String}

impl<'a> Arbitrary<'a> for StringAsciiGraphicOrSpace {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(StringAsciiGraphicOrSpace(
            String::from_utf8_lossy(&VecAsciiGraphicOrSpace::arbitrary(u)?).to_string(),
        ))
    }
}

/// A wrapper around [`String`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates graphic ASCII in the range `!`..=`~` which also is not a
/// special character in bash.
///
/// The allowed characters are
/// `%+,-./0123456789:@ABCDEFGHIJKLMNOPQRSTUVWXYZ^_abcdefghijklmnopqrstuvwxyz`.
#[derive(Debug)]
pub struct StringAsciiGraphicBashSafe(pub String);

impl_deref! {StringAsciiGraphicBashSafe, String}

impl<'a> Arbitrary<'a> for StringAsciiGraphicBashSafe {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(StringAsciiGraphicBashSafe(
            String::from_utf8_lossy(&VecAsciiGraphicBashSafe::arbitrary(u)?).to_string(),
        ))
    }
}
