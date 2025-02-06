use super::impl_deref;
use arbitrary::{Arbitrary, Result, Unstructured};
use std::ops::{Deref, DerefMut};

/// A wrapper around `Vec<u8>` such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates graphic ASCII in the range `!`..=`~`.
#[derive(Debug)]
pub struct VecAsciiGraphic(pub Vec<u8>);

impl_deref! {VecAsciiGraphic, Vec<u8>}

impl<'a> Arbitrary<'a> for VecAsciiGraphic {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let len = u.arbitrary_len::<u8>()?;
        let mut vec = Vec::with_capacity(len);
        for _ in 0..len {
            vec.push(u.int_in_range(b'!'..=b'~')?);
        }
        Ok(VecAsciiGraphic(vec))
    }
}

/// A wrapper around `Vec<u8>` such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates graphic ASCII in the range `!`..=`~` or spaces.
#[derive(Debug)]
pub struct VecAsciiGraphicOrSpace(pub Vec<u8>);

impl_deref! {VecAsciiGraphicOrSpace, Vec<u8>}

impl<'a> Arbitrary<'a> for VecAsciiGraphicOrSpace {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let len = u.arbitrary_len::<u8>()?;
        let mut vec = Vec::with_capacity(len);
        for _ in 0..len {
            vec.push(u.int_in_range(b' '..=b'~')?);
        }
        Ok(VecAsciiGraphicOrSpace(vec))
    }
}

/// A wrapper around `Vec<u8>` such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates graphic ASCII in the range `!`..=`~` which also is not a
/// special character in bash.
///
/// The allowed symbols are
/// `%+,-./0123456789:@ABCDEFGHIJKLMNOPQRSTUVWXYZ^_abcdefghijklmnopqrstuvwxyz`.
#[derive(Debug)]
pub struct VecAsciiGraphicBashSafe(pub Vec<u8>);

impl_deref! {VecAsciiGraphicBashSafe, Vec<u8>}

impl<'a> Arbitrary<'a> for VecAsciiGraphicBashSafe {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        const ALPHA: &[u8] = b"%+,-./0123456789:@ABCDEFGHIJKLMNOPQRSTUVWXYZ^_abcdefghijklmnopqrstuvwxyz";
        Ok(VecAsciiGraphicBashSafe(
            u.arbitrary_iter::<usize>()?
                .flatten()
                .map(|b| ALPHA[b % ALPHA.len()])
                .collect(),
        ))
    }
}
