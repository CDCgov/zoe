use super::{GraphicAsciiByte, impl_deref};
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
            vec.push(GraphicAsciiByte::arbitrary(u)?.0);
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

/// A wrapper around `Vec<T>` such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates vecs with length between `MIN_LEN` and `MAX_LEN`.
#[derive(Debug)]
pub struct VecBounded<const MIN_LEN: usize, const MAX_LEN: usize, T>(pub Vec<T>);

impl<const MIN_LEN: usize, const MAX_LEN: usize, T> Deref for VecBounded<MIN_LEN, MAX_LEN, T> {
    type Target = Vec<T>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const MIN_LEN: usize, const MAX_LEN: usize, T> DerefMut for VecBounded<MIN_LEN, MAX_LEN, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<'a, const MIN_LEN: usize, const MAX_LEN: usize, T> Arbitrary<'a> for VecBounded<MIN_LEN, MAX_LEN, T>
where
    T: Arbitrary<'a>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let vec_len = u.arbitrary_len::<T>()?.min(MAX_LEN).max(MIN_LEN);

        let mut vec = Vec::with_capacity(vec_len);
        for _ in 0..vec_len {
            vec.push(T::arbitrary(u)?);
        }

        Ok(VecBounded(vec))
    }
}

impl<const MIN_LEN: usize, const MAX_LEN: usize, T> AsRef<[T]> for VecBounded<MIN_LEN, MAX_LEN, T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0.as_ref()
    }
}

/// A struct holding a vec of vecs such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates vecs of the same length, bounded between `MIN_LEN` and
/// `MAX_LEN`.
#[derive(Debug)]
pub struct SameSizeVecs<const MIN_LEN: usize, const MAX_LEN: usize, T> {
    pub vecs:    Vec<VecBounded<MIN_LEN, MAX_LEN, T>>,
    pub vec_len: usize,
}

impl<const MIN_LEN: usize, const MAX_LEN: usize, T> SameSizeVecs<MIN_LEN, MAX_LEN, T> {
    #[allow(clippy::missing_errors_doc)]
    pub fn arbitrary_with_len<'a>(u: &mut Unstructured<'a>, vec_len: usize) -> Result<Self>
    where
        T: Arbitrary<'a>, {
        let num_vecs = u.arbitrary_len::<T>()? / vec_len;

        let mut vecs = Vec::with_capacity(num_vecs);
        for _ in 0..num_vecs {
            let mut vec = Vec::with_capacity(vec_len);
            for _ in 0..vec_len {
                vec.push(T::arbitrary(u)?);
            }
            vecs.push(VecBounded(vec));
        }

        Ok(SameSizeVecs { vecs, vec_len })
    }
}

impl<'a, const MIN_LEN: usize, const MAX_LEN: usize, T> Arbitrary<'a> for SameSizeVecs<MIN_LEN, MAX_LEN, T>
where
    T: Arbitrary<'a>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let vec_len = usize::arbitrary(u)? % (MAX_LEN + 1 - MIN_LEN) + MIN_LEN;
        SameSizeVecs::arbitrary_with_len(u, vec_len)
    }
}
