use crate::{
    data::{arbitrary::impl_deref, views::impl_len_for_wrapper},
    prelude::AminoAcids,
};
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a> Arbitrary<'a> for AminoAcids {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(AminoAcids(Vec::<u8>::arbitrary(u)?))
    }
}

/// A wrapper around [`AminoAcids`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates bases in `ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx`.
#[derive(Debug, Clone)]
pub struct AminoAcidsIupacX(pub AminoAcids);

impl_deref! {AminoAcidsIupacX, AminoAcids}
impl_len_for_wrapper! {AminoAcidsIupacX, 0}

impl<'a> Arbitrary<'a> for AminoAcidsIupacX {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        const ALPHA: &[u8] = b"ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx";
        Ok(AminoAcidsIupacX(
            u.arbitrary_iter::<u8>()?
                .flatten()
                .map(|b| ALPHA[b as usize % ALPHA.len()])
                .collect(),
        ))
    }
}
