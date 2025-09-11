use arbitrary::{Arbitrary, Result, Unstructured};

use crate::prelude::AminoAcids;

impl<'a> Arbitrary<'a> for AminoAcids {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(AminoAcids(Vec::<u8>::arbitrary(u)?))
    }
}
