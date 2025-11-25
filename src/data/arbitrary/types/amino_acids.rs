//! An arbitrary implementation and specification struct for [`AminoAcids`].

use crate::{
    data::arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecs, Case, VecSpecs},
    prelude::AminoAcids,
};
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a> Arbitrary<'a> for AminoAcids {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(AminoAcids(Vec::<u8>::arbitrary(u)?))
    }
}

/// Specifications for generating an arbitrary [`AminoAcids`] sequence.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct AminoAcidsSpecs {
    /// The character set to which the `u8` bytes must belong.
    ///
    /// Consider using [`ByteSet::AminoAcidsCanonical`] if needed.
    pub set: ByteSet,

    /// The case of the ASCII characters.
    pub case: Case,
}

impl From<AminoAcidsSpecs> for ByteSpecs {
    #[inline]
    fn from(value: AminoAcidsSpecs) -> Self {
        Self {
            set:  value.set,
            case: value.case,
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for AminoAcidsSpecs {
    type Output = AminoAcids;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let vec_specs = VecSpecs {
            element_specs: ByteSpecs::from(*self),
            ..Default::default()
        };

        vec_specs.make_arbitrary(u).map(AminoAcids)
    }
}
