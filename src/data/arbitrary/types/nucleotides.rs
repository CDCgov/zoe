//! An arbitrary implementations and specification struct for [`Nucleotides`].

use crate::{
    data::arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecs, Case, VecSpecs},
    prelude::Nucleotides,
};
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a> Arbitrary<'a> for Nucleotides {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        u.arbitrary().map(Nucleotides)
    }
}

/// Specifications for generating an arbitrary [`Nucleotides`] sequence.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct NucleotidesSpecs {
    /// The character set to which the `u8` bytes must belong.
    ///
    /// Consider using [`ByteSet::NucleotidesCanonical`] or
    /// [`ByteSet::NucleotidesIupac`] if needed.
    pub set: ByteSet,

    /// The case of the ASCII characters.
    pub case: Case,
}

impl From<NucleotidesSpecs> for ByteSpecs {
    #[inline]
    fn from(value: NucleotidesSpecs) -> Self {
        Self {
            set:  value.set,
            case: value.case,
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for NucleotidesSpecs {
    type Output = Nucleotides;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let vec_specs = VecSpecs {
            element_specs: ByteSpecs::from(*self),
            min_len:       0,
            len:           None,
            max_len:       usize::MAX,
        };

        vec_specs.make_arbitrary(u).map(Nucleotides)
    }
}
