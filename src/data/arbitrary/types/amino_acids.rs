//! An arbitrary implementation and specification struct for [`AminoAcids`].

use crate::{
    data::arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecs, ByteSpecsView, Case, VecSpecs},
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
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct AminoAcidsSpecs {
    /// The character set to which the `u8` bytes must belong.
    ///
    /// Consider using [`ByteSet::AminoAcidsCanonical`] if needed.
    pub set: ByteSet,

    /// The case of the ASCII characters.
    pub case: Case,

    /// The minimum length of the [`AminoAcids`].
    ///
    /// This must be less than or equal to `max_len`.
    pub min_len: usize,

    /// The exact length of the [`AminoAcids`] to generate.
    ///
    /// If set, this ignores the `min_len` and `max_len` fields.
    pub len: Option<usize>,

    /// The maximum length of the [`AminoAcids`].
    ///
    /// This must be greater than or equal to `min_len`.
    pub max_len: usize,
}

impl Default for AminoAcidsSpecs {
    fn default() -> Self {
        Self {
            set:     ByteSet::default(),
            case:    Case::default(),
            min_len: 0,
            len:     None,
            max_len: usize::MAX,
        }
    }
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

    /// Generates an arbitrary [`AminoAcidsSpecs`] conforming to the given
    /// specifications.
    ///
    /// ## Errors
    ///
    /// Any errors from the underlying [`arbitrary`] calls are propagated.
    ///
    /// ## Panics
    ///
    /// `min_len` must be less than or equal to `max_len`.
    ///
    /// [`arbitrary`]: arbitrary::Arbitrary::arbitrary
    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let vec_specs = VecSpecs {
            element_specs: ByteSpecsView {
                set:  &self.set,
                case: self.case,
            },
            min_len:       self.min_len,
            len:           self.len,
            max_len:       self.max_len,
        };

        vec_specs.make_arbitrary(u).map(AminoAcids)
    }
}
