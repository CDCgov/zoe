//! An arbitrary implementations and specification struct for [`QualityScores`].

use crate::{
    data::arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecs, Case, VecSpecs},
    prelude::QualityScores,
};
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a> Arbitrary<'a> for QualityScores {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let specs = VecSpecs {
            element_specs: ByteSpecs {
                set:  ByteSet::AsciiGraphic,
                case: Case::Any,
            },
            min_len:       0,
            len:           None,
            max_len:       usize::MAX,
        };

        // Safety: The bytes will only contain graphic ASCII characters per
        // above.
        specs.make_arbitrary(u).map(QualityScores)
    }
}
