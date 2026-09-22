//! An arbitrary implementations and specification struct for [`QualityScores`].

use crate::{
    data::arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecs, Case, VecSpecs},
    prelude::QualityScores,
};
use arbitrary::{Arbitrary, Result, Unstructured};

/// Specifications for generating arbitrary [`QualityScores`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct QualityScoresSpecs {
    /// The minimum number of quality scores to generate.
    ///
    /// This must be less than or equal to `max_len`.
    pub min_len: usize,

    /// The exact number of quality scores to generate.
    ///
    /// If set, this ignores the `min_len` and `max_len` fields.
    pub len: Option<usize>,

    /// The maximum number of quality scores to generate.
    ///
    /// This must be greater than or equal to `min_len`.
    pub max_len: usize,
}

impl Default for QualityScoresSpecs {
    #[inline]
    fn default() -> Self {
        Self {
            min_len: 0,
            len:     None,
            max_len: usize::MAX,
        }
    }
}

impl QualityScoresSpecs {
    /// Returns a copy of the specifications with a fixed length.
    #[inline]
    #[must_use]
    pub const fn with_len(mut self, len: usize) -> Self {
        self.len = Some(len);
        self
    }
}

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

        // Validity: The bytes will only contain graphic ASCII characters per
        // above.
        Ok(QualityScores::try_from(specs.make_arbitrary(u)?).unwrap())
    }
}

impl<'a> ArbitrarySpecs<'a> for QualityScoresSpecs {
    type Output = QualityScores;

    /// Generates arbitrary [`QualityScores`] conforming to the given
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
        let specs = VecSpecs {
            element_specs: ByteSpecs {
                set:  ByteSet::AsciiGraphic,
                case: Case::Any,
            },
            min_len:       self.min_len,
            len:           self.len,
            max_len:       self.max_len,
        };

        // Validity: The bytes will only contain graphic ASCII characters per
        // above.
        Ok(QualityScores::try_from(specs.make_arbitrary(u)?).unwrap())
    }
}
