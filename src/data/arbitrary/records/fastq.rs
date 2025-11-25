//! Arbitrary implementations and specification structs for the [`FastQ`] record
//! type.

use crate::{
    data::arbitrary::{ArbitrarySpecs, ByteSet, Case, NucleotidesSpecs, StringSpecs},
    prelude::{FastQ, Len, Nucleotides, QualityScores},
};
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a> Arbitrary<'a> for FastQ {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let header = String::arbitrary(u)?;
        let sequence = Nucleotides::arbitrary(u)?;
        let quality = QualityScores::arbitrary(u)?;

        Ok(FastQ {
            header,
            sequence,
            quality,
        })
    }
}

/// Specifications for generating an arbitrary [`FastQ`] record.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastQSpecs {
    /// The specifications for generating the header.
    pub header_specs: StringSpecs,

    /// The specifications for generating the sequence.
    pub sequence_specs: NucleotidesSpecs,

    /// Whether to enforce that the sequence and the quality scores are the same
    /// length.
    pub same_lengths: bool,
}

impl FastQSpecs {
    /// A constructor for [`FastQSpecs`] that ensures the generated [`FastQ`]
    /// can be parsed by [`FastQReader`].
    ///
    /// This includes ensuring that the header only contains graphic ASCII or
    /// spaces, that the sequence only contains graphic ASCII, and that the
    /// sequence and quality scores are the same length.
    ///
    /// [`FastQReader`]: crate::prelude::FastQReader
    #[inline]
    #[must_use]
    pub fn parsable_fastq() -> Self {
        Self {
            header_specs:   StringSpecs {
                set:  ByteSet::AsciiGraphicOrSpace,
                case: Case::Any,
            },
            sequence_specs: NucleotidesSpecs {
                set:  ByteSet::AsciiGraphic,
                case: Case::Any,
            },
            same_lengths:   true,
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for FastQSpecs {
    type Output = FastQ;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let header = self.header_specs.make_arbitrary(u)?;
        let mut sequence = self.sequence_specs.make_arbitrary(u)?;
        let mut quality = QualityScores::arbitrary(u)?;
        if self.same_lengths {
            let min_len = sequence.len().min(quality.len());
            sequence.shorten_to(min_len);
            quality.shorten_to(min_len);
        }
        Ok(FastQ {
            header,
            sequence,
            quality,
        })
    }
}
