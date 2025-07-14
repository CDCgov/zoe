//! Arbitrary implementations and specification structs for the [`SamData`]
//! record type.

use crate::{
    data::{
        arbitrary::{ArbitrarySpecs, CigarSpecs, ClampAlignment, NucleotidesSpecs, StringSpecs},
        cigar::{Cigar, LenInAlignment},
        sam::SamData,
    },
    prelude::{Len, Nucleotides, QualityScores},
};
use arbitrary::{Arbitrary, Unstructured};
use std::num::NonZeroUsize;

impl<'a> Arbitrary<'a> for SamData {
    /// Generates an arbitrary [`SamData`] record from the given unstructured
    /// data.
    ///
    /// This ensures that `rnext` is `*`, `pnext` is 0, `tlen` is 0, and `tags`
    /// is empty.
    #[inline]
    fn arbitrary(u: &mut arbitrary::Unstructured<'a>) -> arbitrary::Result<Self> {
        Ok(SamData::new(
            String::arbitrary(u)?,
            u16::arbitrary(u)?,
            String::arbitrary(u)?,
            usize::arbitrary(u)?,
            u8::arbitrary(u)?,
            Cigar::arbitrary(u)?,
            Nucleotides::arbitrary(u)?,
            QualityScores::arbitrary(u)?,
        ))
    }
}

/// Specifications for generating arbitrary [`SamData`] records.
///
/// All generated records have `rnext` as `*`, `pnext` as 0, `tlen` as 0, and
/// `tags` being empty.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct SamDataSpecs {
    /// The specifications for generating the `qname` field.
    pub qname: StringSpecs,

    /// The specifications for generating the `rname` field.
    pub rname: StringSpecs,

    /// The specifications for generating the `cigar` field.
    pub cigar: CigarSpecs,

    /// The specifications for generating the `seq` field.
    pub seq: NucleotidesSpecs,

    /// Whether to ensure that the `pos` field is nonzero (since it represents a
    /// 1-based position).
    pub nonzero_pos: bool,

    /// Ensures that the `pos` field is generated such that the end position of
    /// the alignment in the reference does not overflow.
    ///
    /// This assumes the match length of the CIGAR string also does not overflow
    /// (and that it does not equal [`usize::MAX`] when `nonzero_pos` is also
    /// specified).
    pub cap_end_pos: bool,

    /// Ensures that the length of the `seq` field agrees with the CIGAR string.
    ///
    /// As a side effect, this will cause the CIGAR string to only contain valid
    /// ciglets (the CIGAR string is iterated over and reconstructed with
    /// [`CigletIterator`], so see its documentation for more details).
    ///
    /// This may also cause the `seq` field to get shrunk. This may
    /// hypothetically invalidate arbitrary assumptions imposed by
    /// [`NucleotidesSpecs`].
    ///
    /// [`CigletIterator`]: crate::data::types::cigar::CigletIterator
    pub correct_seq_len: bool,
}

impl<'a> ArbitrarySpecs<'a> for SamDataSpecs {
    type Output = SamData;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        let mut cigar = self.cigar.make_arbitrary(u)?;
        let mut seq = self.seq.make_arbitrary(u)?;

        if self.correct_seq_len {
            cigar.clamp_query_len(seq.len());
            seq.shorten_to(cigar.query_len_in_alignment());
        }

        let pos = match (self.nonzero_pos, self.cap_end_pos) {
            (false, false) => usize::arbitrary(u)?,
            (false, true) => {
                if let Some(match_len) = cigar.ref_len_in_alignment_checked() {
                    u.int_in_range(0..=(usize::MAX - match_len))?
                } else {
                    usize::arbitrary(u)?
                }
            }
            (true, false) => NonZeroUsize::arbitrary(u)?.into(),
            (true, true) => {
                if let Some(match_len) = cigar.ref_len_in_alignment_checked()
                    && match_len < usize::MAX
                {
                    u.int_in_range(1..=(usize::MAX - match_len))?
                } else {
                    NonZeroUsize::arbitrary(u)?.get()
                }
            }
        };

        Ok(SamData::new(
            self.qname.make_arbitrary(u)?,
            u16::arbitrary(u)?,
            self.rname.make_arbitrary(u)?,
            pos,
            u8::arbitrary(u)?,
            cigar,
            seq,
            QualityScores::arbitrary(u)?,
        ))
    }
}
