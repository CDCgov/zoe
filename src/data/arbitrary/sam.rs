use crate::{
    data::{cigar::Cigar, sam::SamData},
    prelude::{Nucleotides, QualityScores},
};
use arbitrary::Arbitrary;

impl<'a> Arbitrary<'a> for SamData {
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
