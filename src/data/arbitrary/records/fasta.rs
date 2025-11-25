//! An arbitrary implementation for the [`FastaSeq`] record type.

use crate::{data::fasta::FastaSeq, search::ByteSubstringMut};
use arbitrary::{Arbitrary, Result, Unstructured};

// TODO: This needs refactor and specs structs once the fasta API refactor
// occurs.

impl<'a> Arbitrary<'a> for FastaSeq {
    /// Generates an arbitrary [`FastaSeq`] record from the given unstructured
    /// data.
    ///
    /// This ensures that the header and sequence do not contain the `>` symbol.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let name = String::arbitrary(u)?.replace('>', " ");
        let mut sequence = Vec::<u8>::arbitrary(u)?;
        sequence.replace_all_bytes(b'>', b' ');

        Ok(FastaSeq { name, sequence })
    }
}
