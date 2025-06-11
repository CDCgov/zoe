use crate::{alignment::Alignment, data::types::cigar::Cigar, prelude::*};

mod merge_pairs;
mod reader;
#[cfg(test)]
mod test;

pub use merge_pairs::*;
pub use reader::*;

// # NOTICE
// We define `index` to be 1-based and `position` to be 0-based to avoid
// off-by-one errors and encourage better semantics

/// Struct holding the data for a single
/// [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) record.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SamData {
    /// Query name.
    pub qname: String,
    /// SAM flag: strandedness, etc.
    pub flag:  u16,
    /// Reference name.
    pub rname: String,
    /// The 1-based position in the reference to which the start of the query
    /// aligns. This excludes clipped bases.
    pub pos:   usize,
    /// Mystical map quality value.
    pub mapq:  u8,
    /// Old style cigar format that does not include match and mismatch as
    /// separate values.
    pub cigar: Cigar,
    /// Reference name of the mate / next read. Currently not implemented and
    /// set to `*`.
    rnext:     char,
    /// Position of the mate / next read. Currently not implemented and set to
    /// `0`.
    pnext:     u32,
    /// So-called "observed template length." Currently not implemented and
    /// always set to `0`.
    tlen:      i32,
    /// Query sequence.
    pub seq:   Nucleotides,
    /// Query quality scores in ASCII-encoded format with Phred Quality of +33.
    pub qual:  QualityScores,
}

impl SamData {
    /// Constructs a new [`SamData`] record from the corresponding fields.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        qname: String, flag: u16, rname: String, pos: usize, mapq: u8, cigar: Cigar, seq: Nucleotides, qual: QualityScores,
    ) -> Self {
        SamData {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext: '*',
            pnext: 0,
            tlen: 0,
            seq,
            qual,
        }
    }

    /// Creates a new unmapped [`SamData`] record.
    ///
    /// The sequence and quality fields are set to `*`, `POS` is set to 0,
    /// `MAPQ` is set to 255, and the CIGAR string is empty.
    fn unmapped(qname: &str, rname: &str) -> Self {
        // In the context of an unmapped `SamData` record, this should not be
        // misinterpreted
        let seq = Nucleotides::from(b"*");
        let qual = QualityScores::try_from(b"*").unwrap();
        Self::new(qname.to_string(), 4, rname.to_string(), 0, 255, Cigar::new(), seq, qual)
    }

    /// Constructs a new [`SamData`] record from an [`Alignment`] struct as well
    /// as the other provided fields.
    #[inline]
    #[must_use]
    pub fn from_alignment<T>(
        alignment: Alignment<T>, qname: String, flag: u16, rname: String, mapq: u8, seq: Nucleotides, qual: QualityScores,
    ) -> Self {
        // Both SAM and Alignment exclude clipped bases when reporting
        // positions, so we just need to adjust to 1-based
        let pos = alignment.ref_range.start + 1;
        let cigar = alignment.cigar;
        SamData {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext: '*',
            pnext: 0,
            tlen: 0,
            seq,
            qual,
        }
    }

    /// Tests if the [`SamData`] is unmapped.
    ///
    /// A record is considered unmapped if either the `flag` field has 0x4 set,
    /// or if `cigar` has a match length of 0.
    #[inline]
    #[must_use]
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0 || self.cigar.match_length() == 0
    }
}

impl std::fmt::Display for SamData {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let SamData {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
        } = self;

        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
        )
    }
}
