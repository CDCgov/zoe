//! A module for reading and manipulating
//! [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files. Provides some
//! special-case functions used by [IRMA](https://wonder.cdc.gov/amd/flu/irma/).

use crate::{
    data::types::cigar::{Cigar, CigarView},
    prelude::*,
};
use std::hash::Hash;

mod convert;
mod display_sam_fields;
mod flag;
mod get_sam_fields;
mod optional_fields;
mod reader;
mod sort_traits;
mod std_traits;
mod view_traits;

pub use display_sam_fields::*;
pub use flag::*;
pub use get_sam_fields::*;
pub use optional_fields::*;
pub use reader::*;
pub use sort_traits::SamDataSort;

#[cfg(feature = "dev-generic-sam-reader")]
mod generic_reader;
#[cfg(feature = "dev-generic-sam-reader")]
pub use generic_reader::*;

#[cfg(test)]
mod test;

// # NOTICE
// We define `index` to be 1-based and `position` to be 0-based to avoid
// off-by-one errors and encourage better semantics

/// Struct holding the data for a single
/// [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) record.
#[derive(Clone, Debug)]
pub struct SamData {
    /// The name of the query.
    ///
    /// If the query name is not available (represented as `*` in the SAM file),
    /// an empty string is stored.
    pub qname:      String,
    /// The bitflags for the SAM record, including the strand and whether it was
    /// unmapped.
    pub flag:       Flag,
    /// The name of the reference.
    ///
    /// If the reference name is not available (represented as `*` in the SAM
    /// file), an empty string is stored.
    pub rname:      String,
    /// The 1-based position in the reference where the query starts to align.
    /// This excludes clipped bases.
    pub pos:        usize,
    /// The quality of the mapping.
    ///
    /// Traditionally this holds $-10\log_{10}P(\text{mapping position is
    /// wrong})$. 255 indicates that the mapping quality is unavailable.
    pub mapq:       u8,
    /// The CIGAR string holding the alignment.
    pub cigar:      Cigar,
    /// Reference name of the mate / next read. Currently not implemented and
    /// set to `*`.
    rnext:          char,
    /// Position of the mate / next read. Currently not implemented and set to
    /// `0`.
    pnext:          u32,
    /// So-called "observed template length." Currently not implemented and
    /// always set to `0`.
    tlen:           i32,
    /// The query sequence as stored in the SAM file.
    ///
    /// If the sequence is not available (represented as `*` in the SAM file),
    /// an empty [`Nucleotides`] is stored.
    pub seq:        Nucleotides,
    /// The ASCII-encoded quality scores corresponding to the sequence.
    ///
    /// This is meant to hold the Phred quality of with an offset of `+33`. If
    /// the quality scores are not available (represented as `*` in the SAM
    /// file), an empty [`QualityScores`] is stored.
    pub qual:       QualityScores,
    /// Optional fields which can be lazily parsed and accessed.
    pub opt_fields: SamOptRaw,
}

impl Default for SamData {
    #[inline]
    fn default() -> Self {
        Self {
            qname:      String::new(),
            flag:       Flag::default(),
            rname:      String::new(),
            pos:        0,
            mapq:       255,
            cigar:      Cigar::new(),
            rnext:      '*',
            pnext:      0,
            tlen:       0,
            seq:        Nucleotides::new(),
            qual:       QualityScores::new(),
            opt_fields: SamOptRaw::new(),
        }
    }
}

impl PartialEq for SamData {
    /// Tests for `self` and `other` values to be equal, and is used by `==`.
    /// Note that this implementation ignores the `opt_fields` field, which
    /// contains optional SAM values.
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.qname == other.qname
            && self.flag == other.flag
            && self.rname == other.rname
            && self.pos == other.pos
            && self.mapq == other.mapq
            && self.cigar == other.cigar
            && self.rnext == other.rnext
            && self.pnext == other.pnext
            && self.tlen == other.tlen
            && self.seq == other.seq
            && self.qual == other.qual
    }
}

impl Eq for SamData {}

impl Hash for SamData {
    /// Feeds this value into the given `Hasher`. Note that this implementation
    /// ignores the `opt_fields` field, which contains optional SAM values.
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.qname.hash(state);
        self.flag.hash(state);
        self.rname.hash(state);
        self.pos.hash(state);
        self.mapq.hash(state);
        self.cigar.hash(state);
        self.rnext.hash(state);
        self.pnext.hash(state);
        self.tlen.hash(state);
        self.seq.hash(state);
        self.qual.hash(state);
    }
}

/// A view of a [`SamData`] record, where sequence and string types are views
/// (and primitive types are copied).
///
/// See [Views](crate::data#views) for more details. This struct is primarily
/// used for displaying SAM data without requiring ownership.
#[derive(Clone, Debug)]
pub struct SamDataView<'a> {
    /// Query name.
    pub qname:      &'a str,
    /// SAM flag: strandedness, etc.
    pub flag:       Flag,
    /// Reference name.
    pub rname:      &'a str,
    /// The 1-based position in the reference to which the start of the query
    /// aligns. This excludes clipped bases.
    pub pos:        usize,
    /// Mystical map quality value.
    pub mapq:       u8,
    /// Old style cigar format that does not include match and mismatch as
    /// separate values.
    pub cigar:      CigarView<'a>,
    /// Reference name of the mate / next read. Currently not implemented and
    /// set to `*`.
    rnext:          char,
    /// Position of the mate / next read. Currently not implemented and set to
    /// `0`.
    pnext:          u32,
    /// So-called "observed template length." Currently not implemented and
    /// always set to `0`.
    tlen:           i32,
    /// Query sequence.
    pub seq:        NucleotidesView<'a>,
    /// Query quality scores in ASCII-encoded format with Phred Quality of +33.
    pub qual:       QualityScoresView<'a>,
    /// Optional fields which can be lazily parsed and accessed.
    pub opt_fields: SamOptRawView<'a>,
}

impl Default for SamDataView<'_> {
    #[inline]
    fn default() -> Self {
        Self {
            qname:      "",
            flag:       Flag::default(),
            rname:      "",
            pos:        0,
            mapq:       255,
            cigar:      CigarView::new(),
            rnext:      '*',
            pnext:      0,
            tlen:       0,
            seq:        NucleotidesView::new(),
            qual:       QualityScoresView::new(),
            opt_fields: SamOptRawView::new(),
        }
    }
}

impl PartialEq for SamDataView<'_> {
    /// Tests for `self` and `other` values to be equal, and is used by `==`.
    /// Note that this implementation ignores the `opt_fields` field, which
    /// contains optional SAM values.
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.qname == other.qname
            && self.flag == other.flag
            && self.rname == other.rname
            && self.pos == other.pos
            && self.mapq == other.mapq
            && self.cigar == other.cigar
            && self.rnext == other.rnext
            && self.pnext == other.pnext
            && self.tlen == other.tlen
            && self.seq == other.seq
            && self.qual == other.qual
    }
}

impl Eq for SamDataView<'_> {}

impl Hash for SamDataView<'_> {
    /// Feeds this value into the given `Hasher`. Note that this implementation
    /// ignores the `opt_fields` field, which contains optional SAM values.
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.qname.hash(state);
        self.flag.hash(state);
        self.rname.hash(state);
        self.pos.hash(state);
        self.mapq.hash(state);
        self.cigar.hash(state);
        self.rnext.hash(state);
        self.pnext.hash(state);
        self.tlen.hash(state);
        self.seq.hash(state);
        self.qual.hash(state);
    }
}

impl SamData {
    /// Constructs a new [`SamData`] record from the corresponding fields.
    ///
    /// `opt_fields` is set to empty.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        qname: String, flag: Flag, rname: String, pos: usize, mapq: u8, cigar: Cigar, seq: Nucleotides, qual: QualityScores,
    ) -> Self {
        // Validate incoming arguments to ensure they are not missing
        let qname = if is_missing_sam_field(&qname) { String::new() } else { qname };
        let rname = if is_missing_sam_field(&rname) { String::new() } else { rname };
        let seq = if is_missing_sam_field(&seq) { Nucleotides::new() } else { seq };
        let qual = if is_missing_sam_field(&qual) {
            QualityScores::new()
        } else {
            qual
        };

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
            opt_fields: SamOptRaw::new(),
        }
    }

    /// Creates a new unmapped [`SamData`] record.
    ///
    /// POS is set to 0, MAPQ is set to 255, and SEQ, QUAL, and CIGAR are empty
    /// (displayed as `*`).
    #[inline]
    #[must_use]
    pub fn unmapped(qname: &str, rname: &str) -> Self {
        Self::new(
            qname.to_string(),
            Flag::UNMAPPED,
            rname.to_string(),
            0,
            255,
            Cigar::new(),
            Nucleotides::new(),
            QualityScores::new(),
        )
    }
}

impl<'a> SamDataView<'a> {
    /// Constructs a new [`SamDataView`] record from the corresponding fields.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        qname: &'a str, flag: Flag, rname: &'a str, pos: usize, mapq: u8, cigar: CigarView<'a>, seq: NucleotidesView<'a>,
        qual: QualityScoresView<'a>,
    ) -> Self {
        // Validate incoming arguments to ensure they are not missing
        let qname = if is_missing_sam_field(qname) { "" } else { qname };
        let rname = if is_missing_sam_field(rname) { "" } else { rname };
        let seq = if is_missing_sam_field(seq) {
            NucleotidesView::new()
        } else {
            seq
        };
        let qual = if is_missing_sam_field(qual) {
            QualityScoresView::new()
        } else {
            qual
        };

        SamDataView {
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
            opt_fields: SamOptRawView::new(),
        }
    }

    /// Creates a new unmapped [`SamDataView`] record.
    ///
    /// POS is set to 0, MAPQ is set to 255, and SEQ, QUAL, and CIGAR are empty
    /// (displayed as `*`).
    #[inline]
    #[must_use]
    pub fn unmapped(qname: &'a str, rname: &'a str) -> Self {
        Self::new(
            qname,
            Flag::UNMAPPED,
            rname,
            0,
            255,
            CigarView::new(),
            NucleotidesView::new(),
            QualityScoresView::new(),
        )
    }
}

/// Returns whether a SAM field should be treated as missing.
///
/// This supports QNAME, RNAME, CIGAR (after converting to bytes), SEQ, and
/// QUAL.
///
/// *Zoe* accepts both the SAM sentinel `*` and an empty byte sequence as
/// missing for these fields.
pub(crate) fn is_missing_sam_field(field: impl AsRef<[u8]>) -> bool {
    let field = field.as_ref();
    field.is_empty() || field == b"*"
}
