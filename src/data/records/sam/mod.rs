//! A module for reading and manipulating
//! [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files. Provides some
//! special-case functions used by [IRMA](https://wonder.cdc.gov/amd/flu/irma/).

use crate::{
    alignment::{AlignOp, Alignment, AlignmentStates, MaybeAligned, NextCiglet},
    data::{
        cigar::LenInAlignment,
        types::cigar::{Cigar, CigarView},
    },
    math::AnyInt,
    prelude::*,
};
use std::hash::Hash;

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
    /// Query name.
    pub qname:      String,
    /// SAM flag: strandedness, etc.
    pub flag:       Flag,
    /// Reference name.
    pub rname:      String,
    /// The 1-based position in the reference to which the start of the query
    /// aligns. This excludes clipped bases.
    pub pos:        usize,
    /// Mystical map quality value.
    pub mapq:       u8,
    /// Old style cigar format that does not include match and mismatch as
    /// separate values.
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
    /// Query sequence.
    pub seq:        Nucleotides,
    /// Query quality scores in ASCII-encoded format with Phred quality of +33.
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
    /// The sequence and quality fields are set to `*`, `POS` is set to 0,
    /// `MAPQ` is set to 255, and the CIGAR string is empty.
    #[inline]
    #[must_use]
    pub fn unmapped(qname: &str, rname: &str) -> Self {
        // In the context of an unmapped `SamData` record, this should not be
        // misinterpreted
        let seq = Nucleotides::from(b"*");
        // Safety: * is graphic ascii
        let qual = unsafe { QualityScores::from_vec_unchecked(b"*".to_vec()) };
        Self::new(
            qname.to_string(),
            Flag::UNMAPPED,
            rname.to_string(),
            0,
            255,
            Cigar::new(),
            seq,
            qual,
        )
    }

    /// Constructs a new [`SamData`] record from an [`Alignment`] struct as well
    /// as the other provided fields. The score is included as a field under the
    /// tag `AS`.
    ///
    /// For the opposite transformation, see [`SamData::to_alignment`].
    #[inline]
    #[must_use]
    pub fn from_alignment<T: AnyInt + Into<i64>>(
        alignment: &Alignment<T>, qname: String, flag: Flag, rname: String, mapq: u8, seq: Nucleotides, qual: QualityScores,
    ) -> Self {
        // Both SAM and Alignment exclude clipped bases when reporting
        // positions, so we just need to adjust to 1-based
        let pos = alignment.ref_range.start + 1;
        let cigar = alignment.states.to_cigar_unchecked();
        let opt_fields = SamOptRaw::new_with_score(alignment.score);
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
            opt_fields,
        }
    }

    /// Converts the [`SamData`] record into an [`Alignment`] struct (wrapped in
    /// [`MaybeAligned`]).
    ///
    /// If the entire CIGAR string is empty or consists solely of clipping, or
    /// if [`SamData::flag`] has the unmapped flag set, then
    /// [`MaybeAligned::Unmapped`] is returned.
    ///
    /// Any hard clipped bases in the query are not included in the resulting
    /// `query_range` or `query_len`. `query_range` represents the entire `SEQ`
    /// except for soft clipping, while `query_len` represents the length with
    /// soft clipping (which must also equal the length of `SEQ` if present).
    ///
    /// The `score` must be provided as an argument. A ZST like `()` can also be
    /// passed to avoid storing the score. If the score is present under the
    /// `AS` TAG, then it can be retrieved with:
    ///
    /// ```
    /// # use zoe::{
    /// #     data::{cigar::Cigar, sam::{Flag, SamOptValue, SamData}},
    /// #     prelude::{Nucleotides, QualityScores},
    /// # };
    /// #
    /// # let mut sam_data = SamData::new(
    /// #     String::new(),
    /// #     Flag::default(),
    /// #     String::new(),
    /// #     0,
    /// #     255,
    /// #     Cigar::new(),
    /// #     Nucleotides::new(),
    /// #     QualityScores::new(),
    /// # );
    /// #
    /// # sam_data.opt_fields.push("AS", &SamOptValue::Int(0));
    /// #
    /// let score = sam_data
    ///     .opt_fields
    ///     .get("AS")
    ///     .expect("The optional fields must be formatted properly")
    ///     .expect("The score TAG must be present")
    ///     .int()
    ///     .expect("The score VALUE should be an integer");
    /// ```
    ///
    /// For the opposite transformation, see [`SamData::from_alignment`].
    ///
    /// ## Errors
    ///
    /// Unmapped reads are those where the FLAG has the unmapped bit set, the
    /// CIGAR field is not populated, or the CIGAR field is solely clipping.
    /// These are checked in this order, and [`MaybeAligned::Unmapped`] is
    /// returned immediately when a condition is found true. During this
    /// process, an error is returned if the CIGAR field must be read and is
    /// malformed.
    ///
    /// For any mapped read, the following additional errors could occur:
    ///
    /// - The CIGAR field must have no clipping in the interior of the string
    /// - The number of reference bases consumed by the CIGAR string must not
    ///   overflow [`usize`], and similarly for the query
    /// - If the SEQ field is populated, then the sequence length must match the
    ///   number of query bases consumed by the CIGAR string
    /// - The POS field must be populated (i.e., not zero)
    /// - The POS and CIGAR field must yield an ending reference index that does
    ///   not exceed `ref_len`
    #[inline]
    pub fn to_alignment<T>(&self, score: T, ref_len: usize) -> std::io::Result<MaybeAligned<Alignment<T>>> {
        // Handle FLAG
        if self.flag.is_unmapped() {
            return Ok(MaybeAligned::Unmapped);
        }

        // Parse CIGAR to AlignmentStates
        let states = AlignmentStates::try_from(&self.cigar).map_err(std::io::Error::other)?;
        if states.is_empty() {
            return Ok(MaybeAligned::Unmapped);
        }

        // Compute the start of the aligned query by examining clipping in CIGAR
        let mut ciglets = states.as_slice();
        ciglets.next_ciglet_if_op(|op| op == AlignOp::H);
        let soft_clipping_front = ciglets.next_ciglet_if_op(|op| op == AlignOp::S).map_or(0, |c| c.inc);
        ciglets.next_ciglet_back_if_op(|op| op == AlignOp::H);
        let soft_clipping_back = ciglets.next_ciglet_back_if_op(|op| op == AlignOp::S).map_or(0, |c| c.inc);
        if ciglets.is_empty() {
            return Ok(MaybeAligned::Unmapped);
        }
        if ciglets.iter().any(|c| matches!(c.op, AlignOp::S | AlignOp::H)) {
            return Err(std::io::Error::other("Found clipping in the interior of the CIGAR string"));
        }

        let query_range_start = soft_clipping_front;

        // Compute number of residues consumed in CIGAR
        let Some(aligned_ref_len) = states.ref_len_in_alignment_checked() else {
            return Err(std::io::Error::other(
                "The number of reference residues consumed by the CIGAR string overflowed usize",
            ));
        };
        let Some(query_len) = states.query_len_in_alignment_checked() else {
            return Err(std::io::Error::other(
                "The number of query residues consumed by the CIGAR string overflowed usize",
            ));
        };

        // Confirm CIGAR is compatible with SEQ, if SEQ is present
        if !is_missing_sam_field(&self.seq) && self.seq.len() != query_len {
            return Err(std::io::Error::other(format!(
                "The CIGAR field consumed {query_len} residues of the query, but the length of SEQ was different {seq_len}",
                seq_len = self.seq.len()
            )));
        }

        // Parse POS
        let Some(ref_range_start) = self.pos.checked_sub(1) else {
            return Err(std::io::Error::other("The POS field was not populated"));
        };

        // Compute and validate the index for the end of the reference
        let Some(ref_range_end) = aligned_ref_len.checked_add(ref_range_start).filter(|&end| end <= ref_len) else {
            return Err(std::io::Error::other(format!(
                "The alignment implied by POS and CIGAR extends beyond ref_len ({ref_len})"
            )));
        };

        // Compute the end of the aligned query
        let query_range_end = query_len - soft_clipping_back;

        // Form ranges
        let ref_range = ref_range_start..ref_range_end;
        let query_range = query_range_start..query_range_end;

        Ok(MaybeAligned::Some(Alignment {
            score,
            ref_range,
            query_range,
            states,
            ref_len,
            query_len,
        }))
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
    /// The sequence and quality fields are set to `*`, `POS` is set to 0,
    /// `MAPQ` is set to 255, and the CIGAR string is empty.
    #[inline]
    #[must_use]
    pub fn unmapped(qname: &'a str, rname: &'a str) -> Self {
        // In the context of an unmapped `SamData` record, this should not be
        // misinterpreted
        let seq = NucleotidesView::from(b"*");
        // Safety: * is graphic ascii
        let qual = unsafe { QualityScoresView::from_bytes_unchecked(b"*") };
        Self::new(qname, Flag::UNMAPPED, rname, 0, 255, CigarView::new(), seq, qual)
    }
}

/// Returns whether a SAM `SEQ` or `QUAL` field should be treated as missing.
///
/// *Zoe* accepts both the SAM sentinel `*` and an empty byte sequence as
/// missing for these fields.
pub(crate) fn is_missing_sam_field(field: impl AsRef<[u8]>) -> bool {
    let field = field.as_ref();
    field.is_empty() || field == b"*"
}
