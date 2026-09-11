//! A module for reading and manipulating
//! [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files. Provides some
//! special-case functions used by [IRMA](https://wonder.cdc.gov/amd/flu/irma/).

use crate::{
    alignment::{Alignment, AlignmentStates, MaybeAligned, NextCiglet},
    data::{
        cigar::LenInAlignment,
        types::cigar::{Cigar, CigarView, CigarViewMut},
    },
    math::AnyInt,
    prelude::*,
};
use std::{
    fmt::{Display, Formatter},
    hash::Hash,
};

mod optional_fields;
mod reader;
mod sort_traits;
mod std_traits;
mod view_traits;

pub use optional_fields::*;
pub use reader::*;
pub use sort_traits::SamDataSort;

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
    pub flag:       u16,
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
    pub flag:       u16,
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

/// A mutable view of a [`SamData`] record, where sequence and string types are
/// views (and primitive types are copied).
///
/// See [Views](crate::data#views) for more details. This struct is primarily
/// used for displaying SAM data without requiring ownership.
#[derive(Eq, PartialEq, Hash, Debug)]
#[deprecated(
    since = "0.0.33",
    note = "consider using an immutable view or a custom struct instead. This struct will be removed in v0.0.35. Open an issue with a use-case if this struct is required"
)]
pub struct SamDataViewMut<'a> {
    /// Query name.
    pub qname: &'a mut String,
    /// SAM flag: strandedness, etc.
    pub flag:  u16,
    /// Reference name.
    pub rname: &'a mut String,
    /// The 1-based position in the reference to which the start of the query
    /// aligns. This excludes clipped bases.
    pub pos:   usize,
    /// Mystical map quality value.
    pub mapq:  u8,
    /// Old style cigar format that does not include match and mismatch as
    /// separate values.
    pub cigar: CigarViewMut<'a>,
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
    pub seq:   NucleotidesViewMut<'a>,
    /// Query quality scores in ASCII-encoded format with Phred Quality of +33.
    pub qual:  QualityScoresViewMut<'a>,
}

impl SamData {
    /// Constructs a new [`SamData`] record from the corresponding fields.
    ///
    /// `opt_fields` is set to empty.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
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
        Self::new(qname.to_string(), 4, rname.to_string(), 0, 255, Cigar::new(), seq, qual)
    }

    /// Constructs a new [`SamData`] record from an [`Alignment`] struct as well
    /// as the other provided fields. The score is included as a field under the
    /// tag `AS`.
    ///
    /// For the opposite transformation, see [`SamData::to_alignment`].
    #[inline]
    #[must_use]
    pub fn from_alignment<T: AnyInt + Into<i64>>(
        alignment: &Alignment<T>, qname: String, flag: u16, rname: String, mapq: u8, seq: Nucleotides, qual: QualityScores,
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
    /// Any hard clipped bases in the query are not included in the resulting
    /// `query_range` or `query_len`. The `query_len` field is equal to the
    /// length of the incoming `seq` field in the SAM record, and `query_range`
    /// represents the entire `seq` except for soft clipping.
    ///
    /// The `score` must be provided as an argument. If the score is present
    /// under the `AS` TAG, then it can be retrieved with:
    ///
    /// ```
    /// # use zoe::{
    /// #     data::{cigar::Cigar, sam::{SamOptValue, SamData}},
    /// #     prelude::{Nucleotides, QualityScores},
    /// # };
    /// #
    /// # let mut sam_data = SamData::new(
    /// #     String::new(),
    /// #     0,
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
    /// For any mapped read:
    ///
    /// - The `seq` field of the record must be populated (i.e., not `*`).
    /// - The CIGAR operations must be among `M, I, D, N, S, H, P, X, =`
    /// - Every operation in the CIGAR string must have a preceding increment
    /// - Every increment must be followed by an operation
    /// - The increment for each operation must be non-zero and less than
    ///   [`usize::MAX`]
    ///
    /// ## Validity
    ///
    /// For any mapped read, the length of the `seq` field should equal the sum
    /// of the increments for the operations `M`, `I`, `S`, `=`, and `X`. Note
    /// that this includes soft clipped regions but not hard clipped.
    #[inline]
    pub fn to_alignment<T: AnyInt>(&self, score: T, ref_len: usize) -> std::io::Result<MaybeAligned<Alignment<T>>> {
        if self.is_unmapped() {
            return Ok(MaybeAligned::Unmapped);
        }

        if is_missing_sam_field(&self.seq) {
            return Err(std::io::Error::other(
                "The seq field in the SAM data record was not populated.",
            ));
        }

        let query_len = self.seq.len();

        let ref_range_start = self.pos - 1;
        let ref_range_end = ref_range_start + self.cigar.ref_len_in_alignment();
        let ref_range = ref_range_start..ref_range_end;

        let mut ciglets = self.cigar.iter();
        ciglets.next_ciglet_if_op(|op| op == b'H');
        let soft_clipping_front = ciglets.next_ciglet_if_op(|op| op == b'S').map_or(0, |c| c.inc);
        ciglets.next_ciglet_back_if_op(|op| op == b'H');
        let soft_clipping_back = ciglets.next_ciglet_back_if_op(|op| op == b'S').map_or(0, |c| c.inc);

        let soft_clipping = soft_clipping_front + soft_clipping_back;

        let query_range_start = soft_clipping_front;
        let query_range_end = query_range_start + (query_len - soft_clipping);
        let query_range = query_range_start..query_range_end;

        let states = AlignmentStates::try_from(&self.cigar).map_err(std::io::Error::other)?;

        Ok(MaybeAligned::Some(Alignment {
            score,
            ref_range,
            query_range,
            states,
            ref_len,
            query_len,
        }))
    }

    /// Tests if the [`SamData`] is unmapped.
    ///
    /// A record is considered unmapped if either the `flag` field has 0x4 set,
    /// or if `cigar` has a match length of 0.
    #[inline]
    #[must_use]
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0 || self.cigar.ref_len_in_alignment() == 0
    }
}

impl<'a> SamDataView<'a> {
    /// Constructs a new [`SamDataView`] record from the corresponding fields.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        qname: &'a str, flag: u16, rname: &'a str, pos: usize, mapq: u8, cigar: CigarView<'a>, seq: NucleotidesView<'a>,
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
        Self::new(qname, 4, rname, 0, 255, CigarView::new(), seq, qual)
    }
}

#[allow(deprecated)]
impl<'a> SamDataViewMut<'a> {
    /// Constructs a new [`SamDataViewMut`] record from the corresponding
    /// fields.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        qname: &'a mut String, flag: u16, rname: &'a mut String, pos: usize, mapq: u8, cigar: CigarViewMut<'a>,
        seq: NucleotidesViewMut<'a>, qual: QualityScoresViewMut<'a>,
    ) -> Self {
        SamDataViewMut {
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
}

/// Returns whether a SAM `SEQ` or `QUAL` field should be treated as missing.
///
/// *Zoe* accepts both the SAM sentinel `*` and an empty byte sequence as
/// missing for these fields.
pub(crate) fn is_missing_sam_field(field: impl AsRef<[u8]>) -> bool {
    let field = field.as_ref();
    field.is_empty() || field == b"*"
}
