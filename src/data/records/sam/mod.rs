use crate::{
    alignment::{Alignment, AlignmentStates, MaybeAligned, StatesSequence},
    data::{
        cigar::LenInAlignment,
        err::ResultWithErrorContext,
        types::cigar::{Cigar, CigarView, CigarViewMut},
    },
    iter_utils::ProcessResultsExt,
    math::AnyInt,
    prelude::*,
};
use std::{
    fmt::{Display, Formatter},
    hash::Hash,
};

mod reader;
mod sort_traits;
mod std_traits;
mod view_traits;

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
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SamDataView<'a> {
    /// Query name.
    pub qname: &'a str,
    /// SAM flag: strandedness, etc.
    pub flag:  u16,
    /// Reference name.
    pub rname: &'a str,
    /// The 1-based position in the reference to which the start of the query
    /// aligns. This excludes clipped bases.
    pub pos:   usize,
    /// Mystical map quality value.
    pub mapq:  u8,
    /// Old style cigar format that does not include match and mismatch as
    /// separate values.
    pub cigar: CigarView<'a>,
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
    pub seq:   NucleotidesView<'a>,
    /// Query quality scores in ASCII-encoded format with Phred Quality of +33.
    pub qual:  QualityScoresView<'a>,
}

/// A mutable view of a [`SamData`] record, where sequence and string types are
/// views (and primitive types are copied).
///
/// See [Views](crate::data#views) for more details. This struct is primarily
/// used for displaying SAM data without requiring ownership.
#[derive(Eq, PartialEq, Hash, Debug)]
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

        if self.seq == Nucleotides::from(b"*") {
            return Err(std::io::Error::other(
                "The seq field in the SAM data record was not populated.",
            ));
        }

        let query_len = self.seq.len();

        let ref_range_start = self.pos - 1;
        let ref_range_end = ref_range_start + self.cigar.ref_len_in_alignment();
        let ref_range = ref_range_start..ref_range_end;

        let mut ciglets = self.cigar.iter();
        ciglets.next_if_op(|op| op == b'H');
        let soft_clipping_front = ciglets.next_if_op(|op| op == b'S').map_or(0, |c| c.inc);
        ciglets.next_back_if_op(|op| op == b'H');
        let soft_clipping_back = ciglets.next_back_if_op(|op| op == b'S').map_or(0, |c| c.inc);

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

/// Any optional fields stored in a SAM record, lazily parsed on an as-needed
/// basis.
///
/// Each optional field consists of a tag, value type, and value.
#[derive(Clone, Debug, Default)]
pub struct SamOptRaw(Vec<String>);

impl SamOptRaw {
    /// Returns an empty collection of optional fields.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        SamOptRaw(Vec::new())
    }

    /// Returns [`SamOptRaw`] containing just a single field with the alignment
    /// score.
    ///
    /// The score is represented as an integer using the `AS` tag.
    #[inline]
    #[must_use]
    pub fn new_with_score<T: AnyInt + Into<i64>>(score: T) -> Self {
        let mut inner = Vec::with_capacity(1);
        inner.push(format!("AS:i:{score}", score = score.into()));
        SamOptRaw(inner)
    }

    /// Returns whether the optional data is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of optional fields present.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Provides an iterator over the optional fields present (the tag names and
    /// parsed values).
    ///
    /// ## Limitations
    ///
    /// This iterator parses the fields lazily. If the [`SamOptRaw`] struct will
    /// be iterated over many times, consider parsing the fields once and
    /// collecting them.
    ///
    /// ## Errors
    ///
    /// The field in the SAM record must be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// or `B`. `VALUE` must successfully parse into the corresponding type.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = std::io::Result<SamOptField>> {
        self.0.iter().map(|field| {
            let inv_opt_err_msg = || std::io::Error::other(format!("Invalid optional field {field}"));

            let (tag_text, rest) = field.split_once(':').ok_or_else(inv_opt_err_msg)?;
            let (type_text, string_value) = rest.split_once(':').ok_or_else(inv_opt_err_msg)?;

            let tag = SamOptField::parse_tag(tag_text)?;
            let type_code = SamOptField::parse_type(type_text)?;
            let opt_field = SamOptField::parse_value(tag, type_code, string_value)
                .with_context(format!("Failed to parse field '{field}'"))?;
            Ok(opt_field)
        })
    }

    /// Returns the optional data for the provided tag, if it is present.
    ///
    /// ## Limitations
    ///
    /// This struct parses fields lazily and will repeat the computations each
    /// time [`get`] is called. The function runs in $O(n)$ time where $n$ is
    /// the number of fields. Validation of the field format is only performed
    /// where necessary.
    ///
    /// ## Errors
    ///
    /// The field in the SAM record must be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// `B`. `VALUE` must successfully parse into the corresponding type.
    ///
    /// [`get`]: SamOptRaw::get
    pub fn get(&self, tag: &str) -> std::io::Result<Option<SamOptField>> {
        for field in &self.0 {
            let inv_opt_err_msg = || std::io::Error::other(format!("Invalid optional field {field}"));
            let (this_tag, rest) = field.split_once(':').ok_or_else(inv_opt_err_msg)?;

            if this_tag == tag {
                let (type_text, string_value) = rest.split_once(':').ok_or_else(inv_opt_err_msg)?;

                let tag = SamOptField::parse_tag(this_tag)?;
                let type_code = SamOptField::parse_type(type_text)?;

                let opt_field = match SamOptField::parse_value(tag, type_code, string_value) {
                    Ok(opt_field) => opt_field,
                    Err(e) => {
                        return Err(std::io::Error::other(format!(
                            "Failed to parse field '{field}' due to error: {e}"
                        )));
                    }
                };
                return Ok(Some(opt_field));
            }
        }
        Ok(None)
    }

    /// Adds an optional field to the [`SamOptRaw`] struct.
    ///
    /// ## Validity
    ///
    /// The tag name being pushed should not already be present in `self`.
    #[inline]
    pub fn push(&mut self, tag: &str, data: &SamOptValue) {
        self.0.push(format!("{tag}:{data}"));
    }
}

impl FromIterator<String> for SamOptRaw {
    /// Collects an iterator of strings into a [`SamOptRaw`] collection (each
    /// following the SAM file format for an optional field).
    ///
    /// ## Validity
    ///
    /// Each string should conform to the SAM file format for an optional field.
    /// Specifically, each field should be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// or `B`. `VALUE` must successfully parse into the corresponding type.
    /// Furthermore, the tags should be unique.
    #[inline]
    fn from_iter<T: IntoIterator<Item = String>>(iter: T) -> Self {
        SamOptRaw(Vec::from_iter(iter))
    }
}

/// A parsed optional field in the SAM file format.
#[derive(Clone, Debug)]
pub struct SamOptField {
    /// The tag of the optional SAM field.
    pub tag:   [u8; 2],
    /// The value of the optional SAM field.
    pub value: SamOptValue,
}

impl SamOptField {
    /// Parses the tag for the optional SAM field from a string slice.
    fn parse_tag(tag: &str) -> std::io::Result<[u8; 2]> {
        let bytes = tag.as_bytes();
        if bytes.len() != 2 || !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphanumeric() {
            return Err(std::io::Error::other(format!("Invalid SAM optional tag: {tag}")));
        }
        Ok([bytes[0], bytes[1]])
    }

    /// Parses the type for the optional SAM field from a string slice.
    fn parse_type(type_text: &str) -> std::io::Result<char> {
        let mut type_chars = type_text.chars();
        let Some(typ) = type_chars.next() else {
            return Err(std::io::Error::other("Missing optional field type"));
        };
        if type_chars.next().is_some() {
            return Err(std::io::Error::other(format!("Invalid optional field type {type_text}")));
        }

        Ok(typ)
    }

    /// Parses a [`SamOptField`] from a tag, type, and value (as a string).
    ///
    /// ## Errors
    ///
    /// `type_code` must contain a valid character (`A`, `i`, `f`, `Z`, `H`, or
    /// `B`). The `string_value` must successfully parse into the corresponding
    /// type.
    fn parse_value(tag: [u8; 2], type_code: char, string_value: &str) -> std::io::Result<SamOptField> {
        match type_code {
            'A' => {
                let mut chars = string_value.chars();
                let Some(c) = chars.next() else {
                    return Err(std::io::Error::other("'A' field has empty value"));
                };
                if chars.next().is_some() {
                    return Err(std::io::Error::other("'A' field must contain exactly one character"));
                }
                if !c.is_ascii() {
                    return Err(std::io::Error::other("'A' field must be ASCII"));
                }
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::Char(c as u8),
                })
            }
            'i' => {
                let parsed = string_value.parse::<i64>().with_context("Error parsing 'i' field")?;
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::Int(parsed),
                })
            }
            'f' => {
                let parsed = string_value.parse::<f32>().with_context("Error parsing 'f' field")?;
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::Float(parsed),
                })
            }
            'Z' => Ok(SamOptField {
                tag,
                value: SamOptValue::String(String::from(string_value)),
            }),
            'H' => {
                if !string_value.len().is_multiple_of(2) {
                    return Err(std::io::Error::other(format!(
                        "'H' field must contain an even number digits. Found {}",
                        string_value.len()
                    )));
                }
                if !string_value.as_bytes().iter().all(u8::is_ascii_hexdigit) {
                    return Err(std::io::Error::other("'H' field must contain hexadecimal digits"));
                }
                Ok(SamOptField {
                    tag,
                    value: SamOptValue::Hex(string_value.to_ascii_uppercase()),
                })
            }
            'B' => Ok(SamOptField {
                tag,
                value: SamOptValue::parse_opt_array(string_value).with_context("Failed to parse 'B' array")?,
            }),
            _ => Err(std::io::Error::other(format!(
                "Unsupported SAM optional field type {type_code}"
            ))),
        }
    }

    /// Returns the stored character from the [`SamOptField`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn char(self) -> Option<u8> {
        match self.value {
            SamOptValue::Char(c) => Some(c),
            _ => None,
        }
    }

    /// Returns the stored integer from the [`SamOptField`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn int(self) -> Option<i64> {
        match self.value {
            SamOptValue::Int(i) => Some(i),
            _ => None,
        }
    }

    /// Returns the stored floating point number from the [`SamOptField`], or
    /// [`None`] if a different variant is present.
    #[inline]
    #[must_use]
    pub fn float(self) -> Option<f32> {
        match self.value {
            SamOptValue::Float(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored string from the [`SamOptField`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn string(self) -> Option<String> {
        match self.value {
            SamOptValue::String(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored hex string from the [`SamOptField`], or [`None`] if a
    /// different variant is present.
    ///
    /// For example, the six-character hex string "1AE301" represents the byte
    /// array `[0x1a, 0xe3, 0x01]`.
    #[inline]
    #[must_use]
    pub fn hex(self) -> Option<String> {
        match self.value {
            SamOptValue::Hex(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored [`OptArray`] from the [`SamOptField`], or [`None`] if
    /// a different variant is present.
    #[inline]
    #[must_use]
    pub fn array(self) -> Option<OptArray> {
        match self.value {
            SamOptValue::Array(f) => Some(f),
            _ => None,
        }
    }
}

/// The value of an optional field (for the SAM file format).
#[derive(Clone, Debug)]
pub enum SamOptValue {
    /// A printable character (type code `A`).
    Char(u8),
    /// A signed integer (type code `i`).
    Int(i64),
    /// A single-precision floating number (type code `f`).
    Float(f32),
    /// A printable string, including space (type code `Z`).
    String(String),
    /// A byte array in the hex format (type code `H`).
    ///
    /// For example, the six-character hex string "1AE301" represents the byte
    /// array `[0x1a, 0xe3, 0x01]`.
    Hex(String),
    // An integer or numeric array (type code `B`).
    Array(OptArray),
}

impl SamOptValue {
    /// Parses the array of optional SAM fields with type `B`.
    ///
    /// ## Errors
    ///
    /// The first letter in the array indicates the type of numbers in the
    /// following comma-separated array. The letter can be one of `c`, `C`, `s`,
    /// `S`, `i`, `I`, or `f`.
    fn parse_opt_array(string_value: &str) -> std::io::Result<Self> {
        let mut pieces = string_value.split(',');
        let Some(subtype) = pieces.next() else {
            return Err(std::io::Error::other("Missing subtype"));
        };

        match subtype {
            "c" => {
                let values = pieces
                    .map(str::parse::<i8>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'c' subtype (`i8`)")?;

                Ok(Self::Array(OptArray::I8(values)))
            }
            "C" => {
                let values = pieces
                    .map(str::parse::<u8>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'C' subtype (`u8`)")?;
                Ok(Self::Array(OptArray::U8(values)))
            }
            "s" => {
                let values = pieces
                    .map(str::parse::<i16>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 's' subtype (`i16`)")?;
                Ok(Self::Array(OptArray::I16(values)))
            }
            "S" => {
                let values = pieces
                    .map(str::parse::<u16>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'S' subtype (`u16`)")?;
                Ok(Self::Array(OptArray::U16(values)))
            }
            "i" => {
                let values = pieces
                    .map(str::parse::<i32>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'i' subtype (`i32`)")?;
                Ok(Self::Array(OptArray::I32(values)))
            }
            "I" => {
                let values = pieces
                    .map(str::parse::<u32>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'I' subtype (`u32`)")?;
                Ok(Self::Array(OptArray::U32(values)))
            }
            "f" => {
                let values = pieces
                    .map(str::parse::<f32>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'f' subtype (`f32`)")?;
                Ok(Self::Array(OptArray::F32(values)))
            }
            _ => Err(std::io::Error::other(format!("Unsupported subtype {subtype}"))),
        }
    }
}

impl Display for SamOptValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SamOptValue::Char(val) => write!(f, "A:{val}", val = *val as char),
            SamOptValue::Int(val) => write!(f, "i:{val}"),
            SamOptValue::Float(val) => write!(f, "f:{val}"),
            SamOptValue::String(val) => write!(f, "Z:{val}"),
            SamOptValue::Hex(val) => write!(f, "H:{val}"),
            SamOptValue::Array(val) => match val {
                OptArray::I8(vec) => OptArray::fmt_opt_array(f, 'c', vec),
                OptArray::U8(vec) => OptArray::fmt_opt_array(f, 'C', vec),
                OptArray::I16(vec) => OptArray::fmt_opt_array(f, 's', vec),
                OptArray::U16(vec) => OptArray::fmt_opt_array(f, 'S', vec),
                OptArray::I32(vec) => OptArray::fmt_opt_array(f, 'i', vec),
                OptArray::U32(vec) => OptArray::fmt_opt_array(f, 'I', vec),
                OptArray::F32(vec) => OptArray::fmt_opt_array(f, 'f', vec),
            },
        }
    }
}

/// The data array for `B` field data.
#[derive(Debug, Clone)]
pub enum OptArray {
    /// Array subtype code `c`.
    I8(Vec<i8>),
    /// Array subtype code `C`.
    U8(Vec<u8>),
    /// Array subtype code `s`.
    I16(Vec<i16>),
    /// Array subtype code `S`.
    U16(Vec<u16>),
    /// Array subtype code `i`.
    I32(Vec<i32>),
    /// Array subtype code `I`.
    U32(Vec<u32>),
    /// Array subtype code `f`.
    F32(Vec<f32>),
}

impl OptArray {
    fn fmt_opt_array<T: Display>(f: &mut Formatter<'_>, arr_type: char, vals: &[T]) -> std::fmt::Result {
        write!(f, "B:{arr_type}")?;

        let mut iter = vals.iter();
        if let Some(first) = iter.next() {
            write!(f, "{first}")?;
        }
        for v in iter {
            write!(f, ",{v}")?;
        }

        Ok(())
    }
}
