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

mod merge_pairs;
mod reader;
mod std_traits;
mod view_traits;

pub use merge_pairs::*;
pub use reader::*;

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
    /// Optional fields which can be lazily parsed/accessed
    pub aux:   SamAuxRaw,
}

impl PartialEq for SamData {
    /// Tests for `self` and `other` values to be equal, and is used by `==`.
    /// Note that this implementation ignores the `aux` field, which contains
    /// optional SAM values.
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
    /// ignores the `aux` field, which contains optional SAM values.
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
            aux: SamAuxRaw::new(),
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
    /// TAG `AS`.
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
        let aux = SamAuxRaw::new_with_score(alignment.score);
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
            aux,
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
    /// #     data::{cigar::Cigar, sam::{SamAuxValue, SamData}},
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
    /// # sam_data.aux.push("AS", &SamAuxValue::Int(0));
    /// #
    /// let score = sam_data
    ///     .aux
    ///     .get("AS")
    ///     .expect("The auxiliary data must be formatted properly")
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

/// A wrapper around an array for holding optional SAM fields.
#[derive(Clone, Debug, Default)]
pub struct SamAuxRaw(Vec<String>);

impl SamAuxRaw {
    /// Returns a new empty [`SamAuxRaw`] vector.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        SamAuxRaw(Vec::new())
    }

    /// Returns a [`SamAuxRaw`] vector with a single field containing the
    /// alignment score.
    ///
    /// The score is represented as an integer using the `AS` TAG.
    #[inline]
    #[must_use]
    pub fn new_with_score<T: AnyInt + Into<i64>>(score: T) -> Self {
        let mut inner = Vec::with_capacity(1);
        inner.push(format!("AS:i:{score}", score = score.into()));
        SamAuxRaw(inner)
    }

    /// Returns whether the [`SamAuxRaw`] vector is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of auxilary fields present.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Provides an iterator over the field names and parsed values.
    ///
    /// ## Limitations
    ///
    /// This iterator parses the fields lazily. If the [`SamAuxRaw`] struct will
    /// be iterated over many times, consider parsing the fields once and
    /// collecting them.
    ///
    /// ## Errors
    ///
    /// The field in the SAM record must be of the form `TAG:TYPE:VALUE`. `TAG`
    /// cannot contain a colon. `TYPE` must be either `A`, `i`, `f`, `Z`, `H`,
    /// or `B`. `VALUE` must successfully parse into the corresponding type.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = std::io::Result<SamAuxData>> {
        self.0.iter().map(|field| {
            let inv_aux_err_msg = || std::io::Error::other(format!("Invalid optional field {field}"));

            let (tag_text, rest) = field.split_once(':').ok_or_else(inv_aux_err_msg)?;
            let (type_text, value_text) = rest.split_once(':').ok_or_else(inv_aux_err_msg)?;

            let aux_tag = SamAuxData::parse_tag(tag_text)?;
            let aux_type = SamAuxData::parse_type(type_text)?;
            let sam_aux = SamAuxData::parse_value(aux_tag, aux_type, value_text)
                .with_context(format!("Failed to parse field '{field}'"))?;
            Ok(sam_aux)
        })
    }

    /// Retrieves and parses the value for an optional field with a given TAG in
    /// the SAM file format.
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
    /// [`get`]: SamAuxRaw::get
    pub fn get(&self, tag: &str) -> std::io::Result<Option<SamAuxData>> {
        for field in &self.0 {
            let inv_aux_err_msg = || std::io::Error::other(format!("Invalid optional field {field}"));
            let (aux_tag_text, rest) = field.split_once(':').ok_or_else(inv_aux_err_msg)?;

            if aux_tag_text == tag {
                let (type_text, value_text) = rest.split_once(':').ok_or_else(inv_aux_err_msg)?;

                let aux_tag = SamAuxData::parse_tag(aux_tag_text)?;
                let aux_type = SamAuxData::parse_type(type_text)?;

                let sam_aux = match SamAuxData::parse_value(aux_tag, aux_type, value_text) {
                    Ok(sam_aux) => sam_aux,
                    Err(e) => {
                        return Err(std::io::Error::other(format!(
                            "Failed to parse field '{field}' due to error: {e}"
                        )));
                    }
                };
                return Ok(Some(sam_aux));
            }
        }
        Ok(None)
    }

    /// Adds a field to the [`SamAuxRaw`] struct.
    ///
    /// ## Validity
    ///
    /// The tag name being pushed should not already be present in `self`.
    #[inline]
    pub fn push(&mut self, tag: &str, data: &SamAuxValue) {
        self.0.push(format!("{tag}:{data}"));
    }
}

impl FromIterator<String> for SamAuxRaw {
    #[inline]
    fn from_iter<T: IntoIterator<Item = String>>(iter: T) -> Self {
        SamAuxRaw(Vec::from_iter(iter))
    }
}

/// A parsed optional field in the SAM file format.
#[derive(Clone, Debug)]
pub struct SamAuxData {
    /// The tag of the optional SAM field.
    pub tag:   [u8; 2],
    /// The value of the optional SAM field.
    pub value: SamAuxValue,
}

impl SamAuxData {
    /// Parses the tag for the optional SAM field from a string slice.
    fn parse_tag(tag: &str) -> std::io::Result<[u8; 2]> {
        let bytes = tag.as_bytes();
        if bytes.len() != 2 || !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphanumeric() {
            return Err(std::io::Error::other(format!("Invalid SAM optional tag {tag}")));
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

    /// Parses a [`SamAuxData`] from a tag, type, and value (as a string).
    ///
    /// ## Errors
    ///
    /// `aux_type` must contain a single valid character (`A`, `i`, `f`, `Z`,
    /// `H`, or `B`). The `string_value` must successfully parse into the
    /// corresponding type.
    fn parse_value(aux_tag: [u8; 2], aux_type: char, string_value: &str) -> std::io::Result<SamAuxData> {
        match aux_type {
            'A' => {
                let mut chars = string_value.chars();
                let Some(c) = chars.next() else {
                    return Err(std::io::Error::other("'A' field has empty value"));
                };
                if chars.next().is_some() {
                    return Err(std::io::Error::other("'A' field must contain excatly one character"));
                }
                if !c.is_ascii() {
                    return Err(std::io::Error::other("'A' field must be ASCII"));
                }
                Ok(SamAuxData {
                    tag:   aux_tag,
                    value: SamAuxValue::Char(c as u8),
                })
            }
            'i' => {
                let parsed = string_value
                    .parse::<i64>()
                    .map_err(std::io::Error::other)
                    .with_context("Error parsing 'i' field")?;
                Ok(SamAuxData {
                    tag:   aux_tag,
                    value: SamAuxValue::Int(parsed),
                })
            }
            'f' => {
                let parsed = string_value.parse::<f32>().with_context("Error parsing 'f' field")?;
                Ok(SamAuxData {
                    tag:   aux_tag,
                    value: SamAuxValue::Float(parsed),
                })
            }
            'Z' => Ok(SamAuxData {
                tag:   aux_tag,
                value: SamAuxValue::String(String::from(string_value)),
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
                Ok(SamAuxData {
                    tag:   aux_tag,
                    value: SamAuxValue::Hex(string_value.to_ascii_uppercase()),
                })
            }
            'B' => Ok(SamAuxData {
                tag:   aux_tag,
                value: SamAuxValue::parse_aux_array(string_value).with_context("Failed to parse 'B' array")?,
            }),
            _ => Err(std::io::Error::other(format!(
                "Unsupported SAM optional field type {aux_type}"
            ))),
        }
    }

    /// Returns the stored character from the [`SamAuxData`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn char(self) -> Option<u8> {
        match self.value {
            SamAuxValue::Char(c) => Some(c),
            _ => None,
        }
    }

    /// Returns the stored integer from the [`SamAuxData`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn int(self) -> Option<i64> {
        match self.value {
            SamAuxValue::Int(i) => Some(i),
            _ => None,
        }
    }

    /// Returns the stored floating point number from the [`SamAuxData`], or
    /// [`None`] if a different variant is present.
    #[inline]
    #[must_use]
    pub fn float(self) -> Option<f32> {
        match self.value {
            SamAuxValue::Float(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored string from the [`SamAuxData`], or [`None`] if a
    /// different variant is present.
    #[inline]
    #[must_use]
    pub fn string(self) -> Option<String> {
        match self.value {
            SamAuxValue::String(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored Hex string from the [`SamAuxData`], or [`None`] if a
    /// different variant is present.
    ///
    /// For example, the six-character Hex string "1AE301" represents the byte
    /// array `[0x1a, 0xe3, 0x1]`.
    #[inline]
    #[must_use]
    pub fn hex(self) -> Option<String> {
        match self.value {
            SamAuxValue::Hex(f) => Some(f),
            _ => None,
        }
    }

    /// Returns the stored [`AuxArray`] from the [`SamAuxData`], or [`None`] if
    /// a different variant is present.
    #[inline]
    #[must_use]
    pub fn array(self) -> Option<AuxArray> {
        match self.value {
            SamAuxValue::Array(f) => Some(f),
            _ => None,
        }
    }
}

/// Value of the optional SAM field
#[derive(Clone, Debug)]
pub enum SamAuxValue {
    /// Printable character, type code `A`
    Char(u8),
    /// Signed integer, type code `i`
    Int(i64),
    /// Single-precision floating number, type code `f`
    Float(f32),
    /// Printable string, including space, type code `Z`
    String(String),
    /// Byte array in the Hex format, type code `H`
    ///
    /// For example, the six-character Hex string "1AE301" represents the byte
    /// array `[0x1a, 0xe3, 0x1]`.
    Hex(String),
    // Integer or numeric array, type code `B`
    Array(AuxArray),
}

impl SamAuxValue {
    /// Parses the array of optional SAM fields with type `B`.
    ///
    /// ## Errors
    ///
    /// The first letter in the array indicates the type of numbers in the
    /// following comma separated array. The letter can be one of `c`, `C`, `s`,
    /// `S`, `i`, `I`, or `f`.
    fn parse_aux_array(string_value: &str) -> std::io::Result<Self> {
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

                Ok(Self::Array(AuxArray::I8(values)))
            }
            "C" => {
                let values = pieces
                    .map(str::parse::<u8>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'C' subtype (`u8`)")?;
                Ok(Self::Array(AuxArray::U8(values)))
            }
            "s" => {
                let values = pieces
                    .map(str::parse::<i16>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 's' subtype (`i16`)")?;
                Ok(Self::Array(AuxArray::I16(values)))
            }
            "S" => {
                let values = pieces
                    .map(str::parse::<u16>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'S' subtype (`u16`)")?;
                Ok(Self::Array(AuxArray::U16(values)))
            }
            "i" => {
                let values = pieces
                    .map(str::parse::<i32>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'i' subtype (`i32`)")?;
                Ok(Self::Array(AuxArray::I32(values)))
            }
            "I" => {
                let values = pieces
                    .map(str::parse::<u32>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'I' subtype (`u32`)")?;
                Ok(Self::Array(AuxArray::U32(values)))
            }
            "f" => {
                let values = pieces
                    .map(str::parse::<f32>)
                    .process_results(|iter| iter.collect())
                    .with_context("Error parsing 'f' subtype (`f32`)")?;
                Ok(Self::Array(AuxArray::F32(values)))
            }
            _ => Err(std::io::Error::other(format!("Unsupported subtype {subtype}"))),
        }
    }
}

impl Display for SamAuxValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SamAuxValue::Char(val) => write!(f, "A:{val}", val = *val as char),
            SamAuxValue::Int(val) => write!(f, "i:{val}"),
            SamAuxValue::Float(val) => write!(f, "f:{val}"),
            SamAuxValue::String(val) => write!(f, "Z:{val}"),
            SamAuxValue::Hex(val) => write!(f, "H:{val}"),
            SamAuxValue::Array(val) => match val {
                AuxArray::I8(vec) => AuxArray::fmt_aux_array(f, 'c', vec),
                AuxArray::U8(vec) => AuxArray::fmt_aux_array(f, 'C', vec),
                AuxArray::I16(vec) => AuxArray::fmt_aux_array(f, 's', vec),
                AuxArray::U16(vec) => AuxArray::fmt_aux_array(f, 'S', vec),
                AuxArray::I32(vec) => AuxArray::fmt_aux_array(f, 'i', vec),
                AuxArray::U32(vec) => AuxArray::fmt_aux_array(f, 'I', vec),
                AuxArray::F32(vec) => AuxArray::fmt_aux_array(f, 'f', vec),
            },
        }
    }
}

/// Auxiliary data array for `B` field data.
#[derive(Debug, Clone)]
pub enum AuxArray {
    /// Array subtype code `c`
    I8(Vec<i8>),
    /// Array subtype code `C`
    U8(Vec<u8>),
    /// Array subtype code `s`
    I16(Vec<i16>),
    /// Array subtype code `S`
    U16(Vec<u16>),
    /// Array subtype code `i`
    I32(Vec<i32>),
    /// Array subtype code `I`
    U32(Vec<u32>),
    /// Array subtype code `f`
    F32(Vec<f32>),
}

impl AuxArray {
    fn fmt_aux_array<T: Display>(f: &mut Formatter<'_>, aux_type: char, vals: &[T]) -> std::fmt::Result {
        write!(f, "B:{aux_type}")?;

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
