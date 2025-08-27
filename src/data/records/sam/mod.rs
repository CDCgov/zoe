use crate::{
    alignment::Alignment,
    data::types::cigar::{Cigar, CigarView, CigarViewMut},
    math::AnyInt,
    prelude::*,
};
use std::{fmt::Display, hash::Hash};

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
    /// Optional tags which can be lazily parsed/accessed
    pub tags:  SamTags,
}

impl PartialEq for SamData {
    /// Tests for `self` and `other` values to be equal, and is used by `==`.
    /// Note that this implementation ignores the `tags` field, which contains
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
    /// ignores the `tags` field, which contains optional SAM values.
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
            tags: SamTags::new(),
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
    pub fn from_alignment<T: AnyInt + Into<i64>>(
        alignment: &Alignment<T>, qname: String, flag: u16, rname: String, mapq: u8, seq: Nucleotides, qual: QualityScores,
    ) -> Self {
        // Both SAM and Alignment exclude clipped bases when reporting
        // positions, so we just need to adjust to 1-based
        let pos = alignment.ref_range.start + 1;
        let cigar = alignment.states.to_cigar_unchecked();
        let tags = SamTags::new_with_score(alignment.score);
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
            tags,
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

/// A wrapper around an array for holding optional SAM tags.
#[derive(Clone, Debug, Default)]
pub struct SamTags(Vec<String>);

impl SamTags {
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        SamTags(Vec::new())
    }

    #[inline]
    pub fn new_with_score<T: AnyInt + Into<i64>>(score: T) -> Self {
        let mut inner = Vec::with_capacity(1);
        inner.push(format!("AS:i:{score}"));
        SamTags(inner)
    }

    /// Retrieves and parses the value for an optional tag in the SAM file
    /// format.
    ///
    /// ## Limitations
    ///
    /// This struct parses tags lazily and will repeat the computations each
    /// time [`get`] is called. The function runs in $O(n)$ time where $n$ is
    /// the number of tags. Validation of the tag format is only performed where
    /// necessary.
    ///
    /// ## Errors
    ///
    /// The tag in the SAM record must be of the form `TAG:TYPE:VALUE`, with no
    /// other colons. `TYPE` must be supported by *Zoe* (either `A`, `i`, or
    /// `f`). `VALUE` must successfully parse into the corresponding type.
    ///
    /// [`get`]: SamTags::get
    pub fn get(&self, name: &str) -> std::io::Result<Option<SamTagValue>> {
        for tag in &self.0 {
            let Some((tag_name, rest)) = tag.split_once(':') else {
                return Err(std::io::Error::other(format!("Tag is missing first colon: {tag}")));
            };
            if tag_name == name {
                let Some((tag_type, tag_value)) = rest.split_once(':') else {
                    return Err(std::io::Error::other(format!("Tag is missing second colon: {tag}")));
                };
                let tag_value = match SamTagValue::parse(tag_type, tag_value) {
                    Ok(value) => value,
                    Err(e) => {
                        return Err(std::io::Error::other(format!(
                            "Failed to parse tag '{tag}' due to error: {e}"
                        )));
                    }
                };
                return Ok(Some(tag_value));
            }
        }
        Ok(None)
    }

    // /// Adds a tag without checking if it already exists.
    // ///
    // /// ## Validity
    // ///
    // /// In order to conform with the SAM format, the name should consist of two
    // /// uppercase letters.
    // pub fn push(&mut self, name: &str, value: SamTagValue) {
    //     let
    // }
}

impl FromIterator<String> for SamTags {
    #[inline]
    fn from_iter<T: IntoIterator<Item = String>>(iter: T) -> Self {
        SamTags(Vec::from_iter(iter))
    }
}

/// A parsed optional tag in the SAM file format
#[non_exhaustive]
pub enum SamTagValue {
    /// A printable character in the range `! ..= ~`.
    Char(u8),
    /// A signed integer
    Int(i64),
    /// A single-precision floating point number
    Float(f32),
}

impl SamTagValue {
    /// Parses a [`SamTagValue`] from a string for the type and a string for the
    /// value.
    ///
    /// ## Errors
    ///
    /// `tag_type` must contain a single valid character supported by *Zoe*
    /// (currently, `A`, `i`, or `f`). The `value` string must successfully
    /// parse into that type.
    fn parse(tag_type: &str, value: &str) -> std::io::Result<SamTagValue> {
        let Ok([tag_type]) = <[u8; 1]>::try_from(tag_type.as_bytes()) else {
            return Err(std::io::Error::other(format!("Tag type of {tag_type} is not valid")));
        };
        let out = match tag_type {
            b'A' => Self::Char(
                value
                    .parse()
                    .map_err(|_| std::io::Error::other(format!("Failed to parse character value {value}")))?,
            ),
            b'i' => Self::Int(
                value
                    .parse()
                    .map_err(|_| std::io::Error::other(format!("Failed to parse integer value {value}")))?,
            ),
            b'f' => Self::Float(
                value
                    .parse()
                    .map_err(|_| std::io::Error::other(format!("Failed to parse floating point value {value}")))?,
            ),
            _ => {
                return Err(std::io::Error::other(format!(
                    "Tag type {tag_type} is not supported in Zoe",
                    tag_type = tag_type as char
                )));
            }
        };
        Ok(out)
    }
}

impl Display for SamTagValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SamTagValue::Char(val) => write!(f, "A:{val}", val = *val as char),
            SamTagValue::Int(val) => write!(f, "i:{val}"),
            SamTagValue::Float(val) => write!(f, "f:{val}"),
        }
    }
}
