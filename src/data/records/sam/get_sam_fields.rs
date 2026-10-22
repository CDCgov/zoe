//! Code for [`GetSamFields`] and [`GetSamFieldsMut`], getter trait for unifying
//! SAM-like record structs.

use crate::{
    alignment::AlignmentStates,
    data::{
        cigar::{Cigar, CigarView},
        err::{ErrorWithContext, ResultWithErrorContext},
        nucleotides::{Nucleotides, NucleotidesView},
        phred::{QualityScores, QualityScoresView},
        sam::{Flag, SamData, SamDataView, SamOptField, SamOptRaw, SamOptRawView, is_missing_sam_field},
    },
};
use std::{io::ErrorKind, str::FromStr};

/// A trait enabling read access to SAM-like data in generic contexts. This can
/// then be extended by [`DisplaySamFields`] for the purposes of writing SAM
/// data or [`EncodeBamFields`] for the purpose of writing BAM data.
///
/// Implementing this trait involves specifying the type used to store each
/// field as well as providing getters for those fields. The generics allow any
/// types to be used, thereby enabling owned and view types, parsed and unparsed
/// types, and *Zoe* and foreign types to all be compatible. For example,
/// implementors could use [`Cigar`], [`CigarView`], [`AlignmentStates`], or
/// [`String`] for the CIGAR field.
///
/// To support writing to a SAM file, ensure all associated types implement
/// [`Display`] (see [`DisplaySamFields`]). To support writing to a BAM file,
/// ensure the requirements listed in [`EncodeBamFields`] are satisfied.
///
/// Furthermore, holding any of the fields is optional, in which case the getter
/// should return `None`. For non-stored fields, the associated type should be
/// set to some reasonable default (e.g., [`SamOptRaw`] for [`OptFields`] or
/// [`Nucleotides`] for [`Seq`]) to ensure compatibility with
/// [`EncodeBamFields`]. In this case, the display implementation automatically
/// uses the apporiate defaults.
///
/// The methods of this trait should typically not be used directly by users,
/// aside from implementing them. Prefer direct field access instead.
///
/// [`Display`]: std::fmt::Display
/// [`DisplaySamFields`]: crate::data::records::sam::DisplaySamFields
/// [`EncodeBamFields`]: crate::data::records::bam::EncodeBamFields
/// [`Seq`]: GetSamFields::Seq
/// [`OptFields`]: GetSamFields::OptFields
/// [`AlignmentStates`]: crate::alignment::AlignmentStates
pub trait GetSamFields {
    /// The type stored for QNAME, such as [`String`] or [`str`].
    type Qname: ?Sized;

    /// The type stored for FLAG, typically [`Flag`].
    type Flag: ?Sized;

    /// The type stored for RNAME, such as [`String`] or [`str`].
    type Rname: ?Sized;

    /// The type stored for POS, typically `usize`.
    type Pos: ?Sized;

    /// The type stored for MAPQ, typically `u8`.
    type Mapq: ?Sized;

    /// The type stored for CIGAR, such as [`Cigar`], [`CigarView`], or
    /// [`AlignmentStates`].
    ///
    /// [`AlignmentStates`]: crate::alignment::AlignmentStates
    type Cigar: ?Sized;

    /// The type stored for SEQ, such as [`Nucleotides`] or [`NucleotidesView`].
    type Seq: ?Sized;

    /// The type stored for QUAL, such as [`QualityScores`] or
    /// [`QualityScoresView`].
    type Qual: ?Sized;

    /// The type storing the optional fields, such as [`SamOptRaw`],
    /// [`SamOptRawView`], or a custom struct that holds parsed [`SamOptField`].
    ///
    /// [`SamOptField`]: crate::data::records::sam::SamOptField
    type OptFields: ?Sized;

    /// Returns the QNAME field if it is contained in the struct.
    fn qname(&self) -> Option<&Self::Qname>;

    /// Returns the FLAG field if it is contained in the struct.
    fn flag(&self) -> Option<&Self::Flag>;

    /// Returns the RNAME field if it is contained in the struct.
    fn rname(&self) -> Option<&Self::Rname>;

    /// Returns the POS field if it is contained in the struct.
    fn pos(&self) -> Option<&Self::Pos>;

    /// Returns the MAPQ field if it is contained in the struct.
    fn mapq(&self) -> Option<&Self::Mapq>;

    /// Returns the CIGAR field if it is present in the struct.
    fn cigar(&self) -> Option<&Self::Cigar>;

    /// Returns the SEQ field if it is present in the struct.
    fn seq(&self) -> Option<&Self::Seq>;

    /// Returns the QUAL field if it is present in the struct.
    fn qual(&self) -> Option<&Self::Qual>;

    /// Returns any optional SAM fields stored in the struct.
    fn opt_fields(&self) -> Option<&Self::OptFields>;
}

impl<'a, T: GetSamFields> GetSamFields for &'a T
where
    T: 'a,
{
    type Qname = T::Qname;
    type Flag = T::Flag;
    type Rname = T::Rname;
    type Pos = T::Pos;
    type Mapq = T::Mapq;
    type Cigar = T::Cigar;
    type Seq = T::Seq;
    type Qual = T::Qual;
    type OptFields = T::OptFields;

    #[inline]
    fn qname(&self) -> Option<&Self::Qname> {
        (*self).qname()
    }

    #[inline]
    fn flag(&self) -> Option<&Self::Flag> {
        (*self).flag()
    }

    #[inline]
    fn rname(&self) -> Option<&Self::Rname> {
        (*self).rname()
    }

    #[inline]
    fn pos(&self) -> Option<&Self::Pos> {
        (*self).pos()
    }

    #[inline]
    fn mapq(&self) -> Option<&Self::Mapq> {
        (*self).mapq()
    }

    #[inline]
    fn cigar(&self) -> Option<&Self::Cigar> {
        (*self).cigar()
    }

    #[inline]
    fn seq(&self) -> Option<&Self::Seq> {
        (*self).seq()
    }

    #[inline]
    fn qual(&self) -> Option<&Self::Qual> {
        (*self).qual()
    }

    #[inline]
    fn opt_fields(&self) -> Option<&Self::OptFields> {
        (*self).opt_fields()
    }
}

impl GetSamFields for SamData {
    type Qname = String;
    type Flag = Flag;
    type Rname = String;
    type Pos = usize;
    type Mapq = u8;
    type Cigar = Cigar;
    type Seq = Nucleotides;
    type Qual = QualityScores;
    type OptFields = SamOptRaw;

    #[inline]
    fn qname(&self) -> Option<&String> {
        Some(&self.qname)
    }

    #[inline]
    fn flag(&self) -> Option<&Flag> {
        Some(&self.flag)
    }

    #[inline]
    fn rname(&self) -> Option<&String> {
        Some(&self.rname)
    }

    #[inline]
    fn pos(&self) -> Option<&usize> {
        Some(&self.pos)
    }

    #[inline]
    fn mapq(&self) -> Option<&u8> {
        Some(&self.mapq)
    }

    #[inline]
    fn cigar(&self) -> Option<&Cigar> {
        Some(&self.cigar)
    }

    #[inline]
    fn seq(&self) -> Option<&Nucleotides> {
        Some(&self.seq)
    }

    #[inline]
    fn qual(&self) -> Option<&QualityScores> {
        Some(&self.qual)
    }

    #[inline]
    fn opt_fields(&self) -> Option<&SamOptRaw> {
        Some(&self.opt_fields)
    }
}

impl<'a> GetSamFields for SamDataView<'a> {
    type Qname = str;
    type Flag = Flag;
    type Rname = str;
    type Pos = usize;
    type Mapq = u8;
    type Cigar = CigarView<'a>;
    type Seq = NucleotidesView<'a>;
    type Qual = QualityScoresView<'a>;
    type OptFields = SamOptRawView<'a>;

    #[inline]
    fn qname(&self) -> Option<&str> {
        Some(self.qname)
    }

    #[inline]
    fn flag(&self) -> Option<&Flag> {
        Some(&self.flag)
    }

    #[inline]
    fn rname(&self) -> Option<&str> {
        Some(self.rname)
    }

    #[inline]
    fn pos(&self) -> Option<&usize> {
        Some(&self.pos)
    }

    #[inline]
    fn mapq(&self) -> Option<&u8> {
        Some(&self.mapq)
    }

    #[inline]
    fn cigar(&self) -> Option<&CigarView<'a>> {
        Some(&self.cigar)
    }

    #[inline]
    fn seq(&self) -> Option<&NucleotidesView<'a>> {
        Some(&self.seq)
    }

    #[inline]
    fn qual(&self) -> Option<&QualityScoresView<'a>> {
        Some(&self.qual)
    }

    #[inline]
    fn opt_fields(&self) -> Option<&SamOptRawView<'a>> {
        Some(&self.opt_fields)
    }
}

/// A trait enabling write access to SAM-like data in generic contexts. This can
/// then be extended by [`ParseSamFields`] for the purposes of parsing SAM data
/// from a SAM file.
///
/// To implement this, first implement [`GetSamFields`].
///
/// The methods of this trait should typically not be used directly by users,
/// aside from implementing them. Prefer direct field access instead.
pub trait GetSamFieldsMut: GetSamFields {
    /// Returns the QNAME field if it is contained in the struct.
    fn qname_mut(&mut self) -> Option<&mut Self::Qname>;

    /// Returns the FLAG field if it is contained in the struct.
    fn flag_mut(&mut self) -> Option<&mut Self::Flag>;

    /// Returns the RNAME field if it is contained in the struct.
    fn rname_mut(&mut self) -> Option<&mut Self::Rname>;

    /// Returns the POS field if it is contained in the struct.
    fn pos_mut(&mut self) -> Option<&mut Self::Pos>;

    /// Returns the MAPQ field if it is contained in the struct.
    fn mapq_mut(&mut self) -> Option<&mut Self::Mapq>;

    /// Returns the CIGAR field if it is present in the struct.
    fn cigar_mut(&mut self) -> Option<&mut Self::Cigar>;

    /// Returns the SEQ field if it is present in the struct.
    fn seq_mut(&mut self) -> Option<&mut Self::Seq>;

    /// Returns the QUAL field if it is present in the struct.
    fn qual_mut(&mut self) -> Option<&mut Self::Qual>;

    /// Returns any optional SAM fields stored in the struct.
    fn opt_fields_mut(&mut self) -> Option<&mut Self::OptFields>;
}

impl GetSamFieldsMut for SamData {
    #[inline]
    fn qname_mut(&mut self) -> Option<&mut Self::Qname> {
        Some(&mut self.qname)
    }

    #[inline]
    fn flag_mut(&mut self) -> Option<&mut Self::Flag> {
        Some(&mut self.flag)
    }

    #[inline]
    fn rname_mut(&mut self) -> Option<&mut Self::Rname> {
        Some(&mut self.rname)
    }

    #[inline]
    fn pos_mut(&mut self) -> Option<&mut Self::Pos> {
        Some(&mut self.pos)
    }

    #[inline]
    fn mapq_mut(&mut self) -> Option<&mut Self::Mapq> {
        Some(&mut self.mapq)
    }

    #[inline]
    fn cigar_mut(&mut self) -> Option<&mut Self::Cigar> {
        Some(&mut self.cigar)
    }

    #[inline]
    fn seq_mut(&mut self) -> Option<&mut Self::Seq> {
        Some(&mut self.seq)
    }

    #[inline]
    fn qual_mut(&mut self) -> Option<&mut Self::Qual> {
        Some(&mut self.qual)
    }

    #[inline]
    fn opt_fields_mut(&mut self) -> Option<&mut Self::OptFields> {
        Some(&mut self.opt_fields)
    }
}

/// An extension of [`GetSamFieldsMut`] enabling parsing a line of a SAM file
/// into a struct.
///
/// For this to be implemented, [`Default`] must be implemented, and the
/// associated types of [`GetSamFields`] must be able to be parsed from the
/// corresponding unparsed wrapper types.
///
/// For example, if a custom struct for [`OptFields`] is used, implement
/// [`From`], [`Into`], or [`TryFrom`] to provide a conversion from
/// [`UnparsedOptFields`] to the struct. For fallible conversion, the error type
/// can be [`std::io::Error`] or *Zoe*'s [`ErrorWithContext`].
///
/// [`OptFields`]: crate::data::records::sam::GetSamFields::OptFields
pub trait ParseSamFields:
    Default
    + for<'a> GetSamFieldsMut<
        Qname: TryFrom<UnparsedQname<'a>, Error: Into<ErrorWithContext>>,
        Flag: TryFrom<UnparsedFlag<'a>, Error: Into<ErrorWithContext>>,
        Rname: TryFrom<UnparsedRname<'a>, Error: Into<ErrorWithContext>>,
        Pos: TryFrom<UnparsedPos<'a>, Error: Into<ErrorWithContext>>,
        Mapq: TryFrom<UnparsedMapq<'a>, Error: Into<ErrorWithContext>>,
        Cigar: TryFrom<UnparsedCigar<'a>, Error: Into<ErrorWithContext>>,
        Seq: TryFrom<UnparsedSeq<'a>, Error: Into<ErrorWithContext>>,
        Qual: TryFrom<UnparsedQual<'a>, Error: Into<ErrorWithContext>>,
        OptFields: TryFrom<UnparsedOptFields<'a>, Error: Into<ErrorWithContext>>,
    > {
    /// Parses a line of a SAM file (with tab-separated fields) into the struct.
    ///
    /// ## Errors
    ///
    /// At least 11 TSV fields are required. Any errors parsing the fields are
    /// propagated.
    fn parse(line: &str) -> std::io::Result<Self> {
        // 11 required fields and one slot for optional fields
        let parts = line.splitn(12, '\t').collect::<Vec<_>>();

        if parts.len() < 11 {
            return Err(std::io::Error::new(
                ErrorKind::InvalidData,
                format!("At least 11 fields were expected, but found {}", parts.len()),
            ));
        }

        let mut out = Self::default();

        if let Some(qname) = out.qname_mut() {
            *qname = UnparsedQname(parts[0]).try_into().map_err(Into::into)?;
        }

        if let Some(flag) = out.flag_mut() {
            *flag = UnparsedFlag(parts[1]).try_into().map_err(Into::into)?;
        }

        if let Some(rname) = out.rname_mut() {
            *rname = UnparsedRname(parts[2]).try_into().map_err(Into::into)?;
        }

        if let Some(pos) = out.pos_mut() {
            *pos = UnparsedPos(parts[3]).try_into().map_err(Into::into)?;
        }

        if let Some(mapq) = out.mapq_mut() {
            *mapq = UnparsedMapq(parts[4]).try_into().map_err(Into::into)?;
        }

        if let Some(cigar) = out.cigar_mut() {
            *cigar = UnparsedCigar(parts[5]).try_into().map_err(Into::into)?;
        }

        if let Some(seq) = out.seq_mut() {
            *seq = UnparsedSeq(parts[9]).try_into().map_err(Into::into)?;
        }

        if let Some(qual) = out.qual_mut() {
            *qual = UnparsedQual(parts[10]).try_into().map_err(Into::into)?;
        }

        if let Some(opt_fields) = out.opt_fields_mut()
            && let Some(field) = parts.get(11)
        {
            *opt_fields = UnparsedOptFields(field).try_into().map_err(Into::into)?;
        }

        Ok(out)
    }
}

impl<T> ParseSamFields for T where
    T: Default
        + for<'a> GetSamFieldsMut<
            Qname: TryFrom<UnparsedQname<'a>, Error: Into<ErrorWithContext>>,
            Flag: TryFrom<UnparsedFlag<'a>, Error: Into<ErrorWithContext>>,
            Rname: TryFrom<UnparsedRname<'a>, Error: Into<ErrorWithContext>>,
            Pos: TryFrom<UnparsedPos<'a>, Error: Into<ErrorWithContext>>,
            Mapq: TryFrom<UnparsedMapq<'a>, Error: Into<ErrorWithContext>>,
            Cigar: TryFrom<UnparsedCigar<'a>, Error: Into<ErrorWithContext>>,
            Seq: TryFrom<UnparsedSeq<'a>, Error: Into<ErrorWithContext>>,
            Qual: TryFrom<UnparsedQual<'a>, Error: Into<ErrorWithContext>>,
            OptFields: TryFrom<UnparsedOptFields<'a>, Error: Into<ErrorWithContext>>,
        >
{
}

/// An unparsed QNAME field, used as a helper for [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedQname<'a>(&'a str);

/// An unparsed FLAG field, used as a helper for [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedFlag<'a>(&'a str);

/// An unparsed RNAME field, used as a helper for [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedRname<'a>(&'a str);

/// An unparsed POS field, used as a helper for [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedPos<'a>(&'a str);

/// An unparsed MAPQ field, used as a helper for [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedMapq<'a>(&'a str);

/// An unparsed CIGAR field, used as a helper for [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedCigar<'a>(&'a str);

/// An unparsed SEQ field, used as a helper for [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedSeq<'a>(&'a str);

/// An unparsed QUAL field, used as a helper for [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedQual<'a>(&'a str);

/// The unparsed optional fields for a SAM record, used as a helper for
/// [`ParseSamFields`].
#[derive(Copy, Clone)]
pub struct UnparsedOptFields<'a>(&'a str);

impl<'a> UnparsedQname<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl<'a> UnparsedFlag<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl<'a> UnparsedRname<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl<'a> UnparsedPos<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl<'a> UnparsedMapq<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl<'a> UnparsedCigar<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl<'a> UnparsedSeq<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl<'a> UnparsedQual<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl<'a> UnparsedOptFields<'a> {
    /// Returns the inner [`str`], for use when implementing a conversion from
    /// this to a custom type.
    #[inline]
    #[must_use]
    pub fn as_str(self) -> &'a str {
        self.0
    }
}

impl From<UnparsedQname<'_>> for String {
    #[inline]
    fn from(value: UnparsedQname) -> Self {
        if is_missing_sam_field(value.0) {
            String::new()
        } else {
            value.0.to_string()
        }
    }
}

impl From<UnparsedFlag<'_>> for String {
    #[inline]
    fn from(value: UnparsedFlag) -> Self {
        value.0.to_string()
    }
}

impl TryFrom<UnparsedFlag<'_>> for u16 {
    type Error = std::io::Error;

    #[inline]
    fn try_from(value: UnparsedFlag) -> std::io::Result<Self> {
        value
            .0
            .parse()
            .with_context(format!("Failed to parse FLAG field: '{}'", value.0))
            .map_err(Into::into)
    }
}

impl TryFrom<UnparsedFlag<'_>> for Flag {
    type Error = std::io::Error;

    #[inline]
    fn try_from(value: UnparsedFlag) -> std::io::Result<Self> {
        value
            .0
            .parse()
            .with_context(format!("Failed to parse FLAG field: '{}'", value.0))
            .map_err(Into::into)
    }
}

impl From<UnparsedRname<'_>> for String {
    #[inline]
    fn from(value: UnparsedRname) -> Self {
        if is_missing_sam_field(value.0) {
            String::new()
        } else {
            value.0.to_string()
        }
    }
}

impl From<UnparsedPos<'_>> for String {
    #[inline]
    fn from(value: UnparsedPos) -> Self {
        value.0.to_string()
    }
}

impl TryFrom<UnparsedPos<'_>> for usize {
    type Error = std::io::Error;

    #[inline]
    fn try_from(value: UnparsedPos<'_>) -> std::io::Result<Self> {
        value
            .0
            .parse()
            .with_context(format!("Failed to parse POS field: '{}'", value.0))
            .map_err(Into::into)
    }
}

impl From<UnparsedMapq<'_>> for String {
    #[inline]
    fn from(value: UnparsedMapq) -> Self {
        value.0.to_string()
    }
}

impl TryFrom<UnparsedMapq<'_>> for u8 {
    type Error = std::io::Error;

    #[inline]
    fn try_from(value: UnparsedMapq) -> std::io::Result<Self> {
        value
            .0
            .parse()
            .with_context(format!("Failed to parse MAPQ field: '{}'", value.0))
            .map_err(Into::into)
    }
}

impl From<UnparsedCigar<'_>> for String {
    #[inline]
    fn from(value: UnparsedCigar) -> Self {
        let value = value.0.trim_ascii();
        if is_missing_sam_field(value) {
            String::new()
        } else {
            value.to_string()
        }
    }
}

impl From<UnparsedCigar<'_>> for Cigar {
    #[inline]
    fn from(value: UnparsedCigar<'_>) -> Self {
        let value = value.0.as_bytes().trim_ascii();
        // Handles checking for *
        Cigar::from_slice_unchecked(value)
    }
}

impl TryFrom<UnparsedCigar<'_>> for AlignmentStates {
    type Error = std::io::Error;

    #[inline]
    fn try_from(value: UnparsedCigar<'_>) -> std::io::Result<Self> {
        let value = value.0.trim_ascii();
        // Handles checking for *
        value
            .try_into()
            .with_context(format!("Failed to parse CIGAR field: '{value}'"))
            .map_err(Into::into)
    }
}

impl From<UnparsedSeq<'_>> for String {
    #[inline]
    fn from(value: UnparsedSeq) -> Self {
        if is_missing_sam_field(value.0) {
            String::new()
        } else {
            value.0.to_string()
        }
    }
}

impl From<UnparsedSeq<'_>> for Nucleotides {
    #[inline]
    fn from(value: UnparsedSeq<'_>) -> Self {
        if is_missing_sam_field(value.0) {
            Nucleotides::new()
        } else {
            value.0.into()
        }
    }
}

impl From<UnparsedQual<'_>> for String {
    #[inline]
    fn from(value: UnparsedQual) -> Self {
        if is_missing_sam_field(value.0) {
            String::new()
        } else {
            value.0.to_string()
        }
    }
}

impl TryFrom<UnparsedQual<'_>> for QualityScores {
    type Error = std::io::Error;

    #[inline]
    fn try_from(value: UnparsedQual<'_>) -> Result<Self, Self::Error> {
        if is_missing_sam_field(value.0) {
            Ok(QualityScores::new())
        } else {
            value
                .0
                .try_into()
                .with_context(format!("Failed to parse QUAL field: '{}'", value.0))
                .map_err(Into::into)
        }
    }
}

impl From<UnparsedOptFields<'_>> for String {
    #[inline]
    fn from(value: UnparsedOptFields) -> Self {
        value.0.to_string()
    }
}

impl From<UnparsedOptFields<'_>> for SamOptRaw {
    #[inline]
    fn from(value: UnparsedOptFields<'_>) -> Self {
        SamOptRaw(value.0.to_string())
    }
}

impl TryFrom<UnparsedOptFields<'_>> for Vec<SamOptField> {
    type Error = std::io::Error;

    #[inline]
    fn try_from(value: UnparsedOptFields<'_>) -> std::io::Result<Self> {
        value
            .0
            .split('\t')
            .map(FromStr::from_str)
            .collect::<Result<Vec<_>, _>>()
            .with_context(format!("Failed to parse optional fields: '{}'", value.0))
            .map_err(Into::into)
    }
}
