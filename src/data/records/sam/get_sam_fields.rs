//! Code for [`GetSamFields`], a getter trait for unifying SAM-like record
//! structs.

use crate::data::{
    cigar::{Cigar, CigarView},
    nucleotides::{Nucleotides, NucleotidesView},
    phred::{QualityScores, QualityScoresView},
    sam::{Flag, SamData, SamDataView, SamOptRaw, SamOptRawView},
};

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
