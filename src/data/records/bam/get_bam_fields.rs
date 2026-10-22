use crate::data::{
    cigar::ToCigletIterator,
    sam::{Flag, GetSamFields, ToOptFieldsIterator},
};

/// A trait for generic record types that are compatible with [`BamWriter`].
///
/// Specifically, this must implement [`GetSamFields`] with the following
/// guarantees:
///
/// - [`Qname`] and [`Rname`] implement `AsRef<str>`, such as [`String`] or
///   [`str`]
/// - [`GetSamFields::Flag`] is *Zoe*'s [`Flag`] type
/// - [`Pos`] is `usize` and [`Mapq`] is `u8`
/// - [`GetSamFields::Cigar`] implements [`ToCigletIterator`], such as
///   [`Cigar`], [`CigarView`], or [`AlignmentStates`]
/// - [`Seq`] implements `AsRef<[u8]>`, such as [`Nucleotides`],
///   [`NucleotidesView`], `Vec<u8>`, or `[u8]`
/// - [`Qual`] implements `AsRef<[u8]>`, such as [`QualityScores`],
///   [`QualityScoresView`], `Vec<u8>`, or `[u8]`
/// - [`OptFields`] implements [`ToOptFieldsIterator`], such as [`SamOptRaw`],
///   [`SamOptRawView`], or `[SamOptField]`
///
/// [`BamWriter`]: crate::data::records::bam::writer::BamWriter
/// [`Qname`]: GetSamFields::Qname
/// [`Rname`]: GetSamFields::Rname
/// [`Pos`]: GetSamFields::Pos
/// [`Mapq`]: GetSamFields::Mapq
/// [`Cigar`]: crate::data::types::cigar::Cigar
/// [`CigarView`]: crate::data::types::cigar::CigarView
/// [`AlignmentStates`]: crate::alignment::AlignmentStates
/// [`Seq`]: GetSamFields::Seq
/// [`Nucleotides`]: crate::data::types::nucleotides::Nucleotides
/// [`NucleotidesView`]: crate::data::types::nucleotides::NucleotidesView
/// [`Qual`]: GetSamFields::Qual
/// [`QualityScores`]: crate::data::types::phred::QualityScores
/// [`QualityScoresView`]: crate::data::types::phred::QualityScoresView
/// [`OptFields`]: GetSamFields::OptFields
/// [`SamOptRaw`]: crate::data::records::sam::SamOptRaw
/// [`SamOptRawView`]: crate::data::records::sam::SamOptRawView
pub trait EncodeBamFields:
    GetSamFields<
        Qname: AsRef<str>,
        Flag = Flag,
        Rname: AsRef<str>,
        Pos = usize,
        Mapq = u8,
        Cigar: ToCigletIterator,
        Seq: AsRef<[u8]>,
        Qual: AsRef<[u8]>,
        OptFields: ToOptFieldsIterator,
    > {
}

impl<T> EncodeBamFields for T
where
    T: GetSamFields<Flag = Flag, Pos = usize, Mapq = u8>,
    T::Qname: AsRef<str>,
    T::Rname: AsRef<str>,
    T::Cigar: ToCigletIterator,
    T::Seq: AsRef<[u8]>,
    T::Qual: AsRef<[u8]>,
    T::OptFields: ToOptFieldsIterator,
{
}
