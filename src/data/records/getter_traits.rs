//! Getter traits for fields in record types, such as headers and sequences.

use crate::{
    alignment::{LocalProfiles, SharedProfiles},
    data::{
        amino_acids::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
        fasta::{Fasta, FastaView, FastaViewMut},
        fastq::{FastQ, FastQView, FastQViewMut},
        nucleotides::{Nucleotides, NucleotidesView, NucleotidesViewMut},
        views::{AssocViewMutType, AssocViewType},
    },
};

/// Getter trait for structures providing read access to a header/name.
pub trait HeaderReadable {
    /// Gets the header from the record.
    #[must_use]
    fn header(&self) -> &str;
}

impl HeaderReadable for FastQ {
    #[inline]
    fn header(&self) -> &str {
        &self.header
    }
}

impl HeaderReadable for FastQView<'_> {
    #[inline]
    fn header(&self) -> &str {
        self.header
    }
}

impl HeaderReadable for FastQViewMut<'_> {
    #[inline]
    fn header(&self) -> &str {
        self.header
    }
}

impl<S> HeaderReadable for Fasta<S> {
    fn header(&self) -> &str {
        &self.header
    }
}

impl<S: AssocViewType> HeaderReadable for FastaView<'_, S> {
    fn header(&self) -> &str {
        self.header
    }
}

impl<S: AssocViewMutType> HeaderReadable for FastaViewMut<'_, S> {
    fn header(&self) -> &str {
        self.header
    }
}

/// Getter trait for structures providing mutable access to a header/name.
pub trait HeaderMutable {
    /// Gets the header from the record.
    #[must_use]
    fn header_mut(&mut self) -> &mut String;
}

impl HeaderMutable for FastQ {
    #[inline]
    fn header_mut(&mut self) -> &mut String {
        &mut self.header
    }
}

impl HeaderMutable for FastQViewMut<'_> {
    #[inline]
    fn header_mut(&mut self) -> &mut String {
        self.header
    }
}

impl<S> HeaderMutable for Fasta<S> {
    fn header_mut(&mut self) -> &mut String {
        &mut self.header
    }
}

impl<S: AssocViewMutType> HeaderMutable for FastaViewMut<'_, S> {
    fn header_mut(&mut self) -> &mut String {
        self.header
    }
}

/// Getter trait for structures providing read access to a sequence.
///
/// The sequence can be either nucleotides or amino acids, and is returned as a
/// byte slice.
pub trait SequenceReadable {
    /// Get the sequence from the struct as a byte slice.
    #[must_use]
    fn sequence_bytes(&self) -> &[u8];
}

impl SequenceReadable for Nucleotides {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for NucleotidesView<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for NucleotidesViewMut<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for AminoAcids {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for AminoAcidsView<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for AminoAcidsViewMut<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for FastQ {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastQView<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastQViewMut<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl<S: AsRef<[u8]>> SequenceReadable for Fasta<S> {
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl<S> SequenceReadable for FastaView<'_, S>
where
    for<'a> S: AssocViewType<View<'a>: AsRef<[u8]>>,
{
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl<S> SequenceReadable for FastaViewMut<'_, S>
where
    for<'a> S: AssocViewMutType<ViewMut<'a>: AsRef<[u8]>>,
{
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> SequenceReadable for LocalProfiles<'_, M, N, O, S> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        &self.seq
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> SequenceReadable for SharedProfiles<'_, M, N, O, S> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        &self.seq
    }
}
