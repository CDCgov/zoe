use crate::{
    alignment::{LocalProfiles, SharedProfiles},
    data::fasta::{FastaAA, FastaNT, FastaNTAnnot, FastaSeq},
    prelude::{
        AminoAcids, AminoAcidsView, AminoAcidsViewMut, FastQ, FastQView, FastQViewMut, Nucleotides, NucleotidesView,
        NucleotidesViewMut,
    },
};

/// A module for reading and manipulating
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.
pub mod fasta;
/// A module for reading and manipulating
/// [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files.
pub mod fastq;
/// A module for reading and manipulating
/// [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files. Provides some
/// special-case functions used by [IRMA](https://wonder.cdc.gov/amd/flu/irma/).
pub mod sam;

/// Getter trait for structures providing read access to a header/name
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

impl HeaderReadable for FastaSeq {
    #[inline]
    fn header(&self) -> &str {
        &self.name
    }
}

impl HeaderReadable for FastaNT {
    #[inline]
    fn header(&self) -> &str {
        &self.name
    }
}

impl HeaderReadable for FastaAA {
    #[inline]
    fn header(&self) -> &str {
        &self.name
    }
}

impl HeaderReadable for FastaNTAnnot {
    #[inline]
    fn header(&self) -> &str {
        &self.name
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

impl SequenceReadable for FastaSeq {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastaAA {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastaNT {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastaNTAnnot {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> SequenceReadable for LocalProfiles<'_, M, N, O, S> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.seq
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> SequenceReadable for SharedProfiles<'_, M, N, O, S> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.seq
    }
}
