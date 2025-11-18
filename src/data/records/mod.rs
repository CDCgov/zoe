use std::{fs::File, io::ErrorKind, path::Path};

use crate::{
    data::{
        err::WithErrorContext,
        fasta::{FastaAA, FastaNT, FastaNTAnnot, FastaSeq},
    },
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

/// A trait providing helper functions when implementing a reader for a record
/// type (e.g., [`FastQReader`], [`FastaReader`], or [`SAMReader`]).
///
/// [`FastQReader`]: fastq::FastQReader
/// [`FastaReader`]: fasta::FastaReader
/// [`SAMReader`]: sam::SAMReader
trait RecordReader {
    /// The name of the record type, for use in error messages.
    const RECORD_NAME: &str;

    /// Opens a file, checking to ensure that it is non-empty, and providing
    /// context for error messages.
    #[inline]
    fn open_nonempty_file<P: AsRef<Path>>(filename: P) -> std::io::Result<File> {
        let path = filename.as_ref();

        let file = File::open(path).map_err(|err| Self::new_wrapped("file open error", path, err))?;
        let metadata = file
            .metadata()
            .map_err(|err| Self::new_wrapped("metadata error", path, err))?;
        if metadata.len() == 0 {
            return Err(Self::new_kind("file empty", path, ErrorKind::InvalidInput));
        }

        Ok(file)
    }

    /// Creates a new [`std::io::Error`] with a particular [`ErrorKind`]. The
    /// error message will contain a description, the name of the record
    /// [`RECORD_NAME`], and the path to the file for which the error occurred.
    ///
    /// [`RECORD_NAME`]: RecordReader::RECORD_NAME
    #[inline]
    #[must_use]
    fn new_kind(description: &str, path: &Path, kind: ErrorKind) -> std::io::Error {
        std::io::Error::new(
            kind,
            format!(
                "{desc} for {name}: '{path}'",
                desc = description,
                name = Self::RECORD_NAME,
                path = path.display()
            ),
        )
    }

    /// Given a [`std::io::Error`], wrap it in another error with additional
    /// context provided.
    ///
    /// Specifically, the original error `err` is first wrapped in a
    /// [`ErrorWithContext`], with the `source` field pointing to the original
    /// error. This is then put into an [`std::io::Error`], since this is the
    /// main error type used by *Zoe*'s readers.
    ///
    /// The error message will contain a description, the name of the record
    /// [`RECORD_NAME`], and the path to the file for which the error occurred.
    ///
    /// [`RECORD_NAME`]: RecordReader::RECORD_NAME
    /// [`ErrorWithContext`]: crate::data::err::ErrorWithContext
    fn new_wrapped(description: &str, path: &Path, err: std::io::Error) -> std::io::Error {
        std::io::Error::other(err.with_context(format!(
            "{desc} for {name}: '{path}'",
            desc = description,
            name = Self::RECORD_NAME,
            path = path.display()
        )))
    }
}

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

/// Getter trait for structures providing read access to a sequences.
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
