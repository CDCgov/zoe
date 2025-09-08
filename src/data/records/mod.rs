use std::{error::Error, fs::File, io::ErrorKind, path::Path};

use crate::{
    data::fasta::{FastaAA, FastaNT, FastaNTAnnot, FastaSeq},
    prelude::{FastQ, FastQView, FastQViewMut},
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

/// A wrapper around an error with a new message, and the original error
/// accessible via [`Error::source`].
#[derive(Debug)]
struct RecordErrorContext {
    description: String,
    source:      Box<dyn Error + Send + Sync>,
}

impl std::fmt::Display for RecordErrorContext {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.description)
    }
}

impl Error for RecordErrorContext {
    #[inline]
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(self.source.as_ref())
    }
}

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
    /// [`RecordErrorContext`], with the `source` field pointing to the original error.
    /// This is then put into an [`std::io::Error`], since this is the main
    /// error type used by *Zoe*'s readers.
    ///
    /// The error message will contain a description, the name of the record
    /// [`RECORD_NAME`], and the path to the file for which the error occurred.
    ///
    /// [`RECORD_NAME`]: RecordReader::RECORD_NAME
    fn new_wrapped(description: &str, path: &Path, err: std::io::Error) -> std::io::Error {
        std::io::Error::other(RecordErrorContext {
            description: format!(
                "{desc} for {name}: '{path}'",
                desc = description,
                name = Self::RECORD_NAME,
                path = path.display()
            ),
            source:      Box::new(err),
        })
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
