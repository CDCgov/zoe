use std::{error::Error, fs::File, io::ErrorKind, path::Path};

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

#[derive(Debug)]
struct RecordError {
    description: String,
    source:      Box<dyn Error + Send + Sync>,
}

impl std::fmt::Display for RecordError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.description)
    }
}

impl Error for RecordError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(self.source.as_ref())
    }
}

trait RecordReader {
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

    fn new_wrapped(description: &str, path: &Path, err: std::io::Error) -> std::io::Error {
        std::io::Error::other(RecordError {
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
