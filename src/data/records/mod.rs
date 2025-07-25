use std::{fs::File, io::ErrorKind, path::Path};

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

trait RecordReader {
    const RECORD_NAME: &str;

    /// Opens a file, checking to ensure that it is non-empty, and providing
    /// context for error messages.
    #[inline]
    fn open_nonempty_file<P: AsRef<Path>>(filename: P) -> std::io::Result<File> {
        let path = filename.as_ref();
        let file = File::open(path).map_err(|err| Self::diagnose_failed_file(path, &err))?;
        Self::check_nonempty_file(path, &file)?;
        Ok(file)
    }

    /// Wraps an error message with context when a record file fails to open.
    #[inline]
    fn diagnose_failed_file(path: &Path, err: &std::io::Error) -> std::io::Error {
        let kind = err.kind();
        match kind {
            ErrorKind::NotFound => std::io::Error::new(
                kind,
                format!(
                    "Could not find {record_name} file {path} (see error: {err})",
                    record_name = Self::RECORD_NAME,
                    path = path.display()
                ),
            ),
            ErrorKind::PermissionDenied => std::io::Error::new(
                kind,
                format!(
                    "Insufficient permissions to open {record_name} file {path} (see error: {err})",
                    record_name = Self::RECORD_NAME,
                    path = path.display()
                ),
            ),
            _ => std::io::Error::new(
                kind,
                format!(
                    "Failed to open {record_name} file {path} (see error: {err})",
                    record_name = Self::RECORD_NAME,
                    path = path.display()
                ),
            ),
        }
    }

    /// Checks whether a file is non-empty, and if it is empty, return an error
    /// message with context.
    #[inline]
    fn check_nonempty_file(path: &Path, file: &File) -> std::io::Result<()> {
        if file.metadata()?.len() == 0 {
            Err(std::io::Error::new(
                ErrorKind::InvalidData,
                format!(
                    "{record_name} file was empty at path {path}",
                    record_name = Self::RECORD_NAME,
                    path = path.display()
                ),
            ))
        } else {
            Ok(())
        }
    }
}
