use crate::{
    data::{
        err::ResultWithErrorContext,
        records::reader_utils::{BufReadExtension, FusedBufReader},
    },
    prelude::*,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, Error as IOError, ErrorKind, Read},
    iter::FusedIterator,
    path::Path,
};

/// A buffered reader for reading
/// [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format).
///
/// This does not support multiline FASTQ files. In other words, each sequence
/// must be on a single line, and the quality scores must be on a single line.
#[derive(Debug)]
pub struct FastQReader<R: Read> {
    reader: FusedBufReader<R>,
}

impl<R: Read> FastQReader<R> {
    /// Creates an iterator over FASTQ data, wrapping the input in a buffered
    /// reader.
    ///
    /// Unlike [`from_readable`], this does not read any data initially. It also
    /// allows for empty input, in which case the resulting iterator is empty.
    ///
    /// The iterator returns `None` once `inner.read()` returns 0 (signalling
    /// EOF). After this, the iterator will continue to return `None` regardless
    /// of whether `inner.read()` could return new bytes.
    ///
    /// [`from_readable`]: FastQReader::from_readable
    pub fn new(inner: R) -> Self {
        FastQReader {
            reader: FusedBufReader::new(inner),
        }
    }

    /// Creates an iterator over FASTQ data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// The iterator returns `None` once `inner.read()` returns 0 (signalling
    /// EOF). After this, the iterator will continue to return `None` regardless
    /// of whether `inner.read()` could return new bytes.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_readable(read: R) -> std::io::Result<Self> {
        FastQReader::from_bufreader(BufReader::new(read))
    }

    /// Creates an iterator over FASTQ data from a [`BufReader`].
    ///
    /// The iterator returns `None` once `reader.fill_buf()` returns an empty
    /// slice (signalling EOF). After this, the iterator will continue to return
    /// `None` regardless of whether `reader.fill_buf()` could return new bytes.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_bufreader(mut reader: BufReader<R>) -> std::io::Result<Self> {
        if reader.retrying_fill_buf()?.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "No FASTQ data was found!"));
        }

        Ok(FastQReader {
            reader: FusedBufReader::from_bufreader(reader),
        })
    }
}

impl FastQReader<std::fs::File> {
    /// Creates an iterator over the FASTQ data contained in a path, using a
    /// buffered reader.
    ///
    /// The iterator returns `None` once EOF is reached. After this, the
    /// iterator will continue to return `None` regardless of whether additional
    /// bytes are added to the file.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the path does not exist, if there are insufficient
    /// permissions to read from it, or if it contains no data. The path is
    /// included in the error message.
    pub fn from_path<P>(path: P) -> Result<FastQReader<File>, std::io::Error>
    where
        P: AsRef<Path>, {
        let path = path.as_ref();
        let file = File::open(path).with_path_context("Failed to open path", path)?;
        Ok(Self::from_readable(file).with_path_context("Failed to read data at path", path)?)
    }
}

impl<R: Read> Iterator for FastQReader<R> {
    type Item = std::io::Result<FastQ>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_helper().transpose()
    }
}

impl<R: Read> FusedIterator for FastQReader<R> {}

impl<R: Read> FastQReader<R> {
    /// A helper function for [`FastQReader::next`] that returns a [`Result`]
    /// rather than an [`Option`], to allow for a more readable implementation.
    #[inline]
    fn next_helper(&mut self) -> std::io::Result<Option<FastQ>> {
        // Consume "@" at start of header, or skip line and abort
        let buf = self.reader.retrying_fill_buf()?;

        let Some(marker) = buf.first().copied() else {
            // EOF reached
            return Ok(None);
        };

        if marker != b'@' {
            // Skip the entire problematic line to avoid repeated errors (or
            // infinite loops, in the case of an improperly used reader)
            self.reader.skip_until(b'\n')?;

            return Err(IOError::new(
                ErrorKind::InvalidData,
                "Missing '@' symbol at header line beginning! Ensure that the FASTQ file is not multi-line.",
            ));
        }

        self.reader.consume(1);

        // Read header directly into Vec
        let header = self.reader.read_line_owned(0)?;

        if header.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "Missing FASTQ header!"));
        }

        let header = String::from_utf8(header).map_err(|err| IOError::new(ErrorKind::InvalidData, err))?;

        // Read sequence directly into Vec
        let sequence = self.reader.read_line_owned(0)?;

        if sequence.is_empty() {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                format!("Missing FASTQ sequence! See header: {header}"),
            ));
        }

        let sequence = Nucleotides(sequence);

        // Consume "+" line without allocating
        let buffer = self.reader.retrying_fill_buf()?;
        let plus_marker = buffer.first().copied();
        // Skip the line, even in the case when the marker is not present (this
        // avoids repeated errors and infinite loops)
        self.reader.skip_until(b'\n')?;

        if plus_marker != Some(b'+') {
            return Err(IOError::new(
                ErrorKind::InvalidData,
                format!("Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: {header}"),
            ));
        }

        // Read quality line directly into Vec of known capacity, along with
        // space to hold potential \r character before it gets popped. Use
        // wrapping_add for efficiency over saturating_add, which will almost
        // surely never wrap (but if it does, it does not cause a logic error)
        let quality = self.reader.read_line_owned(sequence.len().wrapping_add(1))?;

        if quality.len() != sequence.len() {
            if quality.is_empty() {
                return Err(IOError::new(
                    ErrorKind::InvalidData,
                    format!("Missing FASTQ quality scores! See header: {header}"),
                ));
            }

            return Err(IOError::new(
                ErrorKind::InvalidData,
                format!(
                    "Sequence and quality score length mismatch ({s} ≠ {q})! See: {header}",
                    s = sequence.len(),
                    q = quality.len(),
                ),
            ));
        }

        let quality = QualityScores::try_from(quality)?;

        Ok(Some(FastQ {
            header,
            sequence,
            quality,
        }))
    }
}
