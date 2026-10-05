use crate::{
    data::{err::ResultWithErrorContext, vec_types::ChopLineBreak},
    prelude::*,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, Error as IOError, ErrorKind, Read},
    path::Path,
};

/// A buffered reader for reading
/// [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format).
///
/// This does not support multiline FASTQ files. In other words, each sequence
/// must be on a single line, and the quality scores must be on a single line.
#[derive(Debug)]
pub struct FastQReader<R: Read> {
    reader: BufReader<R>,
}

impl<R: Read> FastQReader<R> {
    /// Creates an iterator over FASTQ data, wrapping the input in a buffered
    /// reader.
    ///
    /// Unlike [`from_readable`], this does not read any data initially. It also
    /// allows for empty input, in which case the resulting iterator is empty.
    ///
    /// [`from_readable`]: FastQReader::from_readable
    pub fn new(inner: R) -> Self {
        FastQReader {
            reader: BufReader::new(inner),
        }
    }

    /// Creates an iterator over FASTQ data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_readable(read: R) -> std::io::Result<Self> {
        FastQReader::from_bufreader(BufReader::new(read))
    }

    /// Creates an iterator over FASTQ data from a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_bufreader(mut reader: BufReader<R>) -> std::io::Result<Self> {
        if reader.retrying_fill_buf()?.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "No FASTQ data was found!"));
        }

        Ok(FastQReader { reader })
    }
}

impl FastQReader<std::fs::File> {
    /// Creates an iterator over the FASTQ data contained in a path, using a
    /// buffered reader.
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

impl<R: Read> FastQReader<R> {
    /// A helper function for [`FastQReader::next`] that returns a [`Result`]
    /// rather than an [`Option`], to allow for a more readable implementation.
    #[inline]
    fn next_helper(&mut self) -> std::io::Result<Option<FastQ>> {
        // Consume "@" at start of header, or skip line and abort
        let buf = self.reader.retrying_fill_buf()?;

        let Some(marker) = buf.first().copied() else { return Ok(None) };
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
        let mut header = read_line_owned(&mut self.reader)?;
        header.chop_line_break();

        if header.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "Missing FASTQ header!"));
        }

        let header = String::from_utf8(header).map_err(|err| IOError::new(ErrorKind::InvalidData, err))?;

        // Read sequence directly into Vec
        let mut sequence = read_line_owned(&mut self.reader)?;
        sequence.chop_line_break();

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

        // Read quality line directly into Vec of known capacity. Use
        // wrapping_add for efficiency over saturating_add, which will almost
        // surely never wrap (but if it does, it does not cause a logic error)
        let mut quality = Vec::with_capacity(sequence.len().wrapping_add(2));
        self.reader.read_until(b'\n', &mut quality)?;
        quality.chop_line_break();

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

/// An extension trait for [`BufRead`] offering a retrying version of
/// [`BufRead::fill_buf`].
trait BufReadExtension: BufRead {
    /// A version of [`BufRead::fill_buf`] that retries upon
    /// [`ErrorKind::Interrupted`].
    ///
    /// Many higher-level [`BufRead`] methods retry on this error like
    /// [`BufRead::read_until`], so this method allows implementing similar
    /// higher-level functionality.
    fn retrying_fill_buf(&mut self) -> std::io::Result<&[u8]> {
        loop {
            match self.fill_buf() {
                Ok(buf) => return Ok(buf),
                Err(err) if err.kind() == ErrorKind::Interrupted => {}
                Err(err) => return Err(err),
            }
        }
    }
}

impl<R: BufRead> BufReadExtension for R {}

/// Reads a line (including any line break) into an exactly-sized [`Vec`] when
/// it is already fully buffered, otherwise falls back to [`read_until`].
///
/// [`read_until`]: BufRead::read_until
#[inline]
fn read_line_owned<R: Read>(reader: &mut BufReader<R>) -> std::io::Result<Vec<u8>> {
    let buf = reader.retrying_fill_buf()?;

    if let Some(i) = buf.find_byte(b'\n') {
        let line = buf[..=i].to_vec();
        reader.consume(i + 1);
        return Ok(line);
    }

    let mut line = Vec::new();
    reader.read_until(b'\n', &mut line)?;
    Ok(line)
}
