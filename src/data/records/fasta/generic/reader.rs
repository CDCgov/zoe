use crate::{
    DEFAULT_SIMD_LANES,
    data::{
        err::ResultWithErrorContext,
        fasta::generic::Fasta,
        vec_types::{ChopLineBreak, StripLineBreak},
    },
    search::ByteSplitIter,
    unwrap_or_return_some_err,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, ErrorKind},
    path::Path,
};

/// A buffered reader for reading
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.
///
/// The sequence type of the resulting [`Fasta`] record is `Vec<u8>`.
///
/// ## Parameters
///
/// `R`: The type of data being read, which is wrapped in a [`BufReader`] before
///  use
#[derive(Debug)]
pub struct FastaReader<R: std::io::Read> {
    reader:       BufReader<R>,
    buffer:       Vec<u8>,
    first_record: bool,
}

impl<R: std::io::Read> FastaReader<R> {
    /// Creates an iterator over FASTA data, wrapping the input in a buffered
    /// reader.
    ///
    /// Unlike [`from_readable`], this does not allocate or read any data
    /// initially. It also allows for empty input, in which case the resulting
    /// iterator is empty.
    ///
    /// [`from_readable`]: FastaReader::from_readable
    pub fn new(inner: R) -> Self {
        FastaReader {
            reader:       BufReader::new(inner),
            buffer:       Vec::new(),
            first_record: true,
        }
    }

    /// Creates an iterator over FASTA data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    ///
    /// [`Read`]: std::io::Read
    pub fn from_readable(read: R) -> std::io::Result<Self> {
        FastaReader::from_bufreader(BufReader::new(read))
    }

    /// Creates an iterator over FASTA data from a `BufReader`.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_bufreader(mut reader: BufReader<R>) -> std::io::Result<Self> {
        if reader.fill_buf()?.is_empty() {
            return Err(std::io::Error::new(ErrorKind::InvalidData, "No FASTA data was found!"));
        }

        Ok(FastaReader {
            reader,
            buffer: Vec::new(),
            first_record: true,
        })
    }

    fn get_error(msg: &str, header: Option<&str>) -> std::io::Result<Fasta> {
        if let Some(header) = header {
            Err(std::io::Error::new(
                ErrorKind::InvalidData,
                format!("{msg} See header: {header}"),
            ))
        } else {
            Err(std::io::Error::new(ErrorKind::InvalidData, msg))
        }
    }

    /// Read the first record, ensuring that the file is not slurped when no `>`
    /// is present.
    ///
    /// If `Some(Ok(_))` is returned, then the leading `>` will already be
    /// consumed for the next record if present, and the buffer will contain the
    /// last sequence and trailing '>' if present).
    fn read_first_record(&mut self) -> Option<std::io::Result<Fasta>> {
        self.first_record = false;

        loop {
            let bytes = unwrap_or_return_some_err!(self.reader.read_until(b'\n', &mut self.buffer));
            if bytes == 0 {
                return Some(Self::get_error("No FASTA data found!", None));
            }

            if let Some(mut header_bytes) = self.buffer.strip_prefix(b">") {
                header_bytes = header_bytes.strip_line_break();

                if header_bytes.is_empty() {
                    return Some(Self::get_error("Missing FASTA header!", None));
                }

                let header = String::from_utf8_lossy(header_bytes).into_owned();

                if header_bytes.contains(&b'>') {
                    return Some(Self::get_error(
                        "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a header!",
                        Some(&header),
                    ));
                }

                self.buffer.clear();
                unwrap_or_return_some_err!(self.reader.read_until(b'>', &mut self.buffer));

                let split = self.buffer.lines_ascii::<{ DEFAULT_SIMD_LANES }>();
                let mut sequence = Vec::with_capacity(split.remaining_len());
                for s in split {
                    sequence.extend_from_slice(s);
                }

                if sequence.ends_with(b">") {
                    sequence.pop();
                }

                if sequence.is_empty() {
                    return Some(Self::get_error("Missing FASTA sequence!", Some(&header)));
                }

                // Check to make sure we read the full sequence
                if !self.buffer.ends_with(b"\n>") {
                    if self.buffer.ends_with(b">") {
                        return Some(Self::get_error(
                            "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a sequence!",
                            Some(&header),
                        ));
                    }
                    // We have finished iteration
                    self.buffer.clear();
                }

                return Some(Ok(Fasta { header, sequence }));
            } else if self.buffer.iter().all(u8::is_ascii_whitespace) {
                // Clear the whitespace
                self.buffer.clear();
            } else {
                return Some(Self::get_error("The FASTA file must start with a '>' symbol!", None));
            }
        }
    }
}

impl FastaReader<std::fs::File> {
    /// Deprecated, use [`from_path`] instead.
    ///
    /// [`from_path`]: FastaReader::from_path
    #[allow(clippy::missing_errors_doc)]
    #[deprecated(
        since = "0.0.27",
        note = "please use `from_path` instead. This function will be removed in v0.0.29"
    )]
    pub fn from_filename<P>(path: P) -> std::io::Result<FastaReader<File>>
    where
        P: AsRef<Path>, {
        Self::from_path(path)
    }

    /// Creates an iterator over the FASTA data contained in a path, using a
    /// buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the path does not exist, if there are insufficient
    /// permissions to read from it, or if it contains no data. The path is
    /// included in the error message.
    pub fn from_path<P>(path: P) -> std::io::Result<FastaReader<File>>
    where
        P: AsRef<Path>, {
        let path = path.as_ref();
        let file = File::open(path).with_path_context("Failed to open path", path)?;
        Ok(Self::from_readable(file).with_path_context("Failed to read data at path", path)?)
    }
}

/// An iterator for buffered reading of
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.
impl<R: std::io::Read> Iterator for FastaReader<R> {
    type Item = std::io::Result<Fasta>;

    fn next(&mut self) -> Option<Self::Item> {
        // Special case logic for first record to ensure we don't slurp file.
        if self.first_record {
            return self.read_first_record();
        }

        // The buffer will always contain the last sequence, except when the end
        // of the file is reached, at which point we clear it.
        if self.buffer.is_empty() {
            return None;
        }

        // Read full record, all the way up to (and including) next '>'
        self.buffer.clear();
        unwrap_or_return_some_err!(self.reader.read_until(b'>', &mut self.buffer));
        let mut split = self.buffer.lines_ascii::<{ DEFAULT_SIMD_LANES }>();

        let Some(header) = split.next() else {
            return Some(Self::get_error("Missing FASTA header!", None));
        };
        let header = String::from_utf8_lossy(header).into_owned();
        if header.is_empty() {
            return Some(Self::get_error("Missing FASTA header!", None));
        }

        let mut sequence = Vec::with_capacity(split.remaining_len());
        for s in split {
            sequence.extend_from_slice(s);
        }
        if sequence.ends_with(b">") {
            sequence.pop();
        }

        if sequence.is_empty() {
            // Terminated due to reaching '>' character
            if let Some(b'>') = self.buffer.last() {
                // Check whether the full header line was actually read
                if self.buffer.contains(&b'\n') {
                    return Some(Self::get_error("Missing FASTA sequence!", Some(&header)));
                }

                // We did not finish reading the header line, so do that and
                // issue error
                unwrap_or_return_some_err!(self.reader.read_until(b'\n', &mut self.buffer));
                self.buffer.chop_line_break();
                let name = String::from_utf8_lossy(&self.buffer).into_owned();
                return Some(Self::get_error(
                    "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a header!",
                    Some(&name),
                ));
            }

            // Terminated due to reaching end of file
            return Some(Self::get_error("Missing FASTA sequence!", Some(&header)));
        }

        // Check to make sure we read the full sequence
        if !self.buffer.ends_with(b"\n>") {
            if self.buffer.ends_with(b">") {
                return Some(Self::get_error(
                    "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a sequence!",
                    Some(&header),
                ));
            }
            // We have finished iteration
            self.buffer.clear();
        }

        Some(Ok(Fasta { header, sequence }))
    }
}
