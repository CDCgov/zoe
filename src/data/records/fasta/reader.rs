use crate::{
    DEFAULT_SIMD_LANES,
    data::{
        err::ResultWithErrorContext,
        fasta::Fasta,
        records::reader_utils::{BufReadExtension, FusedBufReader},
    },
    search::{ByteSubstring, position_by_byte3},
};
use std::{
    fs::File,
    io::{BufRead, BufReader, ErrorKind, Read},
    iter::FusedIterator,
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
pub struct FastaReader<R: Read> {
    reader:              FusedBufReader<R>,
    first_record:        bool,
    parse_header_symbol: bool,
}

impl<R: Read> FastaReader<R> {
    /// Creates an iterator over FASTA data, wrapping the input in a buffered
    /// reader.
    ///
    /// Unlike [`from_readable`], this does not read any data initially. It also
    /// allows for empty input, in which case the resulting iterator is empty.
    ///
    /// The iterator returns `None` once `inner.read()` returns 0 (signalling
    /// EOF). After this, the iterator will continue to return `None` regardless
    /// of whether `inner.read()` could return new bytes.
    ///
    /// [`from_readable`]: FastaReader::from_readable
    pub fn new(inner: R) -> Self {
        FastaReader {
            reader:              FusedBufReader::new(inner),
            first_record:        true,
            parse_header_symbol: false,
        }
    }

    /// Creates an iterator over FASTA data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// The iterator returns `None` once `inner.read()` returns 0 (signalling
    /// EOF). After this, the iterator will continue to return `None` regardless
    /// of whether `inner.read()` could return new bytes.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    ///
    /// [`Read`]: std::io::Read
    pub fn from_readable(read: R) -> std::io::Result<Self> {
        FastaReader::from_bufreader(BufReader::new(read))
    }

    /// Creates an iterator over FASTA data from a [`BufReader`].
    ///
    /// The iterator returns `None` once `reader.fill_buf()` returns an empty
    /// slice (signalling EOF). After this, the iterator will continue to return
    /// `None` regardless of whether `reader.fill_buf()` could return new bytes.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_bufreader(mut reader: BufReader<R>) -> std::io::Result<Self> {
        if reader.fill_buf()?.is_empty() {
            return Err(std::io::Error::new(ErrorKind::InvalidData, "No FASTA data was found!"));
        }

        Ok(FastaReader {
            reader:              FusedBufReader::from_bufreader(reader),
            first_record:        true,
            parse_header_symbol: false,
        })
    }

    fn get_error(msg: &str, header: Option<&str>) -> std::io::Error {
        if let Some(header) = header {
            std::io::Error::new(ErrorKind::InvalidData, format!("{msg} See header: {header}"))
        } else {
            std::io::Error::new(ErrorKind::InvalidData, msg)
        }
    }

    fn consume_first_header_symbol(&mut self) -> std::io::Result<Option<()>> {
        let mut at_line_start = true;

        // Repeatedly fill the buffer, consuming any whitespace encountered.
        // Once a non-whitespace character is encountered, the loop breaks with
        // true if > is found at the start of a line, and false otherwise.
        let header_symbol_found = loop {
            let buf = self.reader.retrying_fill_buf()?;

            // If the buffer is empty, return Ok(None)
            let Some(&marker) = buf.first() else { return Ok(None) };

            // Terminate loop if we have found ">"
            if marker == b'>' {
                // This only represents a valid header if it is at the start of
                // a line
                break at_line_start;
            }

            // If we can find a non-whitespace byte, then the loop will
            // terminate
            if let Some(idx) = buf.iter().position(|c| !c.is_ascii_whitespace()) {
                let header_symbol_found = if buf[idx] == b'>' {
                    // Check if we are at start of line. Validity: idx > 0 since
                    // buf did not start with ">"
                    buf[idx - 1] == b'\n'
                } else {
                    // Found a non-whitespace byte other than ">"
                    false
                };

                // Consume all whitespace bytes
                self.reader.consume(idx);

                break header_symbol_found;
            }

            // buf.last is known to be Some since buf is non-empty, but may as
            // well avoid an unwrap
            at_line_start = buf.last() == Some(&b'\n');

            // Consume the buffer, which is all whitespace
            let len = buf.len();
            self.reader.consume(len);
        };

        if header_symbol_found {
            // Consume ">"
            self.reader.consume(1);
            Ok(Some(()))
        } else {
            // Skip the entire line that does not start with ">"
            self.reader.skip_until(b'\n')?;
            Err(Self::get_error("The FASTA file must start with a '>' symbol!", None))
        }
    }
}

impl FastaReader<File> {
    /// Creates an iterator over the FASTA data contained in a path, using a
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
    pub fn from_path<P>(path: P) -> std::io::Result<FastaReader<File>>
    where
        P: AsRef<Path>, {
        let path = path.as_ref();
        let file = File::open(path).with_path_context("Failed to open path", path)?;
        Ok(Self::from_readable(file).with_path_context("Failed to read data at path", path)?)
    }
}

impl<R: Read> FastaReader<R> {
    fn next_helper(&mut self) -> std::io::Result<Option<Fasta>> {
        if self.reader.is_at_eof() {
            return Ok(None);
        }

        // Consume leading ">" for first record
        if self.first_record {
            self.first_record = false;
            match self.consume_first_header_symbol()? {
                Some(()) => {}
                None => return Ok(None),
            }
        } else if self.parse_header_symbol {
            // Consume ">" at start of header, or skip line and abort
            let buf = self.reader.retrying_fill_buf()?;
            let Some(marker) = buf.first().copied() else {
                // EOF reached
                return Ok(None);
            };
            if marker != b'>' {
                self.reader.skip_until(b'\n')?;
                return Err(Self::get_error("Missing '>' symbol at header line beginning!", None));
            }
            self.reader.consume(1);
            self.parse_header_symbol = false;
        }

        // Read header directly into Vec
        let header = self.reader.read_line_owned(0)?;
        if header.is_empty() {
            self.parse_header_symbol = true;
            return Err(Self::get_error("Missing FASTA header!", None));
        }
        let header = String::from_utf8_lossy_owned(header);
        if header.find_byte(b'>').is_some() {
            self.parse_header_symbol = true;
            return Err(Self::get_error(
                "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a header!",
                Some(&header),
            ));
        }

        let mut sequence = Vec::new();
        let mut at_line_start = true;
        loop {
            let buf = self.reader.retrying_fill_buf()?;
            if buf.is_empty() {
                break;
            }

            if let Some(i) = position_by_byte3::<{ DEFAULT_SIMD_LANES }>(buf, b'\r', b'\n', b'>') {
                let byte = buf[i];

                if byte == b'>' {
                    if at_line_start && i == 0 {
                        self.reader.consume(1);
                        break;
                    }

                    self.reader.skip_until_with_prefix(b'>', b'\n')?;
                    return Err(Self::get_error(
                        "No instance of '>' can occur in the sequence!",
                        Some(&header),
                    ));
                }

                sequence.extend_from_slice(&buf[..i]);
                self.reader.consume(i + 1);
                at_line_start = byte == b'\n';
            } else {
                let len = buf.len();
                sequence.extend_from_slice(buf);
                self.reader.consume(len);
                at_line_start = false;
            }
        }

        if sequence.is_empty() {
            return Err(Self::get_error("Missing FASTA sequence!", Some(&header)));
        }

        Ok(Some(Fasta { header, sequence }))
    }
}

impl<R: Read> Iterator for FastaReader<R> {
    type Item = std::io::Result<Fasta>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_helper().transpose()
    }
}

impl<R: Read> FusedIterator for FastaReader<R> {}
