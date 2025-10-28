use crate::{
    DEFAULT_SIMD_LANES,
    data::{
        id_types::FastaIDs,
        records::RecordReader,
        types::{
            amino_acids::AminoAcids,
            nucleotides::{self, Nucleotides, ToDNA, Translate},
        },
        validation::CheckSequence,
        vec_types::{ChopLineBreak, StripLineBreak},
    },
    search::ByteSplitIter,
    unwrap_or_return_some_err,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, Error as IOError, ErrorKind},
    path::Path,
};

#[cfg(test)]
mod test;

/// Provides a container struct for data from a generic
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaSeq {
    pub name:     String,
    pub sequence: Vec<u8>,
}

/// Similar to [`FastaSeq`] but assumes that the `sequence` contains valid
/// [`Nucleotides`].
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaNT {
    pub name:     String,
    pub sequence: Nucleotides,
}

/// Similar to [`FastaSeq`] but assumes that the `sequence` contains valid
/// [`AminoAcids`].
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaAA {
    pub name:     String,
    pub sequence: AminoAcids,
}

/// Structure for buffered reading of `FASTA` files.
#[derive(Debug)]
pub struct FastaReader<R: std::io::Read> {
    reader:       std::io::BufReader<R>,
    buffer:       Vec<u8>,
    first_record: bool,
}

impl FastaSeq {
    /// Reverse complements the sequence stored in the struct using a new
    /// buffer.
    pub fn reverse_complement(&mut self) {
        self.sequence = nucleotides::reverse_complement(&self.sequence);
    }

    /// Recodes to uppercase IUPAC DNA with corrected gaps, otherwise
    /// mapping to `N`. Returns [`FastaNT`].
    #[inline]
    #[must_use]
    pub fn recode_to_dna(self) -> FastaNT {
        FastaNT {
            name:     self.name,
            sequence: self.sequence.recode_to_dna(),
        }
    }

    /// Filters and recodes to uppercase IUPAC DNA with corrected gaps. Returns
    /// [`FastaNT`].
    #[inline]
    #[must_use]
    pub fn filter_to_dna(self) -> FastaNT {
        FastaNT {
            name:     self.name,
            sequence: self.sequence.filter_to_dna(),
        }
    }

    /// Filters and recodes to uppercase IUPAC DNA without gaps. Returns
    /// [`FastaNT`].
    #[inline]
    #[must_use]
    pub fn filter_to_dna_unaligned(self) -> FastaNT {
        FastaNT {
            name:     self.name,
            sequence: self.sequence.filter_to_dna_uanligned(),
        }
    }

    /// For an annotated `FASTA` with format `id{annotation}` returns a tuple
    /// of the id and annotated taxon.
    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, &str)> {
        self.name.get_id_taxon()
    }

    /// Translates the stored [`Vec<u8>`] to [`AminoAcids`] using a new buffer.
    ///
    /// See [`translate_sequence`] for more details.
    ///
    /// [`translate_sequence`]: nucleotides::translate_sequence
    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     self.name,
            sequence: AminoAcids(nucleotides::translate_sequence(&self.sequence)),
        }
    }
}

impl FastaNT {
    /// Reverse complements the sequence stored in the struct using a new buffer.
    #[inline]
    pub fn reverse_complement(&mut self) {
        self.sequence.make_reverse_complement();
    }

    /// Translates the stored [`Nucleotides`] to [`AminoAcids`] using a new buffer.
    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     self.name,
            sequence: self.sequence.translate(),
        }
    }

    /// For an annotated FASTA file with format `id{annotation}` returns a tuple
    /// of the id and annotated taxon.
    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, &str)> {
        self.name.get_id_taxon()
    }
}

/// A struct for containing an [`FastaNT`] + owned `taxon` [`String`].
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct FastaNTAnnot {
    pub name:     String,
    pub sequence: Nucleotides,
    pub taxon:    String,
}

/// Fallibly converts from a generic [`FastaSeq`] to a [`FastaNTAnnot`], failing
/// if no annotation is found. DNA is filtered in the process.
impl TryFrom<FastaSeq> for FastaNTAnnot {
    type Error = std::io::Error;
    fn try_from(fa: FastaSeq) -> Result<Self, Self::Error> {
        if let Some((id, taxon)) = fa.name.get_id_taxon() {
            Ok(FastaNTAnnot {
                name:     id.to_string(),
                sequence: fa.sequence.filter_to_dna(),
                taxon:    taxon.to_string(),
            })
        } else {
            Err(std::io::Error::new(
                ErrorKind::InvalidData,
                format!("No taxon for: {id}", id = fa.name),
            ))
        }
    }
}

/// Allows converting from [`FastaSeq`] to [`FastaNT`] without checks or
/// filtering.
impl From<FastaSeq> for FastaNT {
    fn from(record: FastaSeq) -> Self {
        FastaNT {
            name:     record.name,
            sequence: record.sequence.into(),
        }
    }
}

impl FastaAA {
    /// For an annotated `FASTA` with format `id{annotation}` returns a tuple
    /// of the id and annotated taxon.
    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, &str)> {
        self.name.get_id_taxon()
    }
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
            reader:       std::io::BufReader::new(inner),
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
        FastaReader::from_bufreader(std::io::BufReader::new(read))
    }

    /// Creates an iterator over FASTA data from a `BufReader`.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_bufreader(mut reader: BufReader<R>) -> std::io::Result<Self> {
        if reader.fill_buf()?.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "No FASTA data was found!"));
        }

        Ok(FastaReader {
            reader,
            buffer: Vec::new(),
            first_record: true,
        })
    }

    fn get_error(msg: &str, header: Option<&str>) -> std::io::Result<FastaSeq> {
        if let Some(header) = header {
            Err(IOError::new(ErrorKind::InvalidData, format!("{msg} See header: {header}")))
        } else {
            Err(IOError::new(ErrorKind::InvalidData, msg))
        }
    }

    /// Read the first record, ensuring that the file is not slurped when no `>`
    /// is present.
    ///
    /// If `Some(Ok(_))` is returned, then the leading `>` will already be
    /// consumed for the next record if present, and the buffer will contain the
    /// last sequence and trailing '>' if present).
    fn read_first_record(&mut self) -> Option<std::io::Result<FastaSeq>> {
        self.first_record = false;

        loop {
            let bytes = unwrap_or_return_some_err!(self.reader.read_until(b'\n', &mut self.buffer));
            if bytes == 0 {
                return Some(Self::get_error("No FASTA data found!", None));
            }

            if let Some(mut header) = self.buffer.strip_prefix(b">") {
                header = header.strip_line_break();

                if header.is_empty() {
                    return Some(Self::get_error("Missing FASTA header!", None));
                }

                let name = String::from_utf8_lossy(header).into_owned();

                if header.contains(&b'>') {
                    return Some(Self::get_error(
                        "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a header!",
                        Some(&name),
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
                    return Some(Self::get_error("Missing FASTA sequence!", Some(&name)));
                }

                // Check to make sure we read the full sequence
                if !self.buffer.ends_with(b"\n>") {
                    if self.buffer.ends_with(b">") {
                        return Some(Self::get_error(
                            "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a sequence!",
                            Some(&name),
                        ));
                    }
                    // We have finished iteration
                    self.buffer.clear();
                }

                return Some(Ok(FastaSeq { name, sequence }));
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
    /// Reads a fasta file into an iterator backed by a buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if file or permissions do not exist, or if the file is
    /// empty. The file path is included in the error message.
    pub fn from_filename<P>(filename: P) -> std::io::Result<FastaReader<File>>
    where
        P: AsRef<Path>, {
        let file = Self::open_nonempty_file(filename)?;
        Ok(FastaReader::new(file))
    }
}

impl<R: std::io::Read> RecordReader for FastaReader<R> {
    const RECORD_NAME: &str = "FASTA";
}

/// An iterator for buffered reading of
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.
impl<R: std::io::Read> Iterator for FastaReader<R> {
    type Item = std::io::Result<FastaSeq>;

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

        let Some(name_line) = split.next() else {
            return Some(Self::get_error("Missing FASTA header!", None));
        };
        let name = String::from_utf8_lossy(name_line).into_owned();
        if name.is_empty() {
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
                    return Some(Self::get_error("Missing FASTA sequence!", Some(&name)));
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
            return Some(Self::get_error("Missing FASTA sequence!", Some(&name)));
        }

        // Check to make sure we read the full sequence
        if !self.buffer.ends_with(b"\n>") {
            if self.buffer.ends_with(b">") {
                return Some(Self::get_error(
                    "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a sequence!",
                    Some(&name),
                ));
            }
            // We have finished iteration
            self.buffer.clear();
        }

        Some(Ok(FastaSeq { name, sequence }))
    }
}

impl std::fmt::Display for FastaSeq {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.sequence.is_ascii_simd::<16>() {
            // SAFETY: we just checked it is ASCII using our fast SIMD function.
            // ASCII is valid UTF8.
            write!(f, ">{}\n{}\n", self.name, unsafe {
                std::str::from_utf8_unchecked(&self.sequence)
            })
        } else {
            write!(f, ">{}\n{}\n", self.name, String::from_utf8_lossy(&self.sequence))
        }
    }
}

impl std::fmt::Display for FastaNT {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, self.sequence)
    }
}

impl std::fmt::Display for FastaAA {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, self.sequence)
    }
}
