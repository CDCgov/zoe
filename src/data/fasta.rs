use crate::{
    data::{
        convert::ToDNA,
        id_types::FastaIDs,
        mappings::TO_UNALIGNED_DNA_UC,
        types::{
            amino_acids::AminoAcids,
            nucleotides::{reverse_complement, Nucleotides},
        },
        vec_types::{CheckSequence, ValidateSequence},
    },
    search::ByteSplitIter,
};
use std::{
    fs::File,
    io::{BufRead, Error as IOError, ErrorKind},
    path::Path,
};

/// Provides a container struct for data from a generic
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file.
#[derive(Debug)]
pub struct FastaSeq {
    pub name:     String,
    pub sequence: Vec<u8>,
}

/// Similar to [`FastaSeq`] but assumes that the `sequence` contains valid [Nucleotides].
pub struct FastaNT {
    pub name:     String,
    pub sequence: Nucleotides,
}

/// Similar to [`FastaSeq`] but assumes that the `sequence` contains valid [`AminoAcids`].
pub struct FastaAA {
    pub name:     String,
    pub sequence: AminoAcids,
}

/// Structure for buffered reading of `FASTA` files.
pub struct FastaReader<R: std::io::Read> {
    fasta_reader: std::io::BufReader<R>,
    fasta_buffer: Vec<u8>,
}

impl FastaSeq {
    /// Reverse complements the sequence stored in the struct using a new buffer.
    pub fn reverse_complement(&mut self) {
        self.sequence = reverse_complement(&self.sequence);
    }

    /// Retains only valid DNA characters, removing alignment characters as well.
    #[must_use]
    pub fn to_dna_unaligned(mut self) -> FastaNT {
        self.sequence.retain_by_recoding(TO_UNALIGNED_DNA_UC);
        FastaNT {
            name:     self.name,
            sequence: Nucleotides(self.sequence),
        }
    }

    /// Trims all invalid state from the sequence and validates name as UTF-8.
    /// Returns `FastaNT`.
    #[inline]
    #[must_use]
    pub fn filter_to_dna(self) -> FastaNT {
        FastaNT {
            name:     self.name,
            sequence: self.sequence.filter_to_dna(),
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
    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     self.name,
            sequence: AminoAcids(super::types::nucleotides::translate_sequence(&self.sequence)),
        }
    }
}

impl FastaNT {
    /// Reverse complements the sequence stored in the struct using a new buffer.
    #[inline]
    pub fn reverse_complement(&mut self) {
        self.sequence = self.sequence.reverse_complement();
    }

    /// Translates the stored [Nucleotides] to [`AminoAcids`] using a new buffer.
    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     self.name,
            sequence: self.sequence.translate(),
        }
    }

    /// For an annotated `FASTA` with format `id{annotation}` returns a tuple
    /// of the id and annotated taxon.
    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, &str)> {
        self.name.get_id_taxon()
    }
}

/// A struct for containing an [`FastaNT`] + owned `taxon` [String].
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
    /// Create a new `FastaReader` for buffered reading.
    ///
    /// # Errors
    ///
    /// Will return `Err` if the first record has an invalid start or otherwise cannot be read.
    pub fn new(inner: R) -> Result<FastaReader<R>, std::io::Error> {
        let mut r = FastaReader {
            fasta_reader: std::io::BufReader::new(inner),
            fasta_buffer: Vec::new(),
        };

        let bytes = r.fasta_reader.read_until(b'>', &mut r.fasta_buffer)?;
        match bytes {
            0 => Err(IOError::new(ErrorKind::Other, "No FASTA data found.")),
            1 => Ok(r),
            _ => {
                if r.fasta_buffer.ends_with(b"\n>") {
                    Ok(r)
                } else {
                    Err(IOError::new(
                        ErrorKind::InvalidData,
                        "FASTA records must start with the '>' symbol!",
                    ))
                }
            }
        }
    }
}

impl FastaReader<std::fs::File> {
    /// Reads a fasta file into an iterator backed by a buffered reader.
    ///
    /// # Errors
    ///
    /// Will return `Err` if file or permissions or the like do not exist or if the first record is invalid.
    pub fn from_filename<P>(filename: P) -> Result<FastaReader<File>, std::io::Error>
    where
        P: AsRef<Path>, {
        let file = File::open(filename)?;
        FastaReader::new(file)
    }
}

/// An iterator for buffered reading of
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.
impl<R: std::io::Read> Iterator for FastaReader<R> {
    type Item = std::io::Result<FastaSeq>;

    fn next(&mut self) -> Option<Self::Item> {
        let valid_start = self.fasta_buffer.ends_with(b"\n>") || self.fasta_buffer == b">";
        self.fasta_buffer.clear();

        let bytes = match self.fasta_reader.read_until(b'>', &mut self.fasta_buffer) {
            Ok(b) => b,
            Err(e) => return Some(Err(e)),
        };

        if bytes > 0 {
            if !valid_start {
                return Some(Err(IOError::new(
                    ErrorKind::InvalidData,
                    "FASTA records must start with the '>' symbol!",
                )));
            }

            let mut split = self.fasta_buffer.lines_ascii::<32>();
            let name = if let Some(h) = split.next()
                && !h.is_empty()
            {
                String::from_utf8_lossy(h).into_owned()
            } else {
                return Some(Err(IOError::new(ErrorKind::InvalidData, "Missing FASTA header!")));
            };

            let mut sequence = Vec::with_capacity(split.remaining_len());
            for s in split {
                sequence.extend_from_slice(s);
            }

            if sequence.ends_with(b">") {
                sequence.pop();
            }

            if sequence.is_empty() {
                Some(Err(IOError::new(ErrorKind::InvalidData, "Missing sequence data!")))
            } else {
                Some(Ok(FastaSeq { name, sequence }))
            }
        } else {
            None
        }
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
        if self.sequence.is_ascii_simd::<16>() {
            // SAFETY: we just checked it is ASCII using our fast SIMD function.
            // ASCII is valid UTF8.
            write!(f, ">{}\n{}\n", self.name, unsafe {
                std::str::from_utf8_unchecked(&self.sequence.0)
            })
        } else {
            write!(f, ">{}\n{}\n", self.name, String::from_utf8_lossy(&self.sequence.0))
        }
    }
}

impl std::fmt::Display for FastaAA {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.sequence.is_ascii_simd::<16>() {
            // SAFETY: we just checked it is ASCII using our fast SIMD function.
            // ASCII is valid UTF8.
            write!(f, ">{}\n{}\n", self.name, unsafe {
                std::str::from_utf8_unchecked(&self.sequence.0)
            })
        } else {
            write!(f, ">{}\n{}\n", self.name, String::from_utf8_lossy(&self.sequence.0))
        }
    }
}
