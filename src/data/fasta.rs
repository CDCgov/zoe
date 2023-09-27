use super::{
    convert::ToDNA,
    mappings::TO_UNALIGNED_DNA_UC,
    types::{
        amino_acids::AminoAcids,
        nucleotides::{reverse_complement, Nucleotides},
    },
    vec_types::ValidateSequence,
};
use std::io::{BufRead, Error as IOError, ErrorKind};
use std::{fs::File, path::Path};

#[derive(Debug)]
pub struct FastaSeq {
    pub name:     Vec<u8>,
    pub sequence: Vec<u8>,
}

pub struct FastaNT {
    pub name:     String,
    pub sequence: Nucleotides,
}

pub struct FastaAA {
    pub name:     String,
    pub sequence: AminoAcids,
}

pub struct FastaReader<R: std::io::Read> {
    fasta_reader: std::io::BufReader<R>,
    fasta_buffer: Vec<u8>,
}

impl FastaSeq {
    pub fn reverse_complement(&mut self) {
        self.sequence = reverse_complement(&self.sequence);
    }

    #[must_use]
    pub fn to_dna_unaligned(mut self) -> FastaNT {
        self.sequence.retain_by_recoding(TO_UNALIGNED_DNA_UC);
        FastaNT {
            name:     String::from_utf8_lossy(&self.name).to_string(),
            sequence: Nucleotides(self.sequence),
        }
    }

    /// Trims all invalid state from the sequence and validates name as UTF-8.
    /// Returns `FastaNT`.
    #[inline]
    #[must_use]
    pub fn to_dna(&self) -> FastaNT {
        FastaNT {
            name:     String::from_utf8_lossy(&self.name).to_string(),
            sequence: self.sequence.to_dna(),
        }
    }

    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     String::from_utf8_lossy(&self.name).to_string(),
            sequence: AminoAcids(super::types::nucleotides::translate_sequence(&self.sequence)),
        }
    }
}

impl FastaNT {
    pub fn reverse_complement(&mut self) {
        self.sequence = self.sequence.reverse_complement();
    }

    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     self.name,
            sequence: self.sequence.translate(),
        }
    }
}

impl From<FastaSeq> for FastaNT {
    fn from(record: FastaSeq) -> Self {
        FastaNT {
            name:     String::from_utf8_lossy(&record.name).to_string(),
            sequence: record.sequence.into(),
        }
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

            let mut lines = self.fasta_buffer.split(|x| *x == b'\n' || *x == b'\r');
            let name = match lines.next() {
                Some(h) if !h.is_empty() => h.to_vec(),
                _ => return Some(Err(IOError::new(ErrorKind::InvalidData, "Missing FASTA header!"))),
            };

            let mut sequence: Vec<u8> = lines.flatten().copied().collect();
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
        write!(
            f,
            ">{}\n{}\n",
            String::from_utf8_lossy(&self.name),
            String::from_utf8_lossy(&self.sequence)
        )
    }
}

impl std::fmt::Display for FastaNT {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", &self.name, String::from_utf8_lossy(&self.sequence.0))
    }
}

impl std::fmt::Display for FastaAA {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, String::from_utf8_lossy(&self.sequence.0))
    }
}
