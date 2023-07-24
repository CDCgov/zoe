use super::{
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
    pub name:     Vec<u8>,
    pub sequence: Nucleotides,
}

pub struct FastaAA {
    pub name:     Vec<u8>,
    pub sequence: AminoAcids,
}

pub struct FastaReader<R: std::io::Read> {
    pub fasta_reader: std::io::BufReader<R>,
    pub fasta_buffer: Vec<u8>,
}

impl FastaSeq {
    pub fn reverse_complement(&mut self) {
        self.sequence = reverse_complement(&self.sequence);
    }

    #[must_use]
    pub fn retain_unaligned_dna_uc(mut self) -> FastaNT {
        self.sequence.retain_by_recoding(TO_UNALIGNED_DNA_UC);
        FastaNT {
            name:     self.name,
            sequence: Nucleotides(self.sequence),
        }
    }

    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     self.name,
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
            sequence: AminoAcids(self.sequence.translate()),
        }
    }
}

impl From<FastaSeq> for FastaNT {
    fn from(record: FastaSeq) -> Self {
        FastaNT {
            name:     record.name,
            sequence: record.sequence.into(),
        }
    }
}

//impl FastAA {}

impl<R: std::io::Read> FastaReader<R> {
    pub fn new(inner: R) -> Self {
        FastaReader {
            fasta_reader: std::io::BufReader::new(inner),
            fasta_buffer: Vec::new(),
        }
    }
}

impl FastaReader<std::fs::File> {
    /// Reads a fasta file into an iterator backed by a buffered reader.
    ///
    /// # Errors
    ///
    /// Will return `Err` if file or permissions or the like do not exist.
    pub fn from_filename<P>(filename: P) -> Result<FastaReader<File>, std::io::Error>
    where
        P: AsRef<Path>, {
        let file = File::open(filename)?;
        Ok(FastaReader::new(file))
    }
}

impl<R: std::io::Read> Iterator for FastaReader<R> {
    type Item = std::io::Result<FastaSeq>;

    fn next(&mut self) -> Option<Self::Item> {
        self.fasta_buffer.clear();

        let bytes = match self.fasta_reader.read_until(b'>', &mut self.fasta_buffer) {
            Ok(b) => b,
            Err(e) => return Some(Err(e)),
        };

        match bytes {
            0 => None,
            1 => self.next(),
            _ => {
                if self.fasta_buffer.ends_with(b">") {
                    self.fasta_buffer.pop();
                }

                let mut lines = self.fasta_buffer.split(|x| *x == b'\n' || *x == b'\r');
                let name = match lines.next() {
                    Some(h) if !h.is_empty() => h.to_vec(),
                    _ => return Some(Err(IOError::new(ErrorKind::InvalidData, "Missing FASTA header!"))),
                };

                let sequence: Vec<u8> = lines.flatten().copied().collect();

                if sequence.is_empty() {
                    Some(Err(IOError::new(ErrorKind::InvalidData, "Missing sequence data!")))
                } else {
                    Some(Ok(FastaSeq { name, sequence }))
                }
            }
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
        write!(
            f,
            ">{}\n{}\n",
            String::from_utf8_lossy(&self.name),
            String::from_utf8_lossy(&self.sequence.0)
        )
    }
}

impl std::fmt::Display for FastaAA {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            ">{}\n{}\n",
            String::from_utf8_lossy(&self.name),
            String::from_utf8_lossy(&self.sequence.0)
        )
    }
}
