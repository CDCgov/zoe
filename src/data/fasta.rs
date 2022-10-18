use super::types::nucleotides::{reverse_complement, Nucleotides};
use std::fs::File;
use std::io::BufRead;

#[derive(Debug)]
pub struct FastaSeq {
    pub name: Vec<u8>,
    pub sequence: Vec<u8>,
}

pub struct FastaReader<R: std::io::Read> {
    pub fasta_reader: std::io::BufReader<R>,
    pub fasta_buffer: Vec<u8>,
}
pub struct FastaNT {
    pub name: String,
    pub sequence: Nucleotides,
}

// FastaAA & AminoAcids

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

impl FastaSeq {
    pub fn reverse_complement(&mut self) {
        self.sequence = reverse_complement(&self.sequence);
    }
}

impl<R: std::io::Read> FastaReader<R> {
    pub fn new(inner: R) -> Self {
        FastaReader {
            fasta_reader: std::io::BufReader::new(inner),
            fasta_buffer: Vec::new(),
        }
    }
}

impl FastaReader<std::fs::File> {
    /// # Errors
    ///
    /// Will return `Err` if file or permissions do not exist.
    pub fn from_filename(filename: &str) -> Result<FastaReader<File>, String> {
        match File::open(filename) {
            Err(why) => Err(format!("Couldn't open fasta file '{filename}': {why}")),
            Ok(file) => Ok(FastaReader::new(file)),
        }
    }
}

impl<R: std::io::Read> Iterator for FastaReader<R> {
    type Item = FastaSeq;

    fn next(&mut self) -> Option<Self::Item> {
        self.fasta_buffer.clear();

        let bytes = self
            .fasta_reader
            .read_until(b'>', &mut self.fasta_buffer)
            .unwrap_or_else(|e| {
                eprintln!("FASTA read failed: {e}\n");
                0
            });

        match bytes {
            0 => None,
            1 => self.next(),
            _ => {
                if self.fasta_buffer.ends_with(b">") {
                    self.fasta_buffer.pop();
                }

                let mut lines = self.fasta_buffer.split(|x| *x == b'\n' || *x == b'\r');
                let name = match lines.next() {
                    Some(t) => t.to_vec(),
                    None => b"UNKNOWN".to_vec(),
                };

                let sequence: Vec<u8> = lines.flatten().copied().collect();

                if sequence.is_empty() {
                    None
                } else {
                    Some(FastaSeq { name, sequence })
                }
            }
        }
    }
}
