use super::types::{nucleotides::Nucleotides, phred::QualityScores};
use std::fs::File;
use std::io::BufRead;

pub struct FastQ {
    pub header: String,
    pub sequence: Nucleotides,
    pub quality: QualityScores,
}

impl std::fmt::Display for FastQ {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}\n{}\n+{}\n", self.header, self.sequence, self.quality)
    }
}

impl FastQ {
    pub fn reverse_complement(&mut self) {
        self.sequence = self.sequence.reverse_complement();
        self.quality = self.quality.reverse();
    }
}

pub struct FastQReader<R: std::io::Read> {
    pub fastq_reader: std::io::BufReader<R>,
    pub fastq_buffer: Vec<u8>,
}

impl<R: std::io::Read> FastQReader<R> {
    pub fn new(inner: R) -> Self {
        FastQReader {
            fastq_reader: std::io::BufReader::new(inner),
            fastq_buffer: Vec::new(),
        }
    }
}

impl FastQReader<std::fs::File> {
    /// # Errors
    ///
    /// Will return `Err` if file or permissions do not exist.
    pub fn from_filename(filename: &str) -> Result<FastQReader<File>, String> {
        match File::open(filename) {
            Err(why) => Err(format!("Couldn't open fasta file '{filename}': {why}")),
            Ok(file) => Ok(FastQReader::new(file)),
        }
    }
}

// Consider making fallible
impl<R: std::io::Read> Iterator for FastQReader<R> {
    type Item = FastQ;

    fn next(&mut self) -> Option<Self::Item> {
        self.fastq_buffer.clear();

        // Read HEADER line
        let bytes = self
            .fastq_reader
            .read_until(b'\n', &mut self.fastq_buffer)
            .unwrap_or_else(|e| {
                eprintln!("FASTQ read header failed: {e}\n");
                0
            });
        if bytes == 0 {
            return None;
        }

        let header = if self.fastq_buffer.len() > 1 {
            if self.fastq_buffer.ends_with(b"\n") {
                self.fastq_buffer.pop();
            }

            match String::from_utf8(self.fastq_buffer.clone()) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("FASTQ header failed to parse:\n{e}\n");
                    return None;
                }
            }
        } else {
            String::new()
        };
        self.fastq_buffer.clear();

        // Read SEQUENCE line
        let bytes = self
            .fastq_reader
            .read_until(b'\n', &mut self.fastq_buffer)
            .unwrap_or_else(|e| {
                eprintln!("FASTQ read failed: {e}\n");
                0
            });
        if bytes == 0 {
            return None;
        }

        let sequence = if self.fastq_buffer.len() > 1 {
            if self.fastq_buffer.ends_with(b"\n") {
                self.fastq_buffer.pop();
            }

            Nucleotides(
                self.fastq_buffer
                    .iter()
                    .copied()
                    .map(|c| c.to_ascii_uppercase())
                    .collect(),
            )
        } else {
            Nucleotides::new()
        };
        self.fastq_buffer.clear();

        // Read "+" line
        let bytes = self
            .fastq_reader
            .read_until(b'\n', &mut self.fastq_buffer)
            .unwrap_or_else(|e| {
                eprintln!("FASTQ read failed: {e}\n");
                0
            });
        if bytes == 0 {
            return None;
        }
        self.fastq_buffer.clear();

        // Read QUALITY line
        let bytes = self
            .fastq_reader
            .read_until(b'\n', &mut self.fastq_buffer)
            .unwrap_or_else(|e| {
                eprintln!("FASTQ read failed: {e}\n");
                0
            });
        if bytes == 0 {
            return None;
        }

        let quality = if self.fastq_buffer.len() > 1 {
            if self.fastq_buffer.ends_with(b"\n") {
                self.fastq_buffer.pop();
            }

            QualityScores(
                self.fastq_buffer
                    .iter()
                    .copied()
                    .map(|c| c.to_ascii_uppercase())
                    .collect(),
            )
        } else {
            QualityScores::new()
        };

        Some(FastQ {
            header,
            sequence,
            quality,
        })
    }
}
