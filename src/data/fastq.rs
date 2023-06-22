use super::types::{nucleotides::Nucleotides, phred::QualityScores};
use std::fs::File;
use std::io::BufRead;
use std::io::{Error as IOError, ErrorKind};

pub struct FastQ {
    pub header:   String,
    pub sequence: Nucleotides,
    pub quality:  QualityScores,
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
    /// Reads a fastq file into an iterator backed by a buffered reader.
    ///
    /// # Errors
    ///
    /// Will return `Err` if file or permissions do not exist.
    pub fn from_filename<P>(filename: &str) -> Result<FastQReader<File>, std::io::Error>
    where
        P: AsRef<P>, {
        let file = File::open(filename)?;
        Ok(FastQReader::new(file))
    }
}

// Consider making fallible
impl<R: std::io::Read> Iterator for FastQReader<R> {
    type Item = std::io::Result<FastQ>;

    fn next(&mut self) -> Option<Self::Item> {
        self.fastq_buffer.clear();

        // Read HEADER line
        match self.fastq_reader.read_until(b'\n', &mut self.fastq_buffer) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }

        let header = if self.fastq_buffer.len() > 1 {
            if self.fastq_buffer.ends_with(b"\n") {
                self.fastq_buffer.pop();
            }

            match String::from_utf8(self.fastq_buffer.clone()) {
                Ok(s) => s,
                Err(e) => return Some(Err(IOError::new(ErrorKind::InvalidData, e))),
            }
        } else {
            String::new()
        };
        self.fastq_buffer.clear();

        // Read SEQUENCE line
        match self.fastq_reader.read_until(b'\n', &mut self.fastq_buffer) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }

        let sequence = if self.fastq_buffer.len() > 1 {
            if self.fastq_buffer.ends_with(b"\n") {
                self.fastq_buffer.pop();
            }

            Nucleotides(self.fastq_buffer.iter().copied().map(|c| c.to_ascii_uppercase()).collect())
        } else {
            Nucleotides::new()
        };
        self.fastq_buffer.clear();

        // Read "+" line
        match self.fastq_reader.read_until(b'\n', &mut self.fastq_buffer) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }

        self.fastq_buffer.clear();

        // Read QUALITY line
        match self.fastq_reader.read_until(b'\n', &mut self.fastq_buffer) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }

        let quality = if self.fastq_buffer.len() > 1 {
            if self.fastq_buffer.ends_with(b"\n") {
                self.fastq_buffer.pop();
            }

            QualityScores(self.fastq_buffer.iter().copied().map(|c| c.to_ascii_uppercase()).collect())
        } else {
            QualityScores::new()
        };

        Some(Ok(FastQ {
            header,
            sequence,
            quality,
        }))
    }
}
