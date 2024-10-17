use crate::data::{
    types::{nucleotides::Nucleotides, phred::QualityScores},
    vec_types::ChopLineBreak,
};

use std::{
    fs::File,
    io::{BufRead, Error as IOError, ErrorKind},
    path::Path,
};

/// Holds data for a single [Read](https://en.wikipedia.org/wiki/Read_(biology))
/// or [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) record.
#[derive(Debug)]
pub struct FastQ {
    pub header:   String,
    pub sequence: Nucleotides,
    pub quality:  QualityScores,
}

impl std::fmt::Display for FastQ {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}\n{}\n+\n{}\n", self.header, self.sequence, self.quality)
    }
}

impl FastQ {
    /// Reverse complements the stored DNA sequence data.
    #[inline]
    pub fn reverse_complement(&mut self) {
        self.sequence = self.sequence.reverse_complement();
        self.quality = self.quality.reverse();
    }

    /// Recodes the stored sequence to an uppercase canonical (ACTG + N) one.
    /// Any non-canonical base becomes N.
    #[inline]
    pub fn recode_any_to_actgn_uc(&mut self) {
        self.sequence.recode_any_to_actgn_uc();
    }

    /// Recodes the stored sequence of valid IUPAC codes to a canonical (ACTG + N)
    /// sequence. Ambiguous bases become N while non-IUPAC bytes are left
    /// unchanged.
    #[inline]
    pub fn recode_iupac_to_actgn(&mut self) {
        self.sequence.recode_iupac_to_actgn();
    }

    /// Recodes the stored sequence of valid IUPAC codes to an uppercase
    /// canonical (ACTG + N) sequence. Ambiguous bases become N while non-IUPAC
    /// bytes are left unchanged.
    #[inline]
    pub fn recode_iupac_to_actg_uc(&mut self) {
        self.sequence.recode_iupac_to_actgn_uc();
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
    pub fn from_filename<P>(filename: P) -> Result<FastQReader<File>, std::io::Error>
    where
        P: AsRef<Path>, {
        let file = File::open(filename)?;
        Ok(FastQReader::new(file))
    }
}

/// An iterator for buffered reading of a
/// [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file. Guarantees quality
/// scores are valid.
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

        if !self.fastq_buffer.starts_with(b"@") {
            return Some(Err(IOError::new(
                ErrorKind::InvalidData,
                "Missing '@' symbol at header line beginning.",
            )));
        }

        self.fastq_buffer.chop_line_break();

        // Since '@' is in the header, header must be > 1 in length
        if self.fastq_buffer.len() <= 1 {
            return Some(Err(IOError::new(ErrorKind::InvalidData, "Missing FastQ header!")));
        }

        let header = match String::from_utf8(self.fastq_buffer.clone()) {
            Ok(s) => s,
            Err(e) => return Some(Err(IOError::new(ErrorKind::InvalidData, e))),
        };

        self.fastq_buffer.clear();

        // Read SEQUENCE line
        if let Err(e) = self.fastq_reader.read_until(b'\n', &mut self.fastq_buffer) {
            return Some(Err(e));
        }

        self.fastq_buffer.chop_line_break();

        if self.fastq_buffer.is_empty() {
            return Some(Err(IOError::new(ErrorKind::InvalidData, "Missing FastQ sequence!")));
        }

        let sequence = Nucleotides(self.fastq_buffer.clone());

        self.fastq_buffer.clear();

        // Read "+" line
        if let Err(e) = self.fastq_reader.read_until(b'\n', &mut self.fastq_buffer) {
            return Some(Err(e));
        }

        if !self.fastq_buffer.starts_with(b"+") {
            return Some(Err(IOError::new(ErrorKind::InvalidData, "Missing '+' line!")));
        }

        self.fastq_buffer.clear();

        // Read QUALITY line
        if let Err(e) = self.fastq_reader.read_until(b'\n', &mut self.fastq_buffer) {
            return Some(Err(e));
        }

        self.fastq_buffer.chop_line_break();

        if self.fastq_buffer.len() != sequence.len() {
            if self.fastq_buffer.is_empty() {
                return Some(Err(IOError::new(ErrorKind::InvalidData, "Missing FastQ quality scores!")));
            }

            return Some(Err(IOError::new(
                ErrorKind::InvalidData,
                "Sequence and quality score length mismatch!",
            )));
        }

        let quality = match QualityScores::try_from(self.fastq_buffer.as_slice()) {
            Ok(s) => s,
            Err(e) => return Some(Err(IOError::new(ErrorKind::InvalidData, e))),
        };

        Some(Ok(FastQ {
            header,
            sequence,
            quality,
        }))
    }
}

#[cfg(test)]
mod test;
