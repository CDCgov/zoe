use std::{
    fs::File,
    io::{BufRead, Error as IOError, ErrorKind},
    path::Path,
};

use crate::{
    data::{cigar::Cigar, sam::SamData},
    unwrap_or_return_some_err,
};

/// Enum for holding either a SAM header or data record.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub enum SamRow {
    Header(String),
    Data(SamData),
}

impl SamRow {
    /// If the row has the variant [`Header`], returns `Some` containing the
    /// line. Otherwise return `None`.
    ///
    /// [`Header`]: SamRow::Header
    #[inline]
    #[must_use]
    pub fn header(self) -> Option<String> {
        match self {
            SamRow::Header(header) => Some(header),
            SamRow::Data(_) => None,
        }
    }

    /// If the row has the variant [`Data`], returns `Some` containing the
    /// [`SamData`]. Otherwise return `None`.
    ///
    /// [`Data`]: SamRow::Data
    #[inline]
    #[must_use]
    pub fn data(self) -> Option<SamData> {
        match self {
            SamRow::Header(_) => None,
            SamRow::Data(data) => Some(data),
        }
    }
}

impl std::fmt::Display for SamRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SamRow::Data(d) => write!(f, "{d}"),
            SamRow::Header(h) => write!(f, "{h}"),
        }
    }
}

/// An iterator for buffered reading of a SAM file. Guarantees quality scores
/// are valid.
#[derive(Debug)]
pub struct SAMReader<R: std::io::Read> {
    pub sam_reader: std::io::Lines<std::io::BufReader<R>>,
}

impl<R: std::io::Read> SAMReader<R> {
    /// Creates an iterator over SAM data, wrapping the input in a buffered
    /// reader.
    ///
    /// Unlike [`from_readable`], this does not allocate or read any data
    /// initially. It also allows for empty input, in which case the resulting
    /// iterator is empty.
    ///
    /// [`from_readable`]: SAMReader::from_readable
    pub fn new(inner: R) -> Self {
        SAMReader {
            sam_reader: std::io::BufReader::new(inner).lines(),
        }
    }

    /// Creates an iterator over SAM data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    ///
    /// [`Read`]: std::io::Read
    pub fn from_readable(read: R) -> std::io::Result<Self> {
        let mut sam_reader = std::io::BufReader::new(read);
        if sam_reader.fill_buf()?.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "No SAM data was found!"));
        }

        Ok(SAMReader {
            sam_reader: sam_reader.lines(),
        })
    }
}

impl<R: std::io::Read> Iterator for SAMReader<R> {
    type Item = std::io::Result<SamRow>;

    fn next(&mut self) -> Option<Self::Item> {
        let line = self.sam_reader.next()?;

        match line {
            Ok(s) if s.starts_with('@') => Some(Ok(SamRow::Header(s))),
            Ok(s) => {
                let r: Vec<&str> = s.split('\t').collect();
                // TODO: account for optional fields
                if r.len() < 11 {
                    return Some(Err(IOError::new(
                        ErrorKind::InvalidData,
                        "SAM file did not have at least 11 fields!",
                    )));
                }

                let read_bit_flags = match r[1].parse::<u16>() {
                    Ok(f) => f,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM bit flag '{}'", &r[1]),
                        )));
                    }
                };

                // Per the SAM spec, the 1-based reference position
                let aligned_reference_position = match r[3].parse::<usize>() {
                    Ok(p) => p,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM reference position '{}'", &r[3]),
                        )));
                    }
                };

                let read_mapping_quality = match r[4].parse::<u8>() {
                    Ok(q) => q,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM mapq '{}'", &r[4]),
                        )));
                    }
                };

                let quality_scores = unwrap_or_return_some_err!(r[10].as_bytes().try_into());

                let cigar = {
                    let out = r[5].trim_ascii().as_bytes();
                    if out == b"*" {
                        Cigar::new()
                    } else {
                        Cigar::from_slice_unchecked(out)
                    }
                };

                let row = SamData {
                    qname: r[0].to_owned(),
                    flag: read_bit_flags,
                    rname: r[2].to_owned(),
                    pos: aligned_reference_position,
                    mapq: read_mapping_quality,
                    cigar,
                    rnext: '*',
                    pnext: 0,
                    tlen: 0,
                    // TODO: Where should casing be handled?
                    seq: r[9].as_bytes().to_ascii_uppercase().into(),
                    qual: quality_scores,
                };
                Some(Ok(SamRow::Data(row)))
            }
            Err(e) => Some(Err(e)),
        }
    }
    // Push down predicates?
}

impl SAMReader<std::fs::File> {
    /// Reads a SAM text file into an iterator backed by a buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if file or permissions do not exist, or if the file is
    /// empty.
    pub fn from_filename<P>(filename: P) -> Result<SAMReader<File>, std::io::Error>
    where
        P: AsRef<Path>, {
        let file = File::open(filename)?;
        if file.metadata()?.len() == 0 {
            return Err(IOError::new(ErrorKind::InvalidData, "SAM file was empty!"));
        }
        Ok(SAMReader::new(file))
    }
}
