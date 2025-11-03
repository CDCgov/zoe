use std::{
    fs::File,
    io::{BufRead, Error as IOError, ErrorKind},
    path::Path,
};

use crate::{
    data::{
        cigar::Cigar,
        records::RecordReader,
        sam::{SamData, SamTags},
    },
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

/// An iterator for buffered reading of a SAM file.
///
/// This iterator guarantees quality scores are valid (see [`QualityScores`] for
/// more details). It does not check for the validity of the CIGAR strings.
///
/// [`QualityScores`]: crate::data::types::phred::QualityScores
#[derive(Debug)]
pub struct SAMReader<R: std::io::Read, const TAGS: bool> {
    sam_reader: std::io::Lines<std::io::BufReader<R>>,
}

impl<R: std::io::Read> SAMReader<R, true> {
    /// Creates an iterator over SAM data, wrapping the input in a buffered
    /// reader.
    ///
    /// Unlike [`from_readable`], this does not allocate or read any data
    /// initially. It also allows for empty input, in which case the resulting
    /// iterator is empty.
    ///
    /// If the optional tags in the SAM records are not being used by the
    /// downstream application, consider using [`new_ignore_tags`] for
    /// efficiency.
    ///
    /// [`from_readable`]: SAMReader::from_readable
    /// [`new_ignore_tags`]: SAMReader::new_ignore_tags
    pub fn new(inner: R) -> Self {
        SAMReader {
            sam_reader: std::io::BufReader::new(inner).lines(),
        }
    }

    /// Creates an iterator over SAM data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// If the optional tags in the SAM records are not being used by the
    /// downstream application, consider using [`from_readable_ignore_tags`] for
    /// efficiency.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    ///
    /// [`Read`]: std::io::Read
    /// [`from_readable_ignore_tags`]: SAMReader::from_readable_ignore_tags
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

impl<R: std::io::Read> SAMReader<R, false> {
    /// Creates an iterator over SAM data, wrapping the input in a buffered
    /// reader.
    ///
    /// This version ignores any optional tags present in the SAM records in
    /// order to reduce allocations and improve efficiency. Use [`new`] to
    /// include the tags in the output.
    ///
    /// Unlike [`from_readable_ignore_tags`], this does not allocate or read any
    /// data initially. It also allows for empty input, in which case the
    /// resulting iterator is empty.
    ///
    /// [`from_readable_ignore_tags`]: SAMReader::from_readable_ignore_tags
    /// [`new`]: SAMReader::new
    pub fn new_ignore_tags(inner: R) -> Self {
        SAMReader {
            sam_reader: std::io::BufReader::new(inner).lines(),
        }
    }

    /// Creates an iterator over SAM data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// This version ignores any optional tags present in the SAM records in
    /// order to reduce allocations and improve efficiency. Use
    /// [`from_readable`] to include the tags in the output.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    ///
    /// [`Read`]: std::io::Read
    /// [`from_readable`]: SAMReader::from_readable
    pub fn from_readable_ignore_tags(read: R) -> std::io::Result<Self> {
        let mut sam_reader = std::io::BufReader::new(read);
        if sam_reader.fill_buf()?.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "No SAM data was found!"));
        }

        Ok(SAMReader {
            sam_reader: sam_reader.lines(),
        })
    }
}

impl<R: std::io::Read, const TAGS: bool> Iterator for SAMReader<R, TAGS> {
    type Item = std::io::Result<SamRow>;

    fn next(&mut self) -> Option<Self::Item> {
        let line = unwrap_or_return_some_err!(self.sam_reader.next()?);

        if line.starts_with('@') {
            Some(Ok(SamRow::Header(line)))
        } else {
            let parts = line.split('\t').collect::<Vec<_>>();

            if parts.len() < 11 {
                return Some(Err(IOError::new(
                    ErrorKind::InvalidData,
                    "SAM file did not have at least 11 fields!",
                )));
            }

            let qname = parts[0].to_string();

            let flag = {
                match parts[1].parse::<u16>() {
                    Ok(f) => f,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM bit flag '{}'", parts[1]),
                        )));
                    }
                }
            };

            let rname = parts[2].to_string();

            let pos = {
                match parts[3].parse::<usize>() {
                    Ok(p) => p,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM reference position '{}'", parts[3]),
                        )));
                    }
                }
            };

            let mapq = {
                match parts[4].parse::<u8>() {
                    Ok(q) => q,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM mapq '{}'", parts[4]),
                        )));
                    }
                }
            };

            let cigar = {
                let cigar_str = parts[5].trim_ascii().as_bytes();
                if cigar_str == b"*" {
                    Cigar::new()
                } else {
                    Cigar::from_slice_unchecked(cigar_str)
                }
            };

            // TODO: Where should casing be handled?
            let seq = parts[9].as_bytes().to_ascii_uppercase().into();

            let qual = unwrap_or_return_some_err!(parts[10].as_bytes().try_into());

            let tags = if TAGS {
                parts[11..].iter().map(ToString::to_string).collect::<SamTags>()
            } else {
                SamTags::new()
            };

            let row = SamData {
                qname,
                flag,
                rname,
                pos,
                mapq,
                cigar,
                rnext: '*',
                pnext: 0,
                tlen: 0,
                seq,
                qual,
                tags,
            };
            Some(Ok(SamRow::Data(row)))
        }
    }
    // Push down predicates?
}

impl SAMReader<File, true> {
    /// Reads a SAM text file into an iterator backed by a buffered reader.
    ///
    /// If the optional tags in the SAM records are not being used by the
    /// downstream application, consider using [`from_filename_ignore_tags`] for
    /// efficiency.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if file or permissions do not exist, or if the file is
    /// empty. The file path is included in the error message.
    ///
    /// [`from_filename_ignore_tags`]: SAMReader::from_filename_ignore_tags
    pub fn from_filename<P>(filename: P) -> Result<Self, std::io::Error>
    where
        P: AsRef<Path>, {
        let file = Self::open_nonempty_file(filename)?;
        Ok(SAMReader::new(file))
    }
}

impl SAMReader<File, false> {
    /// Reads a SAM text file into an iterator backed by a buffered reader.
    ///
    /// This version ignores any optional tags present in the SAM records in
    /// order to reduce allocations and improve efficiency. Use [`from_filename`] to
    /// include the tags in the output.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if file or permissions do not exist, or if the file is
    /// empty. The file path is included in the error message.
    ///
    /// [`from_filename`]: SAMReader::from_filename
    pub fn from_filename_ignore_tags<P>(filename: P) -> Result<Self, std::io::Error>
    where
        P: AsRef<Path>, {
        let file = Self::open_nonempty_file(filename)?;
        Ok(SAMReader::new_ignore_tags(file))
    }
}

impl<R: std::io::Read, const TAGS: bool> RecordReader for SAMReader<R, TAGS> {
    const RECORD_NAME: &str = "SAM";
}
