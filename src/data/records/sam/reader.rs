use crate::{
    data::{
        cigar::Cigar,
        err::ResultWithErrorContext,
        nucleotides::Nucleotides,
        phred::QualityScores,
        sam::{SamData, SamOptRaw, is_missing_sam_field},
    },
    unwrap_or_return_some_err,
};
use std::{
    fs::File,
    io::{BufRead, Error as IOError, ErrorKind},
    path::Path,
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
/// ## Parameters
///
/// - `R`: The type of the underlying reader
/// - `OPT`: Whether optional fields should be parsed
///
/// [`QualityScores`]: crate::data::types::phred::QualityScores
#[derive(Debug)]
pub struct SAMReader<R: std::io::Read, const OPT: bool> {
    sam_reader: std::io::Lines<std::io::BufReader<R>>,
}

impl<R: std::io::Read> SAMReader<R, true> {
    /// Creates an iterator over SAM data, wrapping the input in a buffered
    /// reader.
    ///
    /// Unlike [`from_readable`], this does not read any data initially. It also
    /// allows for empty input, in which case the resulting iterator is empty.
    ///
    /// If the optional data in the SAM records is not being used by the
    /// downstream application, consider using [`new_ignore_opt`] for
    /// efficiency.
    ///
    /// [`from_readable`]: SAMReader::from_readable
    /// [`new_ignore_opt`]: SAMReader::new_ignore_opt
    pub fn new(inner: R) -> Self {
        SAMReader {
            sam_reader: std::io::BufReader::new(inner).lines(),
        }
    }

    /// Creates an iterator over SAM data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// If the optional data in the SAM records is not being used by the
    /// downstream application, consider using [`from_readable_ignore_opt`] for
    /// efficiency.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    ///
    /// [`Read`]: std::io::Read
    /// [`from_readable_ignore_opt`]: SAMReader::from_readable_ignore_opt
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
    /// Deprecated, use [`new_ignore_opt`] instead.
    ///
    /// [`new_ignore_opt`]: SAMReader::new_ignore_opt
    #[deprecated(
        since = "0.0.29",
        note = "please use `new_ignore_opt` instead. This function will be removed in v0.0.31"
    )]
    pub fn new_ignore_tags(inner: R) -> Self {
        SAMReader::new_ignore_opt(inner)
    }

    /// Creates an iterator over SAM data, wrapping the input in a buffered
    /// reader.
    ///
    /// This version ignores any optional data present in the SAM records in
    /// order to reduce allocations and improve efficiency. Use [`new`] to
    /// include the optional data in the output.
    ///
    /// Unlike [`from_readable_ignore_opt`], this does not read any data
    /// initially. It also allows for empty input, in which case the resulting
    /// iterator is empty.
    ///
    /// [`from_readable_ignore_opt`]: SAMReader::from_readable_ignore_opt
    /// [`new`]: SAMReader::new
    pub fn new_ignore_opt(inner: R) -> Self {
        SAMReader {
            sam_reader: std::io::BufReader::new(inner).lines(),
        }
    }

    /// Deprecated, use [`from_readable_ignore_opt`] instead.
    ///
    /// [`from_readable_ignore_opt`]: SAMReader::from_readable_ignore_opt
    #[allow(clippy::missing_errors_doc)]
    #[deprecated(
        since = "0.0.29",
        note = "please use `from_readable_ignore_opt` instead. This function will be removed in v0.0.31"
    )]
    pub fn from_readable_ignore_tags(read: R) -> std::io::Result<Self> {
        SAMReader::from_readable_ignore_opt(read)
    }

    /// Creates an iterator over SAM data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// This version ignores any optional data present in the SAM records in
    /// order to reduce allocations and improve efficiency. Use
    /// [`from_readable`] to include the optional data in the output.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    ///
    /// [`Read`]: std::io::Read
    /// [`from_readable`]: SAMReader::from_readable
    pub fn from_readable_ignore_opt(read: R) -> std::io::Result<Self> {
        let mut sam_reader = std::io::BufReader::new(read);
        if sam_reader.fill_buf()?.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "No SAM data was found!"));
        }

        Ok(SAMReader {
            sam_reader: sam_reader.lines(),
        })
    }
}

impl<R: std::io::Read, const OPT: bool> Iterator for SAMReader<R, OPT> {
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

            let seq: Nucleotides = parts[9].into();

            let qual = if is_missing_sam_field(parts[10]) {
                QualityScores::new()
            } else {
                unwrap_or_return_some_err!(parts[10].as_bytes().try_into())
            };

            let opt_fields = if OPT {
                parts[11..].iter().map(ToString::to_string).collect::<SamOptRaw>()
            } else {
                SamOptRaw::new()
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
                opt_fields,
            };
            Some(Ok(SamRow::Data(row)))
        }
    }
    // Push down predicates?
}

impl SAMReader<File, true> {
    /// Creates an iterator over the SAM data contained in a path, using a
    /// buffered reader.
    ///
    /// If the optional data in the SAM records is not being used by the
    /// downstream application, consider using [`from_path_ignore_opt`] for
    /// efficiency.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the path does not exist, if there are insufficient
    /// permissions to read from it, or if it contains no data. The path is
    /// included in the error message.
    ///
    /// [`from_path_ignore_opt`]: SAMReader::from_path_ignore_opt
    pub fn from_path<P>(path: P) -> Result<Self, std::io::Error>
    where
        P: AsRef<Path>, {
        let path = path.as_ref();
        let file = File::open(path).with_path_context("Failed to open path", path)?;
        Ok(Self::from_readable(file).with_path_context("Failed to read data at path", path)?)
    }
}

impl SAMReader<File, false> {
    /// Deprecated, use [`from_path_ignore_opt`] instead.
    ///
    /// [`from_path_ignore_opt`]: SAMReader::from_path_ignore_opt
    #[allow(clippy::missing_errors_doc)]
    #[deprecated(
        since = "0.0.29",
        note = "please use `from_path_ignore_opt` instead. This function will be removed in v0.0.31"
    )]
    pub fn from_path_ignore_tags<P>(path: P) -> Result<Self, std::io::Error>
    where
        P: AsRef<Path>, {
        Self::from_path_ignore_opt(path)
    }

    /// Creates an iterator over the SAM data contained in a path, using a
    /// buffered reader.
    ///
    /// This version ignores any optional data present in the SAM records in
    /// order to reduce allocations and improve efficiency. Use [`from_path`] to
    /// include the optional data in the output.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the path does not exist, if there are insufficient
    /// permissions to read from it, or if it contains no data. The path is
    /// included in the error message.
    ///
    /// [`from_path`]: SAMReader::from_path
    pub fn from_path_ignore_opt<P>(path: P) -> Result<Self, std::io::Error>
    where
        P: AsRef<Path>, {
        let path = path.as_ref();
        let file = File::open(path).with_path_context("Failed to open path", path)?;
        Ok(Self::from_readable_ignore_opt(file).with_path_context("Failed to read data at path", path)?)
    }
}
