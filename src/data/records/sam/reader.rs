//! Code for reading SAM files.

use crate::{
    data::{
        err::ResultWithErrorContext,
        sam::{ParseSamFields, SamData},
    },
    unwrap_or_return_some_err,
};
use std::{
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, Error as IOError, ErrorKind, Read},
    marker::PhantomData,
    path::Path,
};

/// A type alias to allow for easier migration.
pub type SAMReader<R> = SamReader<R>;

/// Enum for holding either a SAM header or data record.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub enum SamRow<T = SamData> {
    Header(String),
    Data(T),
}

impl<T> SamRow<T> {
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
    /// record. Otherwise return `None`.
    ///
    /// [`Data`]: SamRow::Data
    #[inline]
    #[must_use]
    pub fn data(self) -> Option<T> {
        match self {
            SamRow::Header(_) => None,
            SamRow::Data(data) => Some(data),
        }
    }
}

impl Display for SamRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SamRow::Data(d) => write!(f, "{d}"),
            SamRow::Header(h) => write!(f, "{h}"),
        }
    }
}

/// An iterator for buffered reading of a SAM file.
///
/// ## Parameters
///
/// - `R`: The type of the underlying reader
/// - `T`: The record type being parsed into, by default [`SamData`]
#[derive(Debug)]
pub struct SamReader<R, T = SamData>
where
    R: Read,
    T: ParseSamFields, {
    sam_reader: std::io::Lines<std::io::BufReader<R>>,
    record:     PhantomData<T>,
}

impl<R: Read> SamReader<R> {
    /// Creates an iterator over SAM data, wrapping the input in a buffered
    /// reader.
    ///
    /// Unlike [`from_readable`], this does not read any data initially. It also
    /// allows for empty input, in which case the resulting iterator is empty.
    ///
    /// [`from_readable`]: SAMReader::from_readable
    pub fn new(inner: R) -> Self {
        SamReader {
            sam_reader: std::io::BufReader::new(inner).lines(),
            record:     PhantomData,
        }
    }

    /// Creates an iterator over SAM data from a type implementing [`Read`],
    /// wrapping the input in a buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_readable(read: R) -> std::io::Result<Self> {
        SamReader::from_bufreader(BufReader::new(read))
    }

    /// Creates an iterator over SAM data from a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the input data is empty or an IO error occurs.
    pub fn from_bufreader(mut reader: BufReader<R>) -> std::io::Result<Self> {
        if reader.fill_buf()?.is_empty() {
            return Err(IOError::new(ErrorKind::InvalidData, "No FASTQ data was found!"));
        }

        Ok(SamReader {
            sam_reader: reader.lines(),
            record:     PhantomData,
        })
    }

    /// Changes the yielded record type to be `T` instead of [`SamData`].
    pub fn with_record_type<T: ParseSamFields>(self) -> SamReader<R, T> {
        SamReader {
            sam_reader: self.sam_reader,
            record:     PhantomData,
        }
    }
}

impl<R: Read, T: ParseSamFields> Iterator for SamReader<R, T> {
    type Item = std::io::Result<SamRow>;

    fn next(&mut self) -> Option<Self::Item> {
        let line = unwrap_or_return_some_err!(self.sam_reader.next()?);

        if line.starts_with('@') {
            Some(Ok(SamRow::Header(line)))
        } else {
            Some(SamData::parse(&line).map(SamRow::Data))
        }
    }
}

impl SamReader<File> {
    /// Creates an iterator over the SAM data contained in a path, using a
    /// buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if the path does not exist, if there are insufficient
    /// permissions to read from it, or if it contains no data. The path is
    /// included in the error message.
    pub fn from_path<P>(path: P) -> Result<Self, std::io::Error>
    where
        P: AsRef<Path>, {
        let path = path.as_ref();
        let file = File::open(path).with_path_context("Failed to open path", path)?;
        Ok(Self::from_readable(file).with_path_context("Failed to read data at path", path)?)
    }
}
