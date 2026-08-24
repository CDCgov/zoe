//! Error types returned by BAM serialization.
//!
//! [`BamError`] is the top-level error returned by [`BamWriter`]. It separates
//! lower-level I/O failures from BAM header, writer-state, and record-encoding
//! failures.
//!
//! Header and record failures expose their failure kind through
//! [`BamHeaderError`] and [`BamRecordError`]. Record failures also store record
//! context at the top level. Both domains can wrap [`BamEncodingError`] for
//! shared low-level encoding failures such as size overflows and embedded NUL
//! bytes.
//!
//! Variants that wrap lower-level failures, such as [`BamError::Io`],
//! [`BamEncodingError::Other`], and [`BamRecordError::InvalidCigar`], expose
//! those failures through [`Error::source`].
//!
//! [`BamWriter`]: crate::data::records::bam::writer::BamWriter

use crate::data::{
    cigar::CigarError,
    err::{ErrorWithContext, GetCode},
};
use std::{error::Error, fmt};

/// Error returned when writing BAM output or converting SAM records to BAM.
#[derive(Debug)]
#[non_exhaustive]
pub enum BamError {
    /// A filesystem or output-stream operation failed.
    Io {
        /// The underlying I/O error.
        source: std::io::Error,
    },
    /// A SAM/BAM header failed validation or serialization.
    Header {
        /// The header-specific error.
        source: BamHeaderError,
    },
    /// A SAM header line was added after the BAM header had already been
    /// written.
    HeaderAlreadyWritten,
    /// A write was attempted after the BAM writer had been finalized.
    WriterFinalized,
    /// A SAM alignment record could not be represented as BAM.
    Record {
        /// Record name context.
        qname:  String,
        /// The record-specific error..
        source: BamRecordError,
    },
}

impl BamError {
    /// Wraps a header-domain error as a top-level BAM error.
    pub fn header(source: impl Into<BamHeaderError>) -> Self {
        BamError::Header { source: source.into() }
    }

    /// Wraps a record-domain error with the record name that failed.
    pub fn record(qname: impl Into<String>, source: impl Into<BamRecordError>) -> Self {
        BamError::Record {
            qname:  qname.into(),
            source: source.into(),
        }
    }
}

/// Header-specific BAM validation failures.
#[derive(Debug)]
#[non_exhaustive]
pub enum BamHeaderError {
    /// The SAM header contains two `@SQ` lines with the same `SN` value.
    DuplicateReference {
        /// The duplicated reference name.
        name: String,
    },
    /// A lower-level BAM encoding failure occurred while processing the header.
    Encoding { source: BamEncodingError },
}

/// Record-specific BAM validation failures.
#[derive(Debug)]
#[non_exhaustive]
pub enum BamRecordError {
    /// A record references a sequence name that is absent from the `@SQ` header
    /// dictionary.
    ReferenceNotFound {
        /// The missing reference name.
        name: String,
    },
    /// A CIGAR string failed parsing or BAM-specific CIGAR validation.
    InvalidCigar { source: CigarError },
    /// Alignment coordinates exceed the range supported by BAI binning (`[0,
    /// 2^29)`).
    BinningOutOfRange,
    /// A lower-level BAM encoding failure occurred while processing the record.
    Encoding { source: BamEncodingError },
}

/// Shared BAM encoding failures that can occur in multiple domains.
#[derive(Debug)]
#[non_exhaustive]
pub enum BamEncodingError {
    /// A value does not fit in the integer type or BAM bit width used for the
    /// encoded field.
    SizeOverflow {
        /// The name of the field being encoded.
        field:  &'static str,
        /// Target type or BAM field size that would be exceeded.
        target: &'static str,
    },
    /// A BAM validation or encoding failure without a more specific public
    /// variant.
    Other {
        /// The error message.
        message: String,
        /// An optional underlying error.
        source:  Option<Box<dyn Error + Send + Sync + 'static>>,
    },
}

impl BamEncodingError {
    /// Constructs a generic BAM encoding error without an underlying source.
    pub fn other(message: impl Into<String>) -> Self {
        BamEncodingError::Other {
            message: message.into(),
            source:  None,
        }
    }

    /// Constructs a generic BAM encoding error with an underlying source.
    pub fn other_with_source(message: impl Into<String>, source: impl Error + Send + Sync + 'static) -> Self {
        BamEncodingError::Other {
            message: message.into(),
            source:  Some(Box::new(source)),
        }
    }
}

impl From<std::io::Error> for BamError {
    fn from(source: std::io::Error) -> Self {
        BamError::Io { source }
    }
}

impl From<ErrorWithContext> for BamError {
    fn from(source: ErrorWithContext) -> Self {
        BamError::Io { source: source.into() }
    }
}

impl From<BamHeaderError> for BamError {
    fn from(source: BamHeaderError) -> Self {
        BamError::header(source)
    }
}

impl From<BamEncodingError> for BamHeaderError {
    fn from(source: BamEncodingError) -> Self {
        BamHeaderError::Encoding { source }
    }
}

impl From<BamEncodingError> for BamRecordError {
    fn from(source: BamEncodingError) -> Self {
        BamRecordError::Encoding { source }
    }
}

impl fmt::Display for BamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BamError::Io { .. } => f.write_str("BAM IO error"),
            BamError::Header { .. } => f.write_str("BAM header error"),
            BamError::HeaderAlreadyWritten => {
                f.write_str("Cannot add BAM header lines after the BAM header has been written")
            }
            BamError::WriterFinalized => f.write_str("Cannot write to BAM stream after finalization"),
            BamError::Record { qname, .. } => write!(f, "BAM record error while processing {qname}"),
        }
    }
}

impl fmt::Display for BamHeaderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BamHeaderError::DuplicateReference { name } => write!(f, "Duplicate @SQ SN value {name}"),
            BamHeaderError::Encoding { .. } => f.write_str("BAM header value cannot be encoded"),
        }
    }
}

impl fmt::Display for BamRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BamRecordError::ReferenceNotFound { name } => write!(f, "Reference {name} not found in @SQ header"),
            BamRecordError::InvalidCigar { .. } => f.write_str("CIGAR string cannot be encoded as BAM"),
            BamRecordError::BinningOutOfRange => {
                f.write_str("Alignment coordinates exceed the range supported by BAI binning (`[0, 2^29)`).")
            }
            BamRecordError::Encoding { .. } => f.write_str("BAM record value cannot be encoded"),
        }
    }
}

impl fmt::Display for BamEncodingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BamEncodingError::SizeOverflow { field, target } => write!(f, "{field} does not fit into {target}"),
            BamEncodingError::Other { message, .. } => f.write_str(message),
        }
    }
}

impl Error for BamError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            BamError::Io { source } => Some(source),
            BamError::Header { source, .. } => Some(source),
            BamError::Record { source, .. } => Some(source),
            BamError::HeaderAlreadyWritten | BamError::WriterFinalized => None,
        }
    }
}

impl Error for BamHeaderError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            BamHeaderError::Encoding { source } => Some(source),
            BamHeaderError::DuplicateReference { .. } => None,
        }
    }
}

impl Error for BamRecordError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            BamRecordError::InvalidCigar { source } => Some(source),
            BamRecordError::Encoding { source } => Some(source),
            BamRecordError::ReferenceNotFound { .. } | BamRecordError::BinningOutOfRange => None,
        }
    }
}

impl Error for BamEncodingError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            BamEncodingError::Other { source, .. } => {
                source.as_ref().map(|source| source.as_ref() as &(dyn Error + 'static))
            }
            BamEncodingError::SizeOverflow { .. } => None,
        }
    }
}

impl GetCode for BamError {
    fn get_code(&self) -> i32 {
        match self {
            BamError::Io { source } => source.get_code(),
            BamError::Header { source } => source.get_code(),
            BamError::HeaderAlreadyWritten | BamError::WriterFinalized => 1,
            BamError::Record { source, .. } => source.get_code(),
        }
    }
}

impl GetCode for BamHeaderError {
    fn get_code(&self) -> i32 {
        match self {
            BamHeaderError::Encoding { source } => source.get_code(),
            BamHeaderError::DuplicateReference { .. } => 1,
        }
    }
}

impl GetCode for BamRecordError {
    fn get_code(&self) -> i32 {
        match self {
            BamRecordError::InvalidCigar { source } => source.get_code(),
            BamRecordError::Encoding { source } => source.get_code(),
            BamRecordError::ReferenceNotFound { .. } | BamRecordError::BinningOutOfRange => 1,
        }
    }
}

impl GetCode for BamEncodingError {
    fn get_code(&self) -> i32 {
        match self {
            BamEncodingError::SizeOverflow { .. } | BamEncodingError::Other { .. } => 1,
        }
    }
}
