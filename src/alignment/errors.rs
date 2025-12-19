use crate::data::err::GetCode;
use std::{error::Error, fmt};

/// Errors that can occur when constructing a profile for alignment.
#[derive(Debug)]
pub enum ProfileError {
    /// The provided sequence was empty.
    EmptySequence,
    /// The gap open weight was out of range (must be between -127 and 0).
    GapOpenOutOfRange { gap_open: i8 },
    /// The gap extend weight was out of range (must be between -127 and 0).
    GapExtendOutOfRange { gap_extend: i8 },
    /// The gap open weight must be less than or equal to the gap extend weight.
    BadGapWeights { gap_open: i8, gap_extend: i8 },
}

impl fmt::Display for ProfileError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ProfileError::EmptySequence => write!(f, "The provided sequence was empty"),
            ProfileError::GapOpenOutOfRange { gap_open } => {
                write!(
                    f,
                    "The gap open weight must be between -127 and 0, but {gap_open} was provided"
                )
            }
            ProfileError::GapExtendOutOfRange { gap_extend } => {
                write!(
                    f,
                    "The gap extend weight must be between -127 and 0, but {gap_extend} was provided"
                )
            }
            ProfileError::BadGapWeights { gap_open, gap_extend } => write!(
                f,
                "The gap open weight must be less than or equal to the gap extend weight, but {gap_open} (gap open) and {gap_extend} (gap extend) were provided"
            ),
        }
    }
}

impl Error for ProfileError {}
impl GetCode for ProfileError {}

/// An enum representing errors that can happen when calculating an alignment
/// score for a particular CIGAR string.
#[derive(PartialEq)]
#[non_exhaustive]
pub enum ScoringError {
    /// The CIGAR string produced a negative score
    NegativeScore(i32),
    /// Query ended before the entire CIGAR string was consumed
    QueryEnded,
    /// Reference ended before the entire CIGAR string was consumed
    ReferenceEnded,
    /// Failed to consume the full, provided query, which was expected to
    /// contain no more than what was represented by the CIGAR string
    FullQueryNotUsed,
    /// Failed to consume the full, provided reference, which was expected to
    /// contain only the aligned region of the original reference
    FullReferenceNotUsed,
    /// Unsupported CIGAR opcode used in argument
    InvalidCigarOp(u8),
}

impl fmt::Display for ScoringError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ScoringError::NegativeScore(score) => write!(f, "The alignment produced a negative score: {score}"),
            ScoringError::QueryEnded => write!(f, "The query ended before the entire CIGAR string was consumed!"),
            ScoringError::ReferenceEnded => write!(f, "The reference ended before the entire CIGAR string was consumed!"),
            ScoringError::FullQueryNotUsed => write!(
                f,
                "Failed to consume the full, provided query, which was expected to contain no more than what was represented by the CIGAR string"
            ),
            ScoringError::FullReferenceNotUsed => {
                write!(
                    f,
                    "Failed to consume the full, provided reference, which was expected to contain only the aligned region of the original reference"
                )
            }
            ScoringError::InvalidCigarOp(op) => write!(f, "An unsupported CIGAR opcode was encountered: {op}"),
        }
    }
}

impl fmt::Debug for ScoringError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for ScoringError {}
impl GetCode for ScoringError {}
