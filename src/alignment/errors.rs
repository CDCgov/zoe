use crate::data::err::GetCode;
use std::{error::Error, fmt};

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
    /// Failed to consume the full query
    FullQueryNotUsed,
    /// Failed to consume the full reference
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
            ScoringError::FullQueryNotUsed => write!(f, "Failed to consume the full query!"),
            ScoringError::FullReferenceNotUsed => {
                write!(f, "Failed to consume the full portion of the reference specified!")
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
