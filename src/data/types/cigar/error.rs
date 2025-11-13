use crate::data::err::GetCode;
use std::{error::Error, fmt};

#[non_exhaustive]
#[derive(Clone, Eq, PartialEq, Copy, Hash)]
/// Custom error type for constructing CIGAR strings.
pub enum CigarError {
    /// The CIGAR operation must be one of: `M, I, D, N, S, H, P, X, =`
    InvalidOperation,
    /// The CIGAR increment must be a non-zero positive integer
    IncZero,
    /// The CIGAR increment must be smaller than [`usize::MAX`]
    IncOverflow,
    /// The CIGAR operation must have preceding increment
    MissingInc,
    /// The CIGAR increment must be followed by operation
    MissingOp,
}

impl fmt::Display for CigarError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            CigarError::InvalidOperation => f.write_str("CIGAR operation must be one of: M, I, D, N, S, H, P, X, ="),
            CigarError::IncZero => f.write_str("CIGAR increment must be a non-zero positive integer"),
            CigarError::IncOverflow => write!(f, "CIGAR increment must be smaller than {}", usize::MAX),
            CigarError::MissingInc => f.write_str("CIGAR operation must have preceding increment"),
            CigarError::MissingOp => f.write_str("CIGAR increment must be followed by operator"),
        }
    }
}

impl fmt::Debug for CigarError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for CigarError {}
impl GetCode for CigarError {}
