use crate::data::err::GetCode;
use std::{error::Error, fmt};

#[non_exhaustive]
#[derive(Clone, Eq, PartialEq, Copy, Hash)]
/// Custom error type for constructing CIGAR strings.
pub enum CigarError {
    InvalidOperation,
    IncZero,
    IncOverflow,
    MissingInc,
    MissingOp,
}

impl fmt::Display for CigarError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            CigarError::InvalidOperation => write!(f, "CIGAR operator must be one of: M, I, D, N, S, H, P, X, ="),
            CigarError::IncZero => write!(f, "CIGAR increment must be a non-zero positive integer"),
            CigarError::IncOverflow => write!(f, "CIGAR increment must be smaller than {}", usize::MAX),
            CigarError::MissingInc => write!(f, "CIGAR operator must have preceding increment"),
            CigarError::MissingOp => write!(f, "CIGAR increment must be followed by operator"),
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
