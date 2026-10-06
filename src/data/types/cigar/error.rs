//! Error types for CIGAR strings.

use std::{
    error::Error,
    fmt::{self, Display},
};

// TODO: Could eventually have InvalidOperation hold InvalidCigarOpErr
#[non_exhaustive]
#[derive(Clone, Eq, PartialEq, Copy, Hash)]
/// Custom error type for constructing CIGAR strings.
pub enum CigarError {
    /// The CIGAR operation must be one of: `M, I, D, N, S, H, P, X, =`.
    InvalidOperation,
    /// The CIGAR increment must be a non-zero positive integer.
    IncZero,
    /// The CIGAR increment must be less than or equal to [`usize::MAX`].
    IncOverflow,
    /// The CIGAR increment after merging consecutive ciglets with the same
    /// operation must be less than or equal to [`usize::MAX`].
    MergeIncOverflow,
    /// The CIGAR operation must have preceding increment.
    MissingInc,
    /// The CIGAR increment must be followed by operation.
    MissingOp,
}

impl fmt::Display for CigarError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            CigarError::InvalidOperation => f.write_str("CIGAR operation must be one of: M, I, D, N, S, H, P, X, ="),
            CigarError::IncZero => f.write_str("CIGAR increment must be a non-zero positive integer"),
            CigarError::IncOverflow => write!(f, "CIGAR increment must be less than or equal to {}", usize::MAX),
            CigarError::MergeIncOverflow => write!(
                f,
                "After merging consecutive ciglets with the same operation, the total increment must be less than or equal to {}",
                usize::MAX
            ),
            CigarError::MissingInc => f.write_str("CIGAR operation must have preceding increment"),
            CigarError::MissingOp => f.write_str("CIGAR increment must be followed by operation"),
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

/// An error for an invalid operation while parsing a CIGAR string.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[repr(transparent)]
pub struct InvalidCigarOpErr(pub u8);

impl Display for InvalidCigarOpErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "An invalid CIGAR operation was found: '{}' (expected one of M, I, D, N, S, H, P, X, =)",
            self.0 as char
        )
    }
}

impl Error for InvalidCigarOpErr {}

impl From<InvalidCigarOpErr> for CigarError {
    fn from(_value: InvalidCigarOpErr) -> Self {
        CigarError::InvalidOperation
    }
}
