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
