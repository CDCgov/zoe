use crate::data::err::GetCode;
use std::{error::Error, fmt};

/// A k-mer related error.
#[non_exhaustive]
pub enum KmerError {
    /// The specified k-mer length was invalid (not between `2` and `MAX_LEN`)
    InvalidLength,
}

impl fmt::Display for KmerError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            KmerError::InvalidLength => write!(f, "kmer length must be between 2 and MAX_LEN"),
        }
    }
}

impl fmt::Debug for KmerError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for KmerError {}
impl GetCode for KmerError {}
