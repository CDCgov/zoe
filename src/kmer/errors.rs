use crate::data::err::GetCode;
use std::{error::Error, fmt};

#[non_exhaustive]
pub enum KmerError {
    InvalidLength,
}

impl fmt::Display for KmerError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            KmerError::InvalidLength => write!(f, "kmer length must be between 2 and 21"),
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
