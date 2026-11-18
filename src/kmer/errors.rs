//! K-mer related error types

use crate::data::err::impl_err_conv;
use std::{
    error::Error,
    fmt::{Debug, Display},
};

/// A k-mer related error.
#[non_exhaustive]
pub enum KmerError {
    /// The specified k-mer length was invalid (not between `2` and `MAX_LEN`)
    InvalidLength,
}

impl Display for KmerError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            KmerError::InvalidLength => write!(f, "k-mer length must be between 2 and MAX_LEN"),
        }
    }
}

impl Debug for KmerError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for KmerError {}

impl_err_conv!(KmerError);
