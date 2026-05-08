use crate::data::err::GetCode;
use std::{error::Error, fmt};

/// An enum representing errors that can happen when working with pHMMs.
#[non_exhaustive]
#[derive(Eq, PartialEq)]
pub enum PhmmError {
    /// An error caused by an invalid model
    InvalidModel(InvalidModelError),
    /// No alignment with nonzero probability can be found with the pHMM
    NoAlignmentFound,
    /// An invalid path through a pHMM was specified
    InvalidPath,
    /// Invalid CIGAR string
    InvalidCigar,
    ///  Unsupported CIGAR opcode used in argument
    InvalidCigarOp(u8),
}

impl fmt::Display for PhmmError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            PhmmError::InvalidModel(_) => write!(f, "The pHMM model is invalid!"),
            PhmmError::NoAlignmentFound => write!(f, "No alignment with nonzero probability found!"),
            PhmmError::InvalidPath => write!(f, "The path specified is not valid for the given pHMM!"),
            PhmmError::InvalidCigar => write!(f, "The CIGAR string was invalid!"),
            PhmmError::InvalidCigarOp(op) => {
                write!(f, "An unsupported CIGAR operation was encountered: {op}", op = *op as char)
            }
        }
    }
}

impl fmt::Debug for PhmmError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for PhmmError {
    #[inline]
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            PhmmError::InvalidModel(invalid_model_error) => Some(invalid_model_error),
            _ => None,
        }
    }
}

impl GetCode for PhmmError {}

impl From<InvalidModelError> for PhmmError {
    #[inline]
    fn from(value: InvalidModelError) -> Self {
        PhmmError::InvalidModel(value)
    }
}

/// An enum representing errors that can happen when working with pHMMs.
#[non_exhaustive]
#[derive(Copy, Clone, Eq, PartialEq)]
pub enum InvalidModelError {
    /// The model has no layers
    EmptyModel,
    /// Not enough layers were specified
    TooFewLayers(usize),
    /// Alignment modules were incompatible with core pHMM
    IncompatibleModule,
}

impl fmt::Display for InvalidModelError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvalidModelError::EmptyModel => write!(f, "No layers were present in the pHMM!"),
            InvalidModelError::TooFewLayers(num) => {
                write!(f, "Too few layers were specified for the pHMM! At least {num} are required")
            }
            InvalidModelError::IncompatibleModule => {
                write!(f, "The pHMM's alignment modules are incompatible with the core pHMM!")
            }
        }
    }
}

impl fmt::Debug for InvalidModelError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for InvalidModelError {}
impl GetCode for InvalidModelError {}
