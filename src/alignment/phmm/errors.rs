use crate::data::err::GetCode;
use std::{error::Error, fmt};

/// An enum representing errors that can happen when working with pHMMs.
#[non_exhaustive]
#[derive(Eq, PartialEq)]
pub enum PhmmError {
    /// The model has no layers
    EmptyModel,
    /// No alignment with nonzero probability can be found with the pHMM
    NoAlignmentFound,
    /// An invalid path through a pHMM was specified
    InvalidPath,
    /// Failed to consume the full, provided sequence, which was expected to
    /// contain no more than what was represented by the CIGAR string
    FullSeqNotUsed,
    /// Failed to consume the full range of the model specified, which was
    /// expected to contain only the aligned layers of the original model
    FullModelNotUsed,
    /// Invalid CIGAR string
    InvalidCigar,
    /// Unsupported CIGAR opcode used in argument
    InvalidCigarOp,
    /// Alignment modules were incompatible with core pHMM
    IncompatibleModule,
}

impl fmt::Display for PhmmError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            PhmmError::EmptyModel => write!(f, "No layers were present in the pHMM!"),
            PhmmError::NoAlignmentFound => write!(f, "No alignment with nonzero probability found!"),
            PhmmError::InvalidPath => write!(f, "The path specified is not valid for the given pHMM!"),
            PhmmError::FullSeqNotUsed => write!(
                f,
                "Failed to consume the full, provided sequence, which was expected to contain no more than what was represented by the CIGAR string"
            ),
            PhmmError::FullModelNotUsed => write!(
                f,
                "Failed to consume the full range of the model specified, which was expected to contain only the aligned layers of the original model"
            ),
            PhmmError::InvalidCigar => write!(f, "The CIGAR string was invalid!"),
            PhmmError::InvalidCigarOp => write!(f, "An unsupported CIGAR opcode was encountered!"),
            PhmmError::IncompatibleModule => write!(f, "The pHMM's alignment modules are incompatible with the core pHMM!"),
        }
    }
}

impl fmt::Debug for PhmmError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for PhmmError {}
impl GetCode for PhmmError {}
