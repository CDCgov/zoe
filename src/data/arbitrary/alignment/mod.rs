//! Arbitrary implementations and specification structs for alignment-related
//! structs.

mod output;
#[cfg(feature = "dev-phmm")]
mod phmm;

pub use output::*;
#[cfg(feature = "dev-phmm")]
pub use phmm::*;
