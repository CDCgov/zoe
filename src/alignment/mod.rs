//! ## Functions for aligning sequence data.
//!
//! *Zoe* supports alignment for DNA, protein, or any other
//! sequence data.
//!
//! - [Smith-Waterman]: Optimal local alignment in the [`sw`] module.
//!
//! [Smith-Waterman]: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

#[cfg(feature = "dev-phmm")]
pub mod phmm;
pub mod sw;

mod errors;
mod pairwise;
mod profile;
mod profile_set;
mod std_traits;
mod types;

pub use errors::*;
pub use pairwise::*;
pub use profile::*;
pub use profile_set::*;
pub use types::*;
