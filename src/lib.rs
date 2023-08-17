#![warn(clippy::all, clippy::pedantic)]
#![allow(
    clippy::module_name_repetitions,
    clippy::similar_names,
    clippy::wildcard_imports,
    clippy::enum_glob_use
)]
#![feature(test, portable_simd, const_fn_floating_point_arithmetic, iter_collect_into)]

/// Composition and consensus functions.
pub mod composition;
/// Data import, export, and manipulation functions.
pub mod data;
/// Distance functions, especially for sequence data.
pub mod distance;

/// Generate sequences and other data.
pub(crate) mod generate;
/// Mathematical utilities.
pub(crate) mod math;
/// Sequence search and/or replacement.
pub(crate) mod search;
/// SIMD traits to extend portable SIMD.
pub(crate) mod simd;
/// Sorting functions.
pub(crate) mod sort;

pub mod prelude {
    pub use crate::composition::{AlignmentComposition, CreateConsensus, NucleotideCounts};
    pub use crate::data::{convert::ToDNA, err::OrFail, fasta::FastaReader};
    pub use crate::generate::rand_sequence;
    pub use crate::search::{Subsequence, VectorSubsequence};
}
