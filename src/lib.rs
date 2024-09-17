#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]
#![warn(clippy::all, clippy::pedantic)]
#![allow(
    clippy::module_name_repetitions,
    clippy::similar_names,
    clippy::wildcard_imports,
    clippy::enum_glob_use,
    stable_features
)]
#![feature(test, portable_simd, const_fn_floating_point_arithmetic, let_chains, lazy_cell)]
#![forbid(incomplete_features)]

/// Alignment functions
pub mod alignment;
/// Composition and consensus functions.
pub mod composition;
/// Data import, export, and manipulation functions.
pub mod data;
/// Distance functions, especially for sequence data.
pub mod distance;
/// Sequence search and/or replacement.
pub mod search;

/// Generate sequences and other data.
pub(crate) mod generate;
/// Iterator utilities.
pub(crate) mod iter_utils;
/// Mathematical utilities.
pub(crate) mod math;
/// SIMD traits to extend portable SIMD.
pub(crate) mod simd;
/// Private sorting functions for ASCII use cases.
pub(crate) mod sort;

/// Common structures and traits re-exported
pub mod prelude {
    pub use crate::composition::{AlignmentComposition, CreateConsensus, NucleotideCounts};
    pub use crate::data::types::{amino_acids::AminoAcids, nucleotides::Nucleotides};
    pub use crate::data::vec_types::PairwiseSequence;
    pub use crate::data::{convert::ToDNA, err::OrFail, fasta::FastaReader, fastq::FastQReader};
    pub use crate::generate::rand_sequence;
    pub use crate::search::{ByteSubstring, ByteSubstringMut};
}
