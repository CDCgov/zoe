pub mod err;
pub mod fasta;
pub mod fastq;
pub mod sam;
pub mod types;

pub use vec_types::{BiologicalSequence, ContainsSubsequence};

pub(crate) mod byte_types;
pub(crate) mod constants;
pub(crate) mod vec_types;

pub(crate) use constants::*;
