pub mod fasta;
pub mod fastq;
pub mod sam;

pub mod types;

pub(crate) mod byte_types;
pub(crate) mod matrices;
pub(crate) mod vec_types;

pub mod err;

#[allow(unused_imports)]
pub use vec_types::{BiologicalSequence, ContainsSubsequence};
