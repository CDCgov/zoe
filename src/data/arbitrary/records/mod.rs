//! Arbitrary implementations and specification structs for record types in the
//! [`records`](crate::data::records) module.

mod fasta;
mod fastq;
mod sam;

pub use fasta::*;
pub use fastq::*;
pub use sam::*;
