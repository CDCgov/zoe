//! Arbitrary implementations and specification structs for record types in the
//! [`records`](crate::data::records) module.

mod fasta;
mod fastq;
mod sam;

pub use fastq::*;
pub use sam::*;

#[cfg(feature = "dev-generic-fasta")]
mod generic_fasta;

#[cfg(feature = "dev-generic-fasta")]
pub use generic_fasta::*;
