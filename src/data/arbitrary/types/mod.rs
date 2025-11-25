//! Arbitrary implementations and specification structs for *Zoe* data types in
//! the [`types`](crate::data::types).

mod amino_acids;
mod cigar;
mod nucleotides;
mod phred;

pub use amino_acids::*;
pub use cigar::*;
pub use nucleotides::*;
