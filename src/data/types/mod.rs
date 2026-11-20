//! A module for storing fundamental bioinformatic types, like [`Nucleotides`],
//! [`AminoAcids`], and [`QualityScores`].
//!
//! [`Nucleotides`]: crate::data::types::nucleotides::Nucleotides
//! [`AminoAcids`]: crate::data::types::amino_acids::AminoAcids
//! [`QualityScores`]: crate::data::types::phred::QualityScores

pub mod amino_acids;
pub mod cigar;
pub mod nucleotides;
pub mod phred;

mod traits;

pub use traits::*;
