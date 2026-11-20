//! A module for records types that are structures of other more primitive
//! types. Examples include [`Fasta`], [`FastQ`], and [`SamData`].
//!
//! [`Fasta`]: crate::data::records::fasta::Fasta
//! [`FastQ`]: crate::data::records::fastq::FastQ
//! [`SamData`]: crate::data::records::sam::SamData

pub mod bam;
pub mod fasta;
pub mod fastq;
pub mod sam;

mod getter_traits;
mod reader_utils;

pub use getter_traits::*;
