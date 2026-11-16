pub mod bam;
/// A module for reading and manipulating
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.
pub mod fasta;
/// A module for reading and manipulating
/// [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files.
pub mod fastq;
pub mod sam;

mod getter_traits;
mod reader_utils;

pub use getter_traits::*;
