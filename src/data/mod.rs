/// A module containing traits for convertion our custom types, such as
/// [`Vec<u8>`] to [`Nucleotides`](self::types::nucleotides::Nucleotides).
pub mod convert;
/// A module with error types and convenience traits for handling [`Result`].
pub mod err;

/// A module for reading and manipulating
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.
pub mod fasta;
/// A module for reading and manipulating
/// [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files.
pub mod fastq;
/// A module for reading and manipulating
/// [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files. Provides some
/// special-case functions used by [IRMA](https://wonder.cdc.gov/amd/flu/irma/).
pub mod sam;
/// A module for storing more fundamental types, like
/// [`Nucleotides`](self::types::nucleotides::Nucleotides) and
/// [`AminoAcids`][self::types::amino_acids::AminoAcids].
pub mod types;

pub use crate::{
    data::vec_types::{CheckSequence, PairwiseSequence, StdForSequences},
    search::ByteSubstring,
};
pub use constants::{
    mappings::{ResidueMapping, DNA_RESIDUE_MAPPING},
    matrices::{BiasedWeightMatrix, SimpleWeightMatrix},
};

/// A private module for helper functions on [`u8`].
pub(crate) mod byte_types;
/// A private module for helper alphabets, maps, and matrices that can be used
/// within public methods.
pub(crate) mod constants;
/// A private module for working with IDs.
pub(crate) mod id_types;
/// A private module with helper traits for things like [`Vec<u8>`].
pub(crate) mod vec_types;

pub(crate) use constants::*;
