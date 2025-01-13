/// A module with error types and convenience traits for handling [`Result`].
pub mod err;
/// A module for records types--usually for I/O--that are structures of other
/// more primitive types.
pub mod records;
/// A module for storing more fundamental types, like
/// [`Nucleotides`](self::types::nucleotides::Nucleotides) and
/// [`AminoAcids`][self::types::amino_acids::AminoAcids].
pub mod types;
/// A module containing the traits needed for working with views.
pub mod view_traits;

/// A private module for helper alphabets, maps, and matrices that can be used
/// within public methods.
pub(crate) mod constants;
/// A private module for helper extension traits.
pub(crate) mod extension;

/// Used for type validation
mod validation;

pub use constants::{
    mappings::{ByteIndexMap, DNA_PROFILE_MAP},
    matrices::{BiasedWeightMatrix, SimpleWeightMatrix},
};
pub use records::{fasta, fastq, sam};
pub use types::{amino_acids, cigar, nucleotides, phred};
pub use validation::{CheckSequence, Recode, RetainSequence, StdForSequences};

pub(crate) use constants::{alphas, mappings, matrices};
pub(crate) use extension::{array_types, byte_types, id_types, vec_types};
