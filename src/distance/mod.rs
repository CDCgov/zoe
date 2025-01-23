/// Amino acid distance functions.
///
/// ## Examples:
///
/// Distance functions are provided as both standalone functions and methods via
/// the [`AminoAcidsDistance`] trait. Using standalone functions:
/// ```
/// # use zoe::{prelude::AminoAcids, distance::aa::physiochemical};
///
/// let s1 = b"MANATEE";
/// let s2 = b"MANGAEE";
///
/// let physicochemical_distance = physiochemical(s1, s2).unwrap();
///
/// assert!((physicochemical_distance  - 1.1464952).abs() < 0.01)
/// ```
///
/// Equivalently, using methods:
/// ```
/// # use zoe::{prelude::AminoAcids, distance::aa::AminoAcidsDistance};
///
/// let s1: AminoAcids = b"MANATEE".into();
/// let s2: AminoAcids = b"MANGAEE".into();
///
/// let physicochemical_distance = s1.distance_physiochemical(&s2).unwrap();
///
/// assert!((physicochemical_distance  - 1.1464952).abs() < 0.01)
/// ```
///
/// [`AminoAcidsDistance`]: aa::AminoAcidsDistance
pub mod aa;

/// Various nucleotide substitution models for calculating evolutionary
/// distances between two aligned DNA sequences.
///
/// ## Assumptions:
///
/// * __Alignment:__ both sequences must be aligned.
/// * __Pairwise deletion:__ If a gap or ambiguous (non-ACGT) base is present in
///     any position of the aligned sequences, the bases in that position are
///     excluded from the analysis for both sequences.
/// * __Unequal sequence lengths:__ If sequences are of different lengths, bases
///     in the longer sequence that extend beyond the length of the shorter
///     sequence are disregarded.
///
/// ## Examples:
///
/// Distance functions are provided as both standalone functions and methods via
/// the [`NucleotidesDistance`] trait. Using standalone functions:
/// ```
/// # use zoe::distance::dna::jukes_cantor_69;
/// let seq1: &[u8] = b"ACGT";
/// let seq2: &[u8] = b"ACGA";
///
/// let jc69_distance = jukes_cantor_69(seq1, seq2);
/// ```
///
/// Equivalently, using methods:
/// ```
/// # use zoe::{prelude::Nucleotides, distance::dna::NucleotidesDistance};
/// let seq1: Nucleotides = b"ACGT".into();
/// let seq2: Nucleotides = b"ACGA".into();
///
/// let jc69_distance = seq1.distance_jc69(&seq2);
/// ```
///
/// [`NucleotidesDistance`]: aa::AminoAcidsDistance
pub mod dna;

/// General distance functions for byte strings.
mod general;

pub use general::*;
