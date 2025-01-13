/// Amino acid distance functions.
pub mod aa;
/// Various nucleotide substitution models for calculating evolutionary distances between two aligned DNA sequences.
///
/// ## Assumptions:
///
/// * __Alignment:__ both sequences must be aligned.
/// * __Pairwise deletion:__ If a gap or ambiguous (non-ACGT) base is present in any position of the aligned sequences,
///     the bases in that position are excluded from the analysis for both sequences.
/// * __Unequal sequence lengths:__ If sequences are of different lengths, bases in the longer sequence that extend beyond
///     the length of the shorter sequence are disregarded.
///
/// ## Example:
/// ```
/// # use zoe::distance::dna::jukes_cantor_69;
/// let seq1: &[u8] = b"ACGT";
/// let seq2: &[u8] = b"ACGA";
///
/// let jc69_distance = jukes_cantor_69(seq1, seq2);
/// ```
pub mod dna;

/// General string-based distance functions.
mod general;

pub use general::*;
