use std::slice::ChunksExact;

use crate::data::mappings::StdGeneticCode;

use super::{AminoAcids, getter_traits::NucleotidesReadable};

/// Provides methods for translating nucleotides into amino acids.
pub trait Translate: NucleotidesReadable {
    /// Translate the DNA sequence to [`AminoAcids`].
    #[inline]
    #[must_use]
    fn translate(&self) -> AminoAcids {
        AminoAcids(translate_sequence(self.nucleotide_bytes()))
    }

    /// Creates an iterator for [`AminoAcids`] translation.
    #[inline]
    #[must_use]
    fn to_aa_iter(&self) -> TranslatedNucleotidesIter {
        TranslatedNucleotidesIter {
            codons:        self.nucleotide_bytes().chunks_exact(3),
            has_remainder: true,
        }
    }
}

impl<T: NucleotidesReadable> Translate for T {}

/// Iterator for translating [`Nucleotides`] into [`AminoAcids`].
///
/// [`Nucleotides`]: crate::prelude::Nucleotides
pub struct TranslatedNucleotidesIter<'a> {
    codons:        ChunksExact<'a, u8>,
    has_remainder: bool,
}

impl Iterator for TranslatedNucleotidesIter<'_> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(codon) = self.codons.next() {
            if is_partial_codon(codon) {
                Some(b'~')
            } else {
                Some(StdGeneticCode::get(codon).unwrap_or(b'X'))
            }
        } else if self.has_remainder && is_partial_codon(self.codons.remainder()) {
            self.has_remainder = false;
            Some(b'~')
        } else {
            None
        }
    }
}

/// A codon is considered *partial* if it has fewer than 3 IUPAC bases.
#[inline]
#[must_use]
fn is_partial_codon(codon: &[u8]) -> bool {
    if codon.is_empty() {
        false
    } else if codon.len() < 3 {
        true
    } else {
        let count: u8 = codon.iter().map(|b| u8::from(*b == b'-' || *b == b'~' || *b == b'.')).sum();
        count == 1 || count == 2
    }
}

/// Translates a byte slice into an amino acid byte vector.
#[inline]
#[must_use]
pub fn translate_sequence(s: &[u8]) -> Vec<u8> {
    // TODO: this and the translation iterator can be made into array chunks
    // for further performance, but best to wait until the feature is more
    // likely to be adopted and/or needed by other functions in Zoe.
    let mut codons = s.chunks_exact(3);
    let mut aa_sequence = Vec::with_capacity(s.len() / 3 + 1);

    for codon in codons.by_ref() {
        aa_sequence.push(if is_partial_codon(codon) {
            b'~'
        } else {
            StdGeneticCode::get(codon).unwrap_or(b'X')
        });
    }

    let tail = codons.remainder();
    if is_partial_codon(tail) {
        aa_sequence.push(b'~');
    }

    aa_sequence
}
