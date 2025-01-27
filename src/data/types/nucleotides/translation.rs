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

    #[inline]
    #[must_use]
    fn translate_to_stop(&self) -> AminoAcids {
        let aa_iter = self.to_aa_iter();
        let mut out = Vec::with_capacity(aa_iter.len());
        for aa in aa_iter {
            out.push(aa);
            if aa == b'*' {
                return AminoAcids(out);
            }
        }
        AminoAcids(out)
    }

    /// Creates an iterator for [`AminoAcids`] translation.
    #[inline]
    #[must_use]
    fn to_aa_iter(&self) -> TranslatedNucleotidesIter {
        TranslatedNucleotidesIter::new(self.nucleotide_bytes())
    }
}

impl<T: NucleotidesReadable> Translate for T {}

/// Iterator for translating [`Nucleotides`] into [`AminoAcids`].
///
/// [`Nucleotides`]: crate::prelude::Nucleotides
pub struct TranslatedNucleotidesIter<'a> {
    // TODO: this can be made into array chunks for further performance, but
    // best to wait until the feature is more likely to be adopted and/or needed
    // by other functions in Zoe.
    codons:        ChunksExact<'a, u8>,
    has_remainder: bool,
}

impl<'a> TranslatedNucleotidesIter<'a> {
    #[inline]
    fn new(seq: &'a [u8]) -> Self {
        let codons = seq.chunks_exact(3);
        let has_remainder = !codons.remainder().is_empty();
        Self { codons, has_remainder }
    }
}

impl Iterator for TranslatedNucleotidesIter<'_> {
    type Item = u8;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(codon) = self.codons.next() {
            if is_partial_codon(codon) {
                Some(b'~')
            } else {
                Some(StdGeneticCode::translate_codon(codon))
            }
        } else if self.has_remainder {
            self.has_remainder = false;
            Some(b'~')
        } else {
            None
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = usize::from(self.has_remainder) + self.codons.size_hint().0;
        (size, Some(size))
    }
}

impl ExactSizeIterator for TranslatedNucleotidesIter<'_> {}

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
    // This was shown to be just as efficient as a manual implementation
    TranslatedNucleotidesIter::new(s).collect::<Vec<_>>()
}
