use super::{AminoAcids, getter_traits::NucleotidesReadable};
use crate::data::mappings::StdGeneticCode;

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

    /// Creates an iterator for [`AminoAcids`] translation. Takes an amino acid
    /// encoding for partial codons.
    #[inline]
    #[must_use]
    fn to_aa_iter_with(&self, partial_codon_encoding: u8) -> TranslatedNucleotidesIter {
        TranslatedNucleotidesIter::new_with(self.nucleotide_bytes(), partial_codon_encoding)
    }
}

impl<T: NucleotidesReadable> Translate for T {}

/// Iterator for translating [`Nucleotides`] into [`AminoAcids`].
///
/// [`Nucleotides`]: crate::prelude::Nucleotides
pub struct TranslatedNucleotidesIter<'a> {
    codons:                 std::slice::Iter<'a, [u8; 3]>,
    has_remainder:          bool,
    partial_codon_encoding: u8,
}

impl<'a> TranslatedNucleotidesIter<'a> {
    #[inline]
    fn new(seq: &'a [u8]) -> Self {
        Self::new_with(seq, b'~')
    }

    #[inline]
    fn new_with(seq: &'a [u8], partial_codon_encoding: u8) -> Self {
        let (codons, has_remainder) = {
            let (chunks, remainder) = seq.as_chunks::<3>();
            (chunks.iter(), !remainder.is_empty())
        };

        Self {
            codons,
            has_remainder,
            partial_codon_encoding,
        }
    }
}

impl Iterator for TranslatedNucleotidesIter<'_> {
    type Item = u8;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(codon) = self.codons.next() {
            if is_partial_codon(*codon) {
                Some(self.partial_codon_encoding)
            } else {
                Some(StdGeneticCode::translate_codon(codon))
            }
        } else if self.has_remainder {
            self.has_remainder = false;
            Some(self.partial_codon_encoding)
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
impl std::iter::FusedIterator for TranslatedNucleotidesIter<'_> {}

/// A codon is considered *partial* if it has fewer than 3 IUPAC bases.
#[inline]
#[must_use]
fn is_partial_codon(codon: [u8; 3]) -> bool {
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
