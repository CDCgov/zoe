use super::{
    AminoAcids, Nucleotides, NucleotidesMutable, NucleotidesView, NucleotidesViewMut, getter_traits::NucleotidesReadable,
};
use crate::{data::mappings::StdGeneticCode, private::Sealed};

/// Provides methods for translating nucleotides into amino acids.
pub trait Translate: NucleotidesReadable + Sealed {
    /// Translate the DNA sequence to [`AminoAcids`].
    #[inline]
    #[must_use]
    fn translate(&self) -> AminoAcids {
        AminoAcids(translate_sequence(self.nucleotide_bytes()))
    }

    /// Translates the DNA sequence to [`AminoAcids`] until the first stop codon
    /// in the current reading frame is found.
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

    /// Slides base by base over the sequence until the next translated `aa` is
    /// found and returns that index (or `None` otherwise).
    #[inline]
    #[must_use]
    fn find_next_aa(&self, aa: u8) -> Option<usize> {
        let needle = aa.to_ascii_uppercase();
        self.to_overlapping_aa_iter().position(|aa| aa == needle)
    }

    /// Slides codon by codon over the sequence (reading frame starting from 0)
    /// until the next translated `aa` is found and returns that index (or
    /// `None` otherwise).
    #[inline]
    #[must_use]
    fn find_next_aa_in_frame(&self, aa: u8) -> Option<usize> {
        let needle = aa.to_ascii_uppercase();
        self.to_aa_iter().position(|aa| aa == needle).map(|x| x * 3)
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

    /// Creates an iterator that translates amino acids base by base over all
    /// reading frames.
    #[inline]
    #[must_use]
    fn to_overlapping_aa_iter(&self) -> OverlappingCodonsIter {
        OverlappingCodonsIter::new(self.nucleotide_bytes())
    }
}

impl<T: NucleotidesReadable + Sealed> Translate for T {}

/// Iterator for translating [`Nucleotides`] into [`AminoAcids`] codon by codon
/// with a reading frame starting from 0.
///
/// This is created by [`to_aa_iter`] and [`to_aa_iter_with`]. To translate all
/// codons regardless of the reading frame, see [`OverlappingCodonsIter`].
///
/// [`Nucleotides`]: crate::prelude::Nucleotides
/// [`to_aa_iter`]: Translate::to_aa_iter
/// [`to_aa_iter_with`]: Translate::to_aa_iter_with
pub struct TranslatedNucleotidesIter<'a> {
    codons:                 std::slice::Iter<'a, [u8; 3]>,
    has_remainder:          bool,
    partial_codon_encoding: u8,
}

impl<'a> TranslatedNucleotidesIter<'a> {
    /// Create a new [`TranslatedNucleotidesIter`] using `~` to encode partial
    /// codons.
    #[inline]
    #[must_use]
    fn new(seq: &'a [u8]) -> Self {
        Self::new_with(seq, b'~')
    }

    /// Create a new [`TranslatedNucleotidesIter`], using
    /// `partial_codon_encoding` to encode partial codons.
    #[inline]
    #[must_use]
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

/// Iterator for translating [`Nucleotides`] into [`AminoAcids`] base by base
/// over all reading frames.
///
/// This is created by [`to_overlapping_aa_iter`]. To translate non-overlapping
/// codons in the current reading frame, see [`TranslatedNucleotidesIter`].
///
/// [`Nucleotides`]: crate::prelude::Nucleotides
/// [`to_overlapping_aa_iter`]: Translate::to_overlapping_aa_iter
pub struct OverlappingCodonsIter<'a> {
    codons: std::slice::Windows<'a, u8>,
}

impl<'a> OverlappingCodonsIter<'a> {
    /// Create a new [`OverlappingCodonsIter`] from the provided sequence.
    #[inline]
    #[must_use]
    fn new(seq: &'a [u8]) -> Self {
        Self { codons: seq.windows(3) }
    }
}

impl Iterator for OverlappingCodonsIter<'_> {
    type Item = u8;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.codons.next().map(StdGeneticCode::translate_codon)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.codons.size_hint()
    }
}

impl ExactSizeIterator for OverlappingCodonsIter<'_> {}
impl std::iter::FusedIterator for OverlappingCodonsIter<'_> {}

/// A codon is considered *partial* if it has fewer than 3 IUPAC bases
/// (non-gap).
#[inline]
#[must_use]
fn is_partial_codon(codon: [u8; 3]) -> bool {
    let count: u8 = codon.iter().map(|b| u8::from(*b == b'-' || *b == b'~' || *b == b'.')).sum();
    count == 1 || count == 2
}

/// Translates a byte slice into an amino acid byte vector.
#[inline]
#[must_use]
pub fn translate_sequence(s: &[u8]) -> Vec<u8> {
    // This was shown to be just as efficient as a manual implementation
    TranslatedNucleotidesIter::new(s).collect::<Vec<_>>()
}

/// Provides functionality to retrieve codons from a nucleotides sequences.
pub trait GetCodons: NucleotidesReadable + Sealed {
    /// Gets the bases grouped into codons as a slice of arrays, starting with
    /// the first base. Any trailing bases are included in the second tuple
    /// entry. To use a different reading frame, consider creating a view of the
    /// sequence with [`slice`].
    ///
    /// [`slice`]: crate::data::view_traits::Slice
    #[inline]
    #[must_use]
    fn as_codons(&self) -> (&[[u8; 3]], &[u8]) {
        self.nucleotide_bytes().as_chunks::<3>()
    }

    /// Gets the nth codon in the sequence with the reading frame starting with
    /// [`Self`]. To obtain a mutable reference, see [`nth_codon_mut`]. To avoid
    /// panicking when the index is out of bounds, see [`get_nth_codon`].
    ///
    /// ## Panics
    ///
    /// The function panics if the codon is out of bounds. Specifically, `3 *
    /// index + 3` must be at most `self.len()`.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// The indexing is *0-based*.
    ///
    /// </div>
    ///
    /// [`nth_codon_mut`]: Nucleotides::nth_codon_mut
    /// [`get_nth_codon`]: Nucleotides::get_nth_codon
    #[inline]
    #[must_use]
    fn nth_codon(&self, index: usize) -> [u8; 3] {
        self.nucleotide_bytes()
            .get(3 * index..3 * index + 3)
            .unwrap()
            .try_into()
            .unwrap()
    }

    /// Gets the nth non-overlapping codon in the sequence, returning `None` if
    /// it is out of bounds. To obtain a mutable reference, see
    /// [`get_nth_codon_mut`]. If it is known that the index is in bounds, see
    /// [`nth_codon`].
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// The indexing is *0-based*.
    ///
    /// </div>
    ///
    /// [`get_nth_codon_mut`]: Nucleotides::get_nth_codon_mut
    /// [`nth_codon`]: Nucleotides::nth_codon
    #[inline]
    #[must_use]
    fn get_nth_codon(&self, index: usize) -> Option<[u8; 3]> {
        // try_into will never fail
        self.nucleotide_bytes().get(3 * index..3 * index + 3)?.try_into().ok()
    }
}

impl GetCodons for Nucleotides {}
impl GetCodons for NucleotidesView<'_> {}
impl GetCodons for NucleotidesViewMut<'_> {}

/// Provides functionality to obtain mutable references to codons within
/// nucleotides sequences.
pub trait GetCodonsMut: NucleotidesMutable + Sealed {
    /// Gets the bases grouped into codons as a mutable slice of arrays,
    /// starting with the first base. Any trailing bases are included in the
    /// second tuple entry. To use a different reading frame, consider creating
    /// a view of the sequence with [`slice_mut`].
    ///
    /// [`slice_mut`]: crate::data::view_traits::SliceMut
    #[inline]
    #[must_use]
    fn as_codons_mut(&mut self) -> (&mut [[u8; 3]], &mut [u8]) {
        self.nucleotide_mut_bytes().as_chunks_mut::<3>()
    }

    /// Gets a mutable reference to the nth codon in the sequence with the
    /// reading frame starting with [`Self`]. When mutability is not needed, use
    /// [`nth_codon`]. To avoid panicking when the index is out of bounds, see
    /// [`get_nth_codon_mut`].
    ///
    /// ## Panics
    ///
    /// The function panics if the codon is out of bounds. Specifically, `3 *
    /// index + 3` must be at most `self.len()`.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// The indexing is *0-based*.
    ///
    /// </div>
    ///
    /// [`nth_codon`]: Nucleotides::nth_codon
    /// [`get_nth_codon_mut`]: Nucleotides::get_nth_codon_mut
    #[inline]
    #[must_use]
    fn nth_codon_mut(&mut self, index: usize) -> &mut [u8; 3] {
        self.nucleotide_mut_bytes()
            .get_mut(3 * index..3 * index + 3)
            .unwrap()
            .try_into()
            .unwrap()
    }

    /// Gets a mutable reference to the nth codon in the sequence, returning
    /// `None` if it is out of bounds. When mutability is not needed, use
    /// [`get_nth_codon`]. If it is known that the index is in bounds, see
    /// [`nth_codon_mut`].
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// The indexing is *0-based*.
    ///
    /// </div>
    ///
    /// [`get_nth_codon`]: Nucleotides::get_nth_codon
    /// [`nth_codon_mut`]: Nucleotides::nth_codon_mut
    #[inline]
    #[must_use]
    fn get_nth_codon_mut(&mut self, index: usize) -> Option<&mut [u8; 3]> {
        self.nucleotide_mut_bytes().get_mut(3 * index..3 * index + 3)?.try_into().ok()
    }
}

impl GetCodonsMut for Nucleotides {}
impl GetCodonsMut for NucleotidesViewMut<'_> {}
