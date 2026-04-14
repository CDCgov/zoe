use std::ops::Range;

use super::{
    AminoAcids, Nucleotides, NucleotidesMutable, NucleotidesView, NucleotidesViewMut, getter_traits::NucleotidesReadable,
};
use crate::{
    DEFAULT_SIMD_LANES, data::mappings::StdGeneticCode, private::Sealed, search::find_mapped_match_simd,
    simd::SimdByteFunctions,
};
use std::simd::prelude::*;

/// Provides methods for translating nucleotides into amino acids.
pub trait Translate: NucleotidesReadable + Sealed {
    /// Translates the DNA sequence to [`AminoAcids`].
    ///
    /// This uses a [`TranslatedNucleotidesIter`] for translation. Partial
    /// codons are translated to `~`. See [`TranslatedNucleotidesIter`] and
    /// [`StdGeneticCode`] for more details.
    #[inline]
    #[must_use]
    fn translate(&self) -> AminoAcids {
        AminoAcids(translate_sequence(self.nucleotide_bytes()))
    }

    /// Translates the DNA sequence to [`AminoAcids`] until the first stop codon
    /// in the current reading frame is found.
    ///
    /// The final stop codon (`*`) is included if present. To exclude this, call
    /// [`chop_stop`] on the output.
    ///
    /// Partial codons are translated to `~`. See [`TranslatedNucleotidesIter`]
    /// and [`StdGeneticCode`] for more details.
    ///
    /// [`chop_stop`]: AminoAcids::chop_stop
    #[inline]
    #[must_use]
    fn translate_to_stop(&self) -> AminoAcids {
        let aa_iter = self.to_aa_iter();
        let mut out = Vec::with_capacity(aa_iter.len());
        for aa in aa_iter {
            out.push(aa);
            if aa == b'*' {
                break;
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

    /// Returns the starting index of the first case-insensitive start codon
    /// (`ATG` or `AUG`) or `None` otherwise.
    ///
    /// If the sequence is known to not contain `U` and be of a particular case,
    /// then [`find_substring`] may be faster. This function uses
    /// [`find_mapped_match_simd`] in its implementation.
    ///
    /// [`find_substring`]: crate::search::ByteSubstring::find_substring
    #[inline]
    #[must_use]
    fn find_start_codon(&self) -> Option<Range<usize>> {
        const N: usize = DEFAULT_SIMD_LANES;

        fn simd_map(s: Simd<u8, N>) -> Simd<u8, N> {
            let mut xformed = s.to_ascii_uppercase();
            xformed.if_value_then_replace(b'U', b'T');
            xformed
        }

        fn byte_map(b: u8) -> u8 {
            let mut xformed = b.to_ascii_uppercase();
            if xformed == b'U' {
                xformed = b'T';
            }
            xformed
        }

        find_mapped_match_simd::<N, _, _>(self.nucleotide_bytes(), b"ATG", simd_map, byte_map).map(|s| s..s + 3)
    }

    /// Creates an iterator for [`AminoAcids`] translation.
    ///
    /// Partial codons are translated to `~`. See [`TranslatedNucleotidesIter`]
    /// and [`StdGeneticCode`] for more details.
    #[inline]
    #[must_use]
    fn to_aa_iter(&self) -> TranslatedNucleotidesIter<'_> {
        TranslatedNucleotidesIter::new(self.nucleotide_bytes())
    }

    /// Creates an iterator for [`AminoAcids`] translation, ending iteration
    /// when fewer than three bases remain.
    ///
    /// Partial codons are translated to `~`, except for a partial codon arising
    /// at the end of the sequence due to having fewer than three bases
    /// remaining, which is excluded. See [`TranslatedNucleotidesIter`] and
    /// [`StdGeneticCode`] for more details.
    #[inline]
    #[must_use]
    fn to_aa_iter_exact(&self) -> TranslatedNucleotidesIter<'_> {
        TranslatedNucleotidesIter::new_exact(self.nucleotide_bytes())
    }

    /// Creates an iterator for [`AminoAcids`] translation, using
    /// `partial_codon_encoding` to encode partial codons.
    ///
    /// See [`TranslatedNucleotidesIter`] and [`StdGeneticCode`] for more
    /// details.
    #[inline]
    #[must_use]
    fn to_aa_iter_with(&self, partial_codon_encoding: u8) -> TranslatedNucleotidesIter<'_> {
        TranslatedNucleotidesIter::new_with(self.nucleotide_bytes(), partial_codon_encoding)
    }

    /// Creates an iterator for [`AminoAcids`] translation, ending iteration
    /// when fewer than three bases remain, and using `partial_codon_encoding`
    /// to encode other partial codons.
    ///
    /// See [`TranslatedNucleotidesIter`] and [`StdGeneticCode`] for more
    /// details.
    #[inline]
    #[must_use]
    fn to_aa_iter_exact_with(&self, partial_codon_encoding: u8) -> TranslatedNucleotidesIter<'_> {
        TranslatedNucleotidesIter::new_exact_with(self.nucleotide_bytes(), partial_codon_encoding)
    }

    /// Creates an iterator that translates amino acids base by base over all
    /// reading frames.
    ///
    /// See [`OverlappingCodonsIter`] and [`StdGeneticCode`] for more details.
    #[inline]
    #[must_use]
    fn to_overlapping_aa_iter(&self) -> OverlappingCodonsIter<'_> {
        OverlappingCodonsIter::new(self.nucleotide_bytes())
    }
}

impl<T: NucleotidesReadable + Sealed> Translate for T {}

/// Iterator for translating [`Nucleotides`] into [`AminoAcids`] codon by codon
/// with a reading frame starting from 0.
///
/// This is created by several methods of [`Translate`], such as [`to_aa_iter`].
/// To translate all codons regardless of the reading frame, see
/// [`OverlappingCodonsIter`].
///
/// For details about the genetic code used, see [`StdGeneticCode`]. This
/// iterator will also translate *partial* codons to `~` (or a custom provided
/// character, if [`to_aa_iter_with`] or [`to_aa_iter_exact_with`] is used). A
/// codon is considered partial if at least one of the following is true:
///
/// - The codon does not contain three bases (which will occur if the input
///   sequence contains a number of bases not divisible by three).
/// - The codon contains one or two of the characters `~`, `-`, or `.` among its
///   three bases.
///
/// To avoid the first case of partial codons (and instead end iteration), use
/// [`to_aa_iter_exact`] or [`to_aa_iter_exact_with`].
///
/// [`Nucleotides`]: crate::prelude::Nucleotides
/// [`to_aa_iter`]: Translate::to_aa_iter
/// [`to_aa_iter_exact`]: Translate::to_aa_iter_exact
/// [`to_aa_iter_with`]: Translate::to_aa_iter_with
/// [`to_aa_iter_exact_with`]: Translate::to_aa_iter_exact_with
pub struct TranslatedNucleotidesIter<'a> {
    codons:                 std::slice::Iter<'a, [u8; 3]>,
    has_remainder:          bool,
    partial_codon_encoding: u8,
}

impl<'a> TranslatedNucleotidesIter<'a> {
    /// Creates a new [`TranslatedNucleotidesIter`] using `~` to encode partial
    /// codons.
    ///
    /// See the documentation for [`TranslatedNucleotidesIter`] for more details
    /// on partial codons.
    #[inline]
    #[must_use]
    fn new(seq: &'a [u8]) -> Self {
        Self::new_with(seq, b'~')
    }

    /// Creates a new [`TranslatedNucleotidesIter`] using `~` to encode partial
    /// codons, as well as aborting iteration when fewer than three bases remain
    /// (rather than issuing a partial codon at the end).
    ///
    /// See the documentation for [`TranslatedNucleotidesIter`] for more details
    /// on partial codons.
    #[inline]
    #[must_use]
    fn new_exact(seq: &'a [u8]) -> Self {
        Self::new_exact_with(seq, b'~')
    }

    /// Creates a new [`TranslatedNucleotidesIter`], using
    /// `partial_codon_encoding` to encode partial codons.
    ///
    /// See the documentation for [`TranslatedNucleotidesIter`] for more details
    /// on partial codons.
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

    /// Creates a new [`TranslatedNucleotidesIter`] using
    /// `partial_codon_encoding` to encode partial codons, as well as aborting
    /// iteration when fewer than three bases remain (rather than issuing a
    /// partial codon at the end).
    ///
    /// See the documentation for [`TranslatedNucleotidesIter`] for more details
    /// on partial codons.
    #[inline]
    #[must_use]
    fn new_exact_with(seq: &'a [u8], partial_codon_encoding: u8) -> Self {
        let codons = seq.as_chunks::<3>().0.iter();
        Self {
            codons,
            has_remainder: false,
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
/// For details about the genetic code used, see [`StdGeneticCode`].
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
    pub(crate) fn new(seq: &'a [u8]) -> Self {
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

/// Checks whether a codon is *partial* (if it has one or two bases that are
/// `-`, `~`, or `.`).
#[inline]
#[must_use]
fn is_partial_codon(codon: [u8; 3]) -> bool {
    let count: u8 = codon.iter().map(|b| u8::from(*b == b'-' || *b == b'~' || *b == b'.')).sum();
    count == 1 || count == 2
}

/// Translates a byte slice into an amino acid byte vector.
///
/// This uses a [`TranslatedNucleotidesIter`] for translation. Partial codons
/// are translated to `~`. See [`TranslatedNucleotidesIter`] and
/// [`StdGeneticCode`] for more details.
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
    /// [`slice`]: crate::data::views::Slice
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

    /// Gets the first 3 bases as a codon, returning `None` if no such codon
    /// exists.
    #[inline]
    #[must_use]
    fn get_first_codon(&self) -> Option<[u8; 3]> {
        self.nucleotide_bytes().first_chunk().copied()
    }

    /// Gets the last in-frame codon (from start) for the sequence, returning
    /// `None` if no such codon exists.
    #[inline]
    #[must_use]
    fn get_last_codon(&self) -> Option<[u8; 3]> {
        self.as_codons().0.last().copied()
    }

    /// Gets the last 3 bases as a codon regardless of frame, returning `None`
    /// if no such codon exists.
    #[inline]
    #[must_use]
    fn get_tail_codon(&self) -> Option<[u8; 3]> {
        self.nucleotide_bytes().last_chunk().copied()
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
    /// [`slice_mut`]: crate::data::views::SliceMut
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

    /// Gets a mutable reference to the first 3 bases as a codon, returning
    /// `None` if no such codon exists.
    #[inline]
    #[must_use]
    fn get_first_codon_mut(&mut self) -> Option<&mut [u8; 3]> {
        self.nucleotide_mut_bytes().first_chunk_mut()
    }

    /// Gets a mutable reference to the last in-frame (from start) codon for the
    /// sequence relative to, returning `None` if no such codon exists.
    #[inline]
    #[must_use]
    fn get_last_codon_mut(&mut self) -> Option<&mut [u8; 3]> {
        self.as_codons_mut().0.last_mut()
    }

    /// Gets a mutable reference to the last 3 bases as a codon regardless of
    /// frame, returning `None` if no such codon exists.
    #[inline]
    #[must_use]
    fn get_tail_codon_mut(&mut self) -> Option<&mut [u8; 3]> {
        self.nucleotide_mut_bytes().last_chunk_mut()
    }
}

impl GetCodonsMut for Nucleotides {}
impl GetCodonsMut for NucleotidesViewMut<'_> {}

/// Extension trait for codon-like objects (`[u8; 3]` in *Zoe*) providing common
/// classification methods.
pub trait CodonExtension {
    /// Returns `true` if the codon is a standard stop codon (`TAA`, `TAG`,
    /// `TGA`, or RNA equivalents including some ambiguous forms).
    fn is_std_stop_codon(&self) -> bool;

    /// Returns `true` if the codon translates to an amino acid.
    ///
    /// Ambiguous codons are supported, but they must unambiguously translate to
    /// an amino acid in order to return `true`. Stop codons, partial codons,
    /// deletions (e.g., `---`), and the ambiguous codon `NNN` all return
    /// `false`.
    fn is_amino_acid(&self) -> bool;

    /// Returns `true` if the codon can be resolved to a translation (amino
    /// acid, deletion, missing codon, or stop codon) under the standard genetic
    /// code.
    fn is_resolvable_codon(&self) -> bool;

    /// Returns `true` if the codon is *partial* (contains one or two gap-like
    /// characters: `-`, `~`, or `.`).
    fn is_partial_codon(&self) -> bool;
}

impl CodonExtension for [u8; 3] {
    #[inline]
    fn is_std_stop_codon(&self) -> bool {
        StdGeneticCode::is_stop_codon(self)
    }

    #[inline]
    fn is_amino_acid(&self) -> bool {
        matches!(StdGeneticCode::get(self), Some(aa) if aa.is_ascii_alphabetic() && aa != b'X')
    }

    #[inline]
    fn is_resolvable_codon(&self) -> bool {
        StdGeneticCode::get(self).is_some()
    }

    #[inline]
    fn is_partial_codon(&self) -> bool {
        is_partial_codon(*self)
    }
}
