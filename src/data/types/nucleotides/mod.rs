use crate::{
    alignment::{LocalProfiles, SharedProfiles},
    data::{err::QueryProfileError, matrices::WeightMatrix},
    prelude::*,
};
use std::simd::{LaneCount, SupportedLaneCount};

/// Getters to help with trait implementations.
mod getter_traits;
/// Reverse complement implementations.
mod rev_comp;
/// Traits and methods for sanitizing and validating nucleotide data.
mod sanitize;
/// Standard library implementations.
mod std_traits;
/// Amino acid translation implementations.
mod translation;
/// Nucleotide view types.
mod view_traits;

pub use getter_traits::*;
pub use rev_comp::*;
pub use sanitize::*;
pub use translation::*;

#[cfg(test)]
mod bench;

#[cfg(test)]
mod test;

/// [`Nucleotides`] is a transparent, new-type wrapper around [`Vec<u8>`] that
/// provides DNA-specific functionality and semantics. It may contain either
/// aligned or unaligned valid IUPAC letters.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[repr(transparent)]
pub struct Nucleotides(pub(crate) Vec<u8>);

/// The corresponding immutable view type for [`Nucleotides`]. See
/// [Views](crate::data#views) for more details.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[repr(transparent)]
pub struct NucleotidesView<'a>(pub(crate) &'a [u8]);

/// The corresponding mutable view type for [`Nucleotides`]. See
/// [Views](crate::data#views) for more details.
#[derive(Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[repr(transparent)]
pub struct NucleotidesViewMut<'a>(pub(crate) &'a mut [u8]);

impl Nucleotides {
    // Conversions and indexing

    /// Creates a new [`Nucleotides`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Nucleotides(Vec::new())
    }

    /// Consumes a [`Vec<u8>`] and return [`Nucleotides`] without checking
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_vec_unchecked(v: Vec<u8>) -> Self {
        Nucleotides(v)
    }

    /// Consumes [`Nucleotides`] and returns a [`Vec<u8>`].
    #[inline]
    #[must_use]
    pub fn into_vec(self) -> Vec<u8> {
        self.0
    }

    /// Gets the nucleotides as a byte slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    /// Gets the nucleotides as a mutable byte slice.
    #[inline]
    #[must_use]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        &mut self.0
    }

    /// Gets the base or byte slice at the zero-based index, returning an
    /// [`Option`].
    #[inline]
    #[must_use]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// Creates an iterator over the nucleotides as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }

    /// Creates an iterator over the nucleotides as `&mut u8`.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, u8> {
        self.0.iter_mut()
    }

    // Manipulation

    /// Truncates the length of the sequence to the specified `new_length`. This
    /// is equivalent to 3' trimming up to and including the index.
    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
    }

    /// Cuts the 5' end of the [`Nucleotides`] just prior to the new starting
    /// index (0-based). Be aware that this clones the internal buffer!
    #[inline]
    pub fn cut_to_start(&mut self, new_start: usize) {
        *self = Nucleotides(self.0.drain(new_start..).collect());
    }

    // Nucleotides-specific methods

    /// Creates a [`LocalProfiles`] for alignment.
    ///
    /// # Errors
    ///
    /// Returns an [`QueryProfileError`] if the profile creation fails due to
    /// invalid sequence data or unsupported parameters.
    ///
    /// # Example
    ///
    /// ```
    /// # use zoe::{data::WeightMatrix, prelude::Nucleotides};
    /// let reference: &[u8] = b"GGCCACAGGATTGAG";
    /// let query: Nucleotides = b"CTCAGATTG".into();
    ///
    /// const WEIGHTS: WeightMatrix<u8, 5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = query.into_local_profile::<32, 5>(&WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_u8(reference).unwrap();
    /// assert_eq!(score, 27);
    /// ```
    #[inline]
    pub fn into_local_profile<'a, 'b, const N: usize, const S: usize>(
        &'b self, matrix: &'a WeightMatrix<u8, S>, gap_open: u8, gap_extend: u8,
    ) -> Result<LocalProfiles<'a, N, S>, QueryProfileError>
    where
        LaneCount<N>: SupportedLaneCount,
        'b: 'a, {
        LocalProfiles::new(&self.0, matrix, gap_open, gap_extend)
    }

    /// Creates a [`SharedProfiles`] for alignment.
    ///
    /// # Errors
    ///
    /// Returns an [`QueryProfileError`] if the profile creation fails due to invalid
    /// sequence data or unsupported parameters.
    ///
    /// # Example
    ///
    /// ```
    /// # use zoe::{data::WeightMatrix, prelude::Nucleotides};
    /// let reference: &[u8] = b"GGCCACAGGATTGAG";
    /// let query: Nucleotides = b"CTCAGATTG".into();
    ///
    /// const WEIGHTS: WeightMatrix<u8,5> = WeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = query.into_shared_profile::<32, 5>(&WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_u8(reference).unwrap();
    /// assert_eq!(score, 27);
    /// ```
    #[inline]
    pub fn into_shared_profile<'a, const N: usize, const S: usize>(
        &self, matrix: &'a WeightMatrix<u8, S>, gap_open: u8, gap_extend: u8,
    ) -> Result<SharedProfiles<'a, N, S>, QueryProfileError>
    where
        LaneCount<N>: SupportedLaneCount, {
        SharedProfiles::new(self.as_bytes().into(), matrix, gap_open, gap_extend)
    }

    /// Returns the reverse complement of the sequence as a new record.
    #[inline]
    #[must_use]
    pub fn to_reverse_complement(&self) -> Nucleotides {
        Nucleotides(reverse_complement(&self.0))
    }

    /// Computes the reverse complement of the sequence in-place.
    #[inline]
    pub fn make_reverse_complement(&mut self) {
        make_reverse_complement(&mut self.0);
    }

    // Associated functions

    /// Generates a random DNA sequence of given `length` and using a random
    /// `seed`.  Canonical DNA only contains A, C, G, or T. Requires `rand`
    /// feature to be enabled.
    #[cfg(feature = "rand")]
    #[must_use]
    pub fn generate_random_dna(length: usize, seed: u64) -> Self {
        Nucleotides(crate::generate::rand_sequence(b"AGCT", length, seed))
    }
}

impl<'a> NucleotidesView<'a> {
    // Conversions and indexing

    /// Creates a new [`NucleotidesView`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        NucleotidesView(&[])
    }

    /// Creates a [`NucleotidesView`] from a byte slice without checking
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_bytes_unchecked(v: &'a [u8]) -> Self {
        NucleotidesView(v)
    }

    /// Gets the nucleotides as a byte slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0
    }

    /// Gets the base or byte slice at the zero-based index, returning an
    /// [`Option`].
    #[inline]
    #[must_use]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// Creates an iterator over the nucleotides as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }

    // Nucleotides-specific methods

    /// Returns the reverse complement of the sequence as a new record.
    #[inline]
    #[must_use]
    pub fn to_reverse_complement(&self) -> Nucleotides {
        Nucleotides(reverse_complement(self.0))
    }
}

impl<'a> NucleotidesViewMut<'a> {
    // Conversions and indexing

    /// Creates a new [`NucleotidesViewMut`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        NucleotidesViewMut(&mut [])
    }

    /// Creates a [`NucleotidesViewMut`] from a byte slice without checking
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_bytes_unchecked(v: &'a mut [u8]) -> Self {
        NucleotidesViewMut(v)
    }

    /// Gets the nucleotides as a byte slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0
    }

    /// Gets the nucleotides as a mutable byte slice.
    #[inline]
    #[must_use]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        self.0
    }

    /// Gets the base or byte slice at the zero-based index, returning an
    /// [`Option`].
    #[inline]
    #[must_use]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.as_ref().get(index)
    }

    /// Create an iterator over the nucleotides as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }

    /// Create an iterator over the nucleotides as `&mut u8`.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, u8> {
        self.0.iter_mut()
    }

    // Nucleotides-specific methods

    /// Returns the reverse complement of the sequence as a new record.
    #[inline]
    #[must_use]
    pub fn to_reverse_complement(&self) -> Nucleotides {
        Nucleotides(reverse_complement(self.0))
    }

    /// Computes the reverse complement of the sequence in-place.
    #[inline]
    pub fn make_reverse_complement(&mut self) {
        make_reverse_complement(self.0);
    }
}
