use crate::data::types::impl_std_traits_for_sequence;

mod getter_traits;
mod view_traits;

pub use getter_traits::*;

/// [`AminoAcids`] is a transparent, new-type wrapper around [`Vec<u8>`] that
/// provides protein-specific functionality and semantics.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[repr(transparent)]
pub struct AminoAcids(pub(crate) Vec<u8>);

/// The corresponding immutable view type for [`AminoAcids`]. See
/// [Views](crate::data#views) for more details.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[repr(transparent)]
pub struct AminoAcidsView<'a>(pub(crate) &'a [u8]);

/// The corresponding mutable view type for [`AminoAcids`]. See
/// [Views](crate::data#views) for more details.
#[derive(Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[repr(transparent)]
pub struct AminoAcidsViewMut<'a>(pub(crate) &'a mut [u8]);

impl AminoAcids {
    // Conversions and indexing

    /// Creates a new [`AminoAcids`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        AminoAcids(Vec::new())
    }

    /// Consumes a [`Vec<u8>`] and return [`AminoAcids`] without checking for
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_vec_unchecked(v: Vec<u8>) -> Self {
        AminoAcids(v)
    }

    /// Consumes [`AminoAcids`] and returns a [`Vec<u8>`].
    #[inline]
    #[must_use]
    pub fn into_vec(self) -> Vec<u8> {
        self.0
    }

    /// Gets the amino acids as a byte slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    /// Gets the amino acids as a mutable byte slice.
    #[inline]
    #[must_use]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        &mut self.0
    }

    /// Gets the amino acids or byte slice at the zero-based index, returning an
    /// [`Option`].
    #[inline]
    #[must_use]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// Creates an iterator over the amino acids as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }

    /// Creates an iterator over the amino acids as `&mut u8`.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, u8> {
        self.0.iter_mut()
    }

    // Manipulation

    /// Truncates the length of the sequence to the specified `new_length`. This
    /// is equivalent to C-terminus trimming up to and including the index.
    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
    }

    /// Cuts the N-terminus of the [`AminoAcids`] just prior to the new starting
    /// index (0-based).
    ///
    /// This method clones the internal buffer.
    #[inline]
    pub fn cut_to_start(&mut self, new_start: usize) {
        *self = Self(self.0.drain(new_start..).collect());
    }

    /// Extends amino acids sequence by the given slice.
    #[inline]
    pub fn extend_from_slice(&mut self, slice: impl AsRef<[u8]>) {
        self.0.extend_from_slice(slice.as_ref());
    }

    /// Prepends the amino acid sequence with the given slice.
    ///
    /// This method allocates.
    #[inline]
    pub fn prepend_from_slice(&mut self, slice: impl AsRef<[u8]>) {
        let slice = slice.as_ref();
        if slice.is_empty() {
            return;
        }
        let mut new = Vec::with_capacity(slice.len() + self.0.len());
        new.extend_from_slice(slice);
        new.extend_from_slice(&self.0);
        *self = AminoAcids(new);
    }

    /// Pads the C-terminus (end) of the [`AminoAcids`] by `aa` for `count`
    /// times.
    #[inline]
    pub fn pad_end(&mut self, aa: u8, count: usize) {
        self.0.extend(std::iter::repeat_n(aa, count));
    }

    /// Pads the N-terminus (start) of the [`AminoAcids`] by `aa` for `count`
    /// times.
    ///
    /// This method allocates.
    #[inline]
    pub fn pad_start(&mut self, aa: u8, count: usize) {
        if count == 0 {
            return;
        }
        let mut new = Vec::with_capacity(count + self.0.len());
        new.extend(std::iter::repeat_n(aa, count));
        new.extend_from_slice(&self.0);
        *self = AminoAcids(new);
    }

    /// Pads the [`AminoAcids`] on both flanking sides by `aa` for
    /// `amount_before` and `amount_after` respectively.
    ///
    /// This method allocates.
    #[inline]
    pub fn pad_both_sides(&mut self, aa: u8, amount_before: usize, amount_after: usize) {
        if amount_before == 0 && amount_after == 0 {
            return;
        }
        let mut new = Vec::with_capacity(amount_before + amount_after + self.0.len());
        new.extend(std::iter::repeat_n(aa, amount_before));
        new.extend_from_slice(&self.0);
        new.extend(std::iter::repeat_n(aa, amount_after));
        *self = AminoAcids(new);
    }

    /// Clears the sequence.
    #[inline]
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Removes the final stop codon if present.
    ///
    /// This method takes and gives ownership to support chaining.
    #[inline]
    #[must_use]
    pub fn chop_stop(mut self) -> Self {
        if let Some(&b'*') = self.0.last() {
            self.0.pop();
            self
        } else {
            self
        }
    }

    // Associated functions

    /// Generates a random amino acid sequence from the 20 uppercase canonical
    /// amino acids.
    ///
    /// The sequence is of the specified length, and a seed is passed for
    /// reproducibility.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *rand* feature in your `Cargo.toml` to use this
    /// function.
    ///
    /// </div>
    #[must_use]
    #[cfg(feature = "rand")]
    pub fn generate_random_aa(length: usize, seed: u64) -> Self {
        AminoAcids(crate::generate::rand_sequence(
            crate::data::constants::alphas::AA_IUPAC_NO_GAPS_UC,
            length,
            seed,
        ))
    }
}

impl<'a> AminoAcidsView<'a> {
    // Conversions and indexing

    /// Creates a new [`AminoAcidsView`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        AminoAcidsView(&[])
    }

    /// Creates a [`AminoAcidsView`] from a byte slice without checking
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_bytes_unchecked(v: &'a [u8]) -> Self {
        AminoAcidsView(v)
    }

    /// Obtains the bytes as a slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &'a [u8] {
        self.0
    }

    /// Gets the amino acids or byte slice at the zero-based index, returning an
    /// [`Option`].
    #[inline]
    pub fn get<I>(&self, index: I) -> Option<&'a I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// Creates an iterator over the amino acids as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'a, u8> {
        self.0.iter()
    }
}

impl<'a> AminoAcidsViewMut<'a> {
    // Conversions and indexing

    /// Creates a new [`AminoAcidsViewMut`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        AminoAcidsViewMut(&mut [])
    }

    /// Creates a [`AminoAcidsViewMut`] from a byte slice without checking
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_bytes_unchecked(v: &'a mut [u8]) -> Self {
        AminoAcidsViewMut(v)
    }

    /// Gets the amino acids as a byte slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0
    }

    /// Gets the amino acids as a mutable byte slice.
    #[inline]
    #[must_use]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        self.0
    }

    /// Gets the amino acids or byte slice at the zero-based index, returning an
    /// [`Option`].
    #[inline]
    #[must_use]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// Create an iterator over the amino acids as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }

    /// Create an iterator over the amino acids as `&mut u8`.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, u8> {
        self.0.iter_mut()
    }
}

impl_std_traits_for_sequence!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);

impl std::fmt::Debug for AminoAcids {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl std::fmt::Debug for AminoAcidsView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl std::fmt::Debug for AminoAcidsViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}
