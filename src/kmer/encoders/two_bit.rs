//! A module providing a k-mer encoder using two bits per base.
//!
//! See [`TwoBitKmerEncoder`] for more details. This also provides the type
//! aliases [`TwoBitKmerSet`] and [`TwoBitKmerCounter`].

use crate::{
    data::mappings::TWO_BIT_MAPPING,
    kmer::{Kmer, KmerCounter, KmerEncoder, KmerError, KmerLen, KmerSet, MaxLenToType, SupportedKmerLen},
    math::{AnyInt, Uint},
};
use std::hash::{Hash, Hasher, RandomState};

/// A type alias for a [`KmerLen`] struct with the [`TwoBitKmerEncoder`] as its
/// encoder.
pub(crate) type TwoBitKmerLen<const MAX_LEN: usize> = KmerLen<MAX_LEN, TwoBitKmerEncoder<MAX_LEN>>;

/// A type alias for [`MaxLenToType`] with the [`TwoBitKmerEncoder`] as the
/// encoder.
pub(crate) type TwoBitMaxLenToType<const MAX_LEN: usize> = MaxLenToType<MAX_LEN, TwoBitKmerEncoder<MAX_LEN>>;

/// A type alias for a [`KmerSet`] with the [`TwoBitKmerEncoder`] as its
/// encoder.
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
pub type TwoBitKmerSet<const MAX_LEN: usize, S = RandomState> = KmerSet<MAX_LEN, TwoBitKmerEncoder<MAX_LEN>, S>;

/// A type alias for a [`KmerCounter`] with the [`TwoBitKmerEncoder`] as its
/// encoder.
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
pub type TwoBitKmerCounter<const MAX_LEN: usize, S = RandomState> = KmerCounter<MAX_LEN, TwoBitKmerEncoder<MAX_LEN>, S>;

/// An encoded k-mer using the [`TwoBitKmerEncoder`].
///
/// In most use cases, encoded k-mers need not be handled directly;
/// [`TwoBitKmerSet`] and [`TwoBitKmerCounter`] provide many methods for
/// accomplishing common tasks. If no suitable methods are present, then
/// handling encoded k-mers directly and later decoding them may be appropriate.
///
/// [`TwoBitKmerSet`]: crate::kmer::encoders::two_bit::TwoBitKmerSet
/// [`TwoBitKmerCounter`]: crate::kmer::encoders::two_bit::TwoBitKmerCounter
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
#[repr(transparent)]
pub struct TwoBitEncodedKmer<const MAX_LEN: usize>(TwoBitMaxLenToType<MAX_LEN>)
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen;

impl<const MAX_LEN: usize, T: Uint> From<T> for TwoBitEncodedKmer<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    #[inline]
    fn from(value: T) -> Self {
        Self(value)
    }
}

/// A [`KmerEncoder`] using two bits to represent each base.
///
/// This allows for `A`, `C`, `G`, and `T` to all be represented. This encoder
/// does not preserve case or the distinction between `T` and `U`. Any other
/// base is substituted with `A`, so it is important to validate any
/// bases/sequences before encoding them.
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
///
/// ## Acknowledgements
///
/// This k-mer encoding is inspired by the encoding found in
/// [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
/// from the [BBMap](https://sourceforge.net/projects/bbmap/) project.
#[derive(Copy, Clone, Debug)]
pub struct TwoBitKmerEncoder<const MAX_LEN: usize>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    /// The length of the k-mers that this encoder can handle, at most
    /// `MAX_LEN`.
    kmer_length: usize,
    /// A bitmask containing `1` only in bits used to store the k-mer, and `0`s
    /// in unused bit positions
    kmer_mask:   TwoBitMaxLenToType<MAX_LEN>,
}

impl<T: Uint, const MAX_LEN: usize> TwoBitKmerEncoder<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    /// Encodes a single base.
    ///
    /// Four unique bases are represented by the encoder: `Aa`, `Cc`, `Gg`,
    /// and `TUtu`. Any other base is converted to `A`.
    #[inline]
    #[must_use]
    pub fn encode_base(base: u8) -> T {
        T::from(TWO_BIT_MAPPING[base])
    }

    /// Decodes a single base.
    ///
    /// The base must have been generated using [`TwoBitKmerEncoder`], otherwise
    /// this function may panic or have unexpected behavior.
    #[inline]
    #[must_use]
    pub fn decode_base(encoded_base: T) -> u8 {
        // Validity: as_usize is valid since encoded_base will be in `0..=3`,
        // even if `T` is `u128`
        b"ACGT"[encoded_base.cast_as::<usize>()]
    }
}

impl<const MAX_LEN: usize> PartialEq for TwoBitKmerEncoder<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        // Validity: kmer_mask is directly generated from kmer_length and hence
        // does not need to be compared
        self.kmer_length == other.kmer_length
    }
}

impl<const MAX_LEN: usize> Eq for TwoBitKmerEncoder<MAX_LEN> where TwoBitKmerLen<MAX_LEN>: SupportedKmerLen {}

impl<const MAX_LEN: usize> Hash for TwoBitKmerEncoder<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Validity: kmer_mask is directly generated from kmer_length and hence
        // does not need to be hashed
        self.kmer_length.hash(state);
    }
}

impl<const MAX_LEN: usize, T: Uint> KmerEncoder<MAX_LEN> for TwoBitKmerEncoder<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    type EncodedKmer = TwoBitEncodedKmer<MAX_LEN>;
    type SeqIter<'a> = TwoBitKmerIterator<'a, MAX_LEN>;
    type SeqIntoIter = TwoBitKmerIntoIterator<MAX_LEN>;
    type SeqIterRev<'a> = TwoBitKmerIteratorRev<'a, MAX_LEN>;
    type SeqIntoIterRev = TwoBitKmerIntoIteratorRev<MAX_LEN>;

    /// Creates a new [`TwoBitKmerEncoder`] with the specified k-mer length.
    ///
    /// ## Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than the specified `MAX_LEN`.
    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    fn new(kmer_length: usize) -> Result<Self, KmerError> {
        if kmer_length <= MAX_LEN && kmer_length > 1 {
            // Cast is valid since kmer_length <= MAX_LEN < 256. Unbounded shl
            // is required since 2 * kmer_length may equal the number of bits in
            // T, in which case wrapping subtraction is needed to convert 0 to a
            // bitmask with all bits set
            let kmer_mask = (T::ONE.unbounded_shl(2 * kmer_length as u32)).wrapping_sub(T::ONE);
            Ok(Self { kmer_length, kmer_mask })
        } else {
            Err(KmerError::InvalidLength)
        }
    }

    #[inline]
    fn kmer_length(&self) -> usize {
        self.kmer_length
    }

    #[inline]
    fn encode_kmer(&self, kmer: impl AsRef<[u8]>) -> Self::EncodedKmer {
        let mut encoded_kmer = T::ZERO;

        for &base in kmer.as_ref() {
            encoded_kmer = (encoded_kmer << 2) | Self::encode_base(base);
        }

        TwoBitEncodedKmer(encoded_kmer)
    }

    fn decode_kmer(&self, mut encoded_kmer: Self::EncodedKmer) -> Kmer<MAX_LEN> {
        let mut buffer = [0; MAX_LEN];
        for i in (0..self.kmer_length).rev() {
            let encoded_base = encoded_kmer.0 & T::from_literal(0b11);
            buffer[i] = Self::decode_base(encoded_base);
            encoded_kmer.0 >>= 2;
        }

        // Safety: The buffer only contains valid ASCII because it is
        // initialized with 0 and TwoBitKmerEncoder::decode_base always returns
        // a char in b"ACGT"
        unsafe { Kmer::new_unchecked(self.kmer_length, buffer) }
    }

    fn iter_from_sequence<'a, S: AsRef<[u8]> + ?Sized>(&self, seq: &'a S) -> Self::SeqIter<'a> {
        TwoBitKmerIterator::new(self, seq.as_ref())
    }

    fn iter_consuming_seq<S>(&self, seq: S) -> Self::SeqIntoIter
    where
        S: Into<Vec<u8>>,
        for<'a> &'a S: AsRef<Vec<u8>>, {
        TwoBitKmerIntoIterator::new(self, seq.into())
    }

    fn iter_from_sequence_rev<'a, S: AsRef<[u8]> + ?Sized>(&self, seq: &'a S) -> Self::SeqIterRev<'a> {
        TwoBitKmerIteratorRev::new(self, seq.as_ref())
    }

    fn iter_consuming_seq_rev<S>(&self, seq: S) -> Self::SeqIntoIterRev
    where
        S: Into<Vec<u8>>,
        for<'a> &'a S: AsRef<Vec<u8>>, {
        TwoBitKmerIntoIteratorRev::new(self, seq.into())
    }
}

/// An iterator over the two-bit encoded overlapping k-mers in a sequence, from
/// left to right.
pub struct TwoBitKmerIterator<'a, const MAX_LEN: usize>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    current_kmer: TwoBitEncodedKmer<MAX_LEN>,
    remaining:    std::slice::Iter<'a, u8>,
    kmer_mask:    TwoBitMaxLenToType<MAX_LEN>,
}

impl<'a, const MAX_LEN: usize> TwoBitKmerIterator<'a, MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    /// Creates a new [`TwoBitKmerIterator`].
    #[inline]
    fn new(encoder: &TwoBitKmerEncoder<MAX_LEN>, seq: &'a [u8]) -> Self {
        if seq.len() < encoder.kmer_length {
            Self::empty()
        } else {
            let index = encoder.kmer_length - 1;
            // Validity: This code does not conform to the assumptions of the kmer
            // API, since it uses `encode_kmer` on a kmer that is one base too
            // short. However, each iteration of the loop adds another base as
            // the first operation, so this is corrected when the iterator
            // starts.
            let pre_first_kmer = encoder.encode_kmer(&seq[..index]);

            TwoBitKmerIterator {
                current_kmer: pre_first_kmer,
                remaining:    seq[index..].iter(),
                kmer_mask:    encoder.kmer_mask,
            }
        }
    }

    /// Creates an empty [`TwoBitKmerIterator`].
    ///
    /// The only field that matters is `remaining`, which is initialized to an
    /// iterator over an empty slice. Any call to [`next`] will return `None` as
    /// a result. The remaining fields are initialized here to 0, since they are
    /// not needed.
    ///
    /// [`next`]: Iterator::next
    #[inline]
    #[must_use]
    fn empty() -> Self {
        TwoBitKmerIterator {
            current_kmer: TwoBitEncodedKmer(TwoBitMaxLenToType::<MAX_LEN>::ZERO),
            remaining:    [].iter(),
            kmer_mask:    TwoBitMaxLenToType::<MAX_LEN>::ZERO,
        }
    }

    /// Helper function to shift a kmer to the left and add a new base.
    ///
    /// This encapsulates shared behavior between [`TwoBitKmerIterator`] and
    /// [`TwoBitKmerIntoIterator`].
    #[inline]
    #[must_use]
    pub(crate) fn shift_and_add_base(
        new_base: u8, kmer: TwoBitEncodedKmer<MAX_LEN>, kmer_mask: TwoBitMaxLenToType<MAX_LEN>,
    ) -> TwoBitEncodedKmer<MAX_LEN> {
        let encoded_base = TwoBitKmerEncoder::encode_base(new_base);
        // Shift the k-mer to the left to make room for the new base
        let shifted_kmer = kmer.0 << 2;
        // Add the new base and remove the leftmost base (to preserve the
        // length of the k-mer)
        ((shifted_kmer | encoded_base) & kmer_mask).into()
    }
}

impl<const MAX_LEN: usize> Iterator for TwoBitKmerIterator<'_, MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    type Item = TwoBitEncodedKmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let base = *self.remaining.next()?;
        self.current_kmer = TwoBitKmerIterator::shift_and_add_base(base, self.current_kmer, self.kmer_mask);
        Some(self.current_kmer)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.remaining.size_hint()
    }
}

impl<const MAX_LEN: usize> ExactSizeIterator for TwoBitKmerIterator<'_, MAX_LEN> where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen
{
}

/// An iterator over the two-bit encoded overlapping k-mers in a sequence, from
/// left to right, which consumes/stores the sequence.
///
/// For a non-consuming version, use [`TwoBitKmerIterator`].
pub struct TwoBitKmerIntoIterator<const MAX_LEN: usize>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    seq:          Vec<u8>,
    index:        usize,
    current_kmer: TwoBitEncodedKmer<MAX_LEN>,
    kmer_mask:    TwoBitMaxLenToType<MAX_LEN>,
}

impl<const MAX_LEN: usize> TwoBitKmerIntoIterator<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    /// Creates a new [`TwoBitKmerIterator`].
    #[inline]
    fn new(encoder: &TwoBitKmerEncoder<MAX_LEN>, seq: Vec<u8>) -> Self {
        if seq.len() < encoder.kmer_length {
            Self::empty()
        } else {
            let index = encoder.kmer_length - 1;
            // Validity: This code does not conform to the assumptions of the
            // kmer API, since it uses `encode_kmer` on a kmer that is one base
            // too short. However, each iteration of the loop adds another base
            // as the first operation, so this is corrected when the iterator
            // starts.
            let pre_first_kmer = encoder.encode_kmer(&seq[..index]);

            TwoBitKmerIntoIterator {
                seq,
                index,
                current_kmer: pre_first_kmer,
                kmer_mask: encoder.kmer_mask,
            }
        }
    }

    /// Creates an empty [`TwoBitKmerIntoIterator`].
    ///
    /// The only field that matters is `seq`, which is initialized to an empty
    /// vector. Any call to [`next`] will return `None` as a result. The
    /// remaining fields are initialized here to 0, since they are not needed.
    ///
    /// [`next`]: Iterator::next
    #[inline]
    #[must_use]
    fn empty() -> Self {
        Self {
            // This is okay, since it does not allocate
            seq:          Vec::new(),
            index:        0,
            current_kmer: TwoBitEncodedKmer(TwoBitMaxLenToType::<MAX_LEN>::ZERO),
            kmer_mask:    TwoBitMaxLenToType::<MAX_LEN>::ZERO,
        }
    }
}

impl<const MAX_LEN: usize> Iterator for TwoBitKmerIntoIterator<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    type Item = TwoBitEncodedKmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Extract the next base to add on the right end of the k-mer
        let base = *self.seq.get(self.index)?;
        self.index += 1;
        self.current_kmer = TwoBitKmerIterator::shift_and_add_base(base, self.current_kmer, self.kmer_mask);
        Some(self.current_kmer)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = self.seq.len() - self.index;
        (size, Some(size))
    }
}

impl<const MAX_LEN: usize> ExactSizeIterator for TwoBitKmerIntoIterator<MAX_LEN> where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen
{
}

/// An iterator over the two-bit encoded overlapping k-mers in a sequence, from
/// right to left.
pub struct TwoBitKmerIteratorRev<'a, const MAX_LEN: usize>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    current_kmer: TwoBitEncodedKmer<MAX_LEN>,
    remaining:    std::slice::Iter<'a, u8>,
    kmer_length:  usize,
}

impl<'a, const MAX_LEN: usize> TwoBitKmerIteratorRev<'a, MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    fn new(encoder: &TwoBitKmerEncoder<MAX_LEN>, seq: &'a [u8]) -> Self {
        if seq.len() < encoder.kmer_length {
            Self::empty()
        } else {
            let index = seq.len() + 1 - encoder.kmer_length;
            // Validity: This code does not conform to the assumptions of the kmer
            // API, since it uses `encode_kmer` on a kmer that is one base too
            // short. However, each iteration of the loop adds another base as
            // the first operation, so this is corrected when the iterator
            // starts.
            let pre_last_kmer = encoder.encode_kmer(&seq[index..]).0 << 2;

            TwoBitKmerIteratorRev {
                current_kmer: pre_last_kmer.into(),
                remaining:    seq[..index].iter(),
                kmer_length:  encoder.kmer_length,
            }
        }
    }

    /// Creates an empty [`TwoBitKmerIteratorRev`].
    ///
    /// The only field that matters is `remaining`, which is initialized to an
    /// iterator over an empty slice. Any call to [`next`] will return `None` as
    /// a result. The remaining fields are initialized here to 0, since they are
    /// not needed.
    ///
    /// [`next`]: Iterator::next
    #[inline]
    #[must_use]
    fn empty() -> Self {
        Self {
            current_kmer: TwoBitEncodedKmer(TwoBitMaxLenToType::<MAX_LEN>::ZERO),
            remaining:    [].iter(),
            kmer_length:  0,
        }
    }

    /// Helper function to shift a kmer to the right and add a new base.
    ///
    /// This encapsulates shared behavior between [`TwoBitKmerIteratorRev`] and
    /// [`TwoBitKmerIntoIteratorRev`].
    #[inline]
    #[must_use]
    pub(crate) fn shift_and_add_base_rev(
        new_base: u8, kmer: TwoBitEncodedKmer<MAX_LEN>, kmer_length: usize,
    ) -> TwoBitEncodedKmer<MAX_LEN> {
        let encoded_base = TwoBitKmerEncoder::encode_base(new_base);
        // Shift the new base to the proper position
        let shifted_encoded_base = encoded_base << (2 * (kmer_length - 1));
        // Shift the k-mer to the right to make room for the new base and remove
        // the rightmost base, then add the new base
        ((kmer.0 >> 2) | shifted_encoded_base).into()
    }
}

impl<const MAX_LEN: usize> Iterator for TwoBitKmerIteratorRev<'_, MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    type Item = TwoBitEncodedKmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Extract the next base to add on the left end of the k-mer
        let base = *self.remaining.next_back()?;
        self.current_kmer = TwoBitKmerIteratorRev::shift_and_add_base_rev(base, self.current_kmer, self.kmer_length);
        Some(self.current_kmer)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.remaining.size_hint()
    }
}

impl<const MAX_LEN: usize> ExactSizeIterator for TwoBitKmerIteratorRev<'_, MAX_LEN> where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen
{
}

/// An iterator over the two-bit encoded overlapping k-mers in a sequence, from
/// right to left, which consumes/stores the sequence.
///
/// For a non-consuming version, use [`TwoBitKmerIteratorRev`].
pub struct TwoBitKmerIntoIteratorRev<const MAX_LEN: usize>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    seq:          Vec<u8>,
    current_kmer: TwoBitEncodedKmer<MAX_LEN>,
    kmer_length:  usize,
}

impl<const MAX_LEN: usize> TwoBitKmerIntoIteratorRev<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    fn new(encoder: &TwoBitKmerEncoder<MAX_LEN>, mut seq: Vec<u8>) -> Self {
        if seq.len() < encoder.kmer_length {
            Self::empty()
        } else {
            let index = seq.len() + 1 - encoder.kmer_length;
            // Validity: This code does not conform to the assumptions of the kmer
            // API, since it uses `encode_kmer` on a kmer that is one base too
            // short. However, each iteration of the loop adds another base as
            // the first operation, so this is corrected when the iterator
            // starts.
            let pre_last_kmer = encoder.encode_kmer(&seq[index..]).0 << 2;

            seq.truncate(index);

            TwoBitKmerIntoIteratorRev {
                seq,
                current_kmer: pre_last_kmer.into(),
                kmer_length: encoder.kmer_length,
            }
        }
    }

    /// Creates an empty [`TwoBitKmerIntoIteratorRev`].
    ///
    /// The only field that matters is `index`, which is initialized to 0. Any
    /// call to [`next`] will return `None` as a result. The remaining fields
    /// are initialized here to 0 (and an empty vector for `seq`), since they
    /// are not needed.
    ///
    /// [`next`]: Iterator::next
    #[inline]
    #[must_use]
    fn empty() -> Self {
        Self {
            // This is okay since it does not allocate
            seq:          Vec::new(),
            current_kmer: TwoBitEncodedKmer(TwoBitMaxLenToType::<MAX_LEN>::ZERO),
            kmer_length:  0,
        }
    }
}

impl<const MAX_LEN: usize> Iterator for TwoBitKmerIntoIteratorRev<MAX_LEN>
where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    type Item = TwoBitEncodedKmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Extract the next base to add on the left end of the k-mer
        let base = self.seq.pop()?;
        self.current_kmer = TwoBitKmerIteratorRev::shift_and_add_base_rev(base, self.current_kmer, self.kmer_length);
        Some(self.current_kmer)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = self.seq.len();
        (size, Some(size))
    }
}

impl<const MAX_LEN: usize> ExactSizeIterator for TwoBitKmerIntoIteratorRev<MAX_LEN> where
    TwoBitKmerLen<MAX_LEN>: SupportedKmerLen
{
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_decoding() {
        assert_eq!(TwoBitKmerEncoder::<21>::decode_base(0), b'A');
        assert_eq!(TwoBitKmerEncoder::<21>::decode_base(1), b'C');
        assert_eq!(TwoBitKmerEncoder::<21>::decode_base(2), b'G');
        assert_eq!(TwoBitKmerEncoder::<21>::decode_base(3), b'T');
    }

    #[test]
    #[should_panic(expected = "index out of bounds: the len is 4 but the index is 4")]
    fn test_invalid_decoding() {
        let _ = TwoBitKmerEncoder::<21>::decode_base(4);
    }
}
