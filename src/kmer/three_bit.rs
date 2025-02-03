use crate::{
    data::mappings::THREE_BIT_MAPPING,
    kmer::{Kmer, KmerEncoder, KmerError, KmerLen, KmerSet, MaxLenToType, SupportedKmerLen},
    math::Uint,
    prelude::KmerCounter,
};
use std::hash::{Hash, Hasher, RandomState};

use super::{MismatchNumber, SupportedMismatchNumber};

/// A type alias for a [`KmerLen`] struct with the [`ThreeBitKmerEncoder`] as
/// its encoder.
pub(crate) type ThreeBitKmerLen<const MAX_LEN: usize> = KmerLen<MAX_LEN, ThreeBitKmerEncoder<MAX_LEN>>;

/// A type alias for [`MaxLenToType`] with the [`ThreeBitKmerEncoder`] as the
/// encoder.
pub(crate) type ThreeBitMaxLenToType<const MAX_LEN: usize> = MaxLenToType<MAX_LEN, ThreeBitKmerEncoder<MAX_LEN>>;

/// A type alias for a [`KmerSet`] with the [`ThreeBitKmerEncoder`] as its
/// encoder.
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
pub type ThreeBitKmerSet<const MAX_LEN: usize, S = RandomState> = KmerSet<MAX_LEN, ThreeBitKmerEncoder<MAX_LEN>, S>;

/// A type alias for a [`KmerCounter`] with the [`ThreeBitKmerEncoder`] as its
/// encoder.
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
pub type ThreeBitKmerCounter<const MAX_LEN: usize, S = RandomState> = KmerCounter<MAX_LEN, ThreeBitKmerEncoder<MAX_LEN>, S>;

/// An encoded k-mer using the [`ThreeBitKmerEncoder`]. In most use cases,
/// encoded k-mers need not be handled directly; [`ThreeBitKmerSet`] and
/// [`ThreeBitKmerCounter`] provide many methods for accomplishing common tasks.
/// If no suitable methods are present, then handling encoded k-mers directly
/// and later decoding them may be appropriate.
///
/// To decode a [`ThreeBitKmerEncoder`], call the `decode_kmer` method of the
/// encoder. For the potential case where the encoder is not available,
/// `Display` is implemented for [`ThreeBitEncodedKmer`], which will output the
/// decoded k-mer. However, it is more performant to decode first.
///
/// Note that the ordering given by [`Ord`] and [`PartialOrd`] is different
/// between [`ThreeBitEncodedKmer`] and the corresponding decoded [`Kmer`].
///
/// [`ThreeBitKmerSet`]: super::ThreeBitKmerSet
/// [`ThreeBitKmerCounter`]: super::ThreeBitKmerCounter
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
#[repr(transparent)]
pub struct ThreeBitEncodedKmer<const MAX_LEN: usize>(ThreeBitMaxLenToType<MAX_LEN>)
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen;

impl<const MAX_LEN: usize, T: Uint> From<T> for ThreeBitEncodedKmer<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    #[inline]
    fn from(value: T) -> Self {
        Self(value)
    }
}

impl<const MAX_LEN: usize, T: Uint> std::fmt::Display for ThreeBitEncodedKmer<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut kmer = self.0;
        let mut buffer = [0; MAX_LEN];
        let mut start = MAX_LEN;
        while kmer != T::ZERO && start > 0 {
            start -= 1;
            let encoded_base = kmer & T::from_literal(0b111);
            buffer[start] = ThreeBitKmerEncoder::decode_base(encoded_base);
            kmer >>= 3;
        }
        // Safety: The buffer only contains valid ASCII because it is
        // initialized with 0 and ThreeBitKmerEncoder::decode_base always
        // returns a char in b"000NACGT"
        f.write_str(unsafe { std::str::from_utf8_unchecked(&buffer[start..]) })
    }
}

/// A [`KmerEncoder`] using three bits to represent each base. This allows for
/// `A`, `C`, `G`, `T`, and `N` to all be represented. This encoder does not
/// preserve case or the distinction between `T` and `U`. `N` is used as a
/// catch-all for bases that are not `ACGTUNacgtun`.
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
///
/// ### Acknowledgements
///
/// This k-mer encoding is inspired by the 2-bit encoding found in
/// [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
/// from the [BBMap](https://sourceforge.net/projects/bbmap/) project.
///
/// [`get_variants_one_mismatch`]: KmerEncoder::get_variants_one_mismatch
/// [`encode_kmer`]: KmerEncoder::encode_kmer
/// [`decode_kmer`]: KmerEncoder::decode_kmer
/// [`iter_from_sequence`]: KmerEncoder::iter_from_sequence
#[derive(Clone, Debug)]
pub struct ThreeBitKmerEncoder<const MAX_LEN: usize>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    kmer_length: usize,
    kmer_mask:   ThreeBitMaxLenToType<MAX_LEN>,
}

impl<const MAX_LEN: usize> PartialEq for ThreeBitKmerEncoder<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.kmer_length == other.kmer_length
    }
}

impl<const MAX_LEN: usize> Eq for ThreeBitKmerEncoder<MAX_LEN> where ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen {}

impl<const MAX_LEN: usize> Hash for ThreeBitKmerEncoder<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.kmer_length.hash(state);
    }
}

impl<const MAX_LEN: usize, T: Uint> KmerEncoder<MAX_LEN> for ThreeBitKmerEncoder<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    type EncodedBase = T;
    type EncodedKmer = ThreeBitEncodedKmer<MAX_LEN>;
    type SeqIter<'a> = ThreeBitKmerIterator<'a, MAX_LEN>;
    type SeqIterRev<'a> = ThreeBitKmerIteratorRev<'a, MAX_LEN>;
    type OneMismatchIter = ThreeBitOneMismatchIter<MAX_LEN>;

    /// Creates a new [`ThreeBitKmerEncoder`] with the specified k-mer length.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than the specified `MAX_LEN`.
    #[inline]
    fn new(kmer_length: usize) -> Result<Self, KmerError> {
        if kmer_length <= MAX_LEN && kmer_length > 1 {
            Ok(Self {
                kmer_length,
                kmer_mask: (T::ONE << (3 * kmer_length)) - T::ONE,
            })
        } else {
            Err(KmerError::InvalidLength)
        }
    }

    /// Retrieve the k-mer length associated with this [`ThreeBitKmerEncoder`].
    #[inline]
    fn kmer_length(&self) -> usize {
        self.kmer_length
    }

    /// Encode a single base. The encoder can handle any u8 input without panic.
    /// Five unique bases are represented by the encoder: `Aa`, `Cc`, `Gg`,
    /// `TUtu`, or `N` as a catch-all for any other input.
    #[inline]
    fn encode_base(base: u8) -> Self::EncodedBase {
        T::from(THREE_BIT_MAPPING[base])
    }

    /// The [`ThreeBitKmerEncoder`] can handle all `u8` inputs, so this function
    /// is equivalent to [`encode_base`]. It will never return `None`.
    ///
    /// [`encode_base`]: ThreeBitKmerEncoder::encode_base
    #[inline]
    fn encode_base_checked(base: u8) -> Option<Self::EncodedBase> {
        Some(Self::encode_base(base))
    }

    /// Decode a single base. The base must be in `3..=8` (an invariant upheld
    /// by [`ThreeBitKmerEncoder`]), and may panic or have unexpected behavior
    /// otherwise. Consider [`decode_base_checked`] when it is not known whether
    /// the base will be valid.
    ///
    /// # Panics
    ///
    /// The [`ThreeBitKmerEncoder`] represents each base by a `u8` in `3..=8`.
    /// Other inputs may panic or cause unexpected behavior.
    ///
    /// [`decode_base_checked`]: ThreeBitKmerEncoder::decode_base_checked
    #[inline]
    fn decode_base(encoded_base: Self::EncodedBase) -> u8 {
        // as_usize is valid since encoded_base will be in `3..=8`
        b"000NACGT"[encoded_base.as_usize()]
    }

    /// Decode a single base. If the base is not in `3..=8`, then `None` is
    /// returned.
    #[inline]
    fn decode_base_checked(encoded_base: Self::EncodedBase) -> Option<u8> {
        if (T::from(3)..=T::from(7)).contains(&encoded_base) {
            Some(Self::decode_base(encoded_base))
        } else {
            None
        }
    }

    /// Encode a k-mer. The length of `kmer` must match the `kmer_length` of the
    /// [`ThreeBitKmerEncoder`], otherwise unexpected results can occur.
    /// Consider [`encode_kmer_checked`] when it is not known whether `kmer` has
    /// the appropriate length.
    ///
    /// [`encode_kmer_checked`]: ThreeBitKmerEncoder::encode_kmer_checked
    #[inline]
    fn encode_kmer<S: AsRef<[u8]>>(&self, kmer: S) -> Self::EncodedKmer {
        let mut encoded_kmer = T::ZERO;

        for &base in kmer.as_ref() {
            encoded_kmer = (encoded_kmer << 3) | Self::encode_base(base);
        }

        ThreeBitEncodedKmer(encoded_kmer)
    }

    /// Encode a k-mer. If the length of `kmer` does not match the `kmer_length`
    /// of the [`ThreeBitKmerEncoder`], then `None` is returned.
    #[inline]
    fn encode_kmer_checked<S: AsRef<[u8]>>(&self, kmer: S) -> Option<Self::EncodedKmer> {
        if kmer.as_ref().len() == self.kmer_length {
            Some(self.encode_kmer(kmer))
        } else {
            None
        }
    }

    /// Decode a k-mer. The bases and k-mer length are assumed to be valid for
    /// the given [`ThreeBitKmerEncoder`]. If an invalid base is encountered (a
    /// three-bit value outside `3..=8`), then this function may panic or have
    /// unexpected behavior. If the length of the encoded k-mer does not agree
    /// with the `kmer_length` of the [`ThreeBitKmerEncoder`], then either some
    /// bases may be truncated, or erroneous extra bases may appear in the
    /// output. Consider [`decode_kmer_checked`] when it is not known whether
    /// the bases and k-mer length will be valid.
    ///
    /// # Panics
    ///
    /// If an invalid base is encountered (a three-bit value outside `3..=8`),
    /// then this function will panic or cause unexpected behavior.
    ///
    /// [`decode_kmer_checked`]: ThreeBitKmerEncoder::decode_kmer_checked
    fn decode_kmer(&self, mut encoded_kmer: Self::EncodedKmer) -> Kmer<MAX_LEN> {
        let mut buffer = [0; MAX_LEN];
        for i in (0..self.kmer_length).rev() {
            let encoded_base = encoded_kmer.0 & T::from_literal(0b111);
            buffer[i] = Self::decode_base(encoded_base);
            encoded_kmer.0 >>= 3;
        }
        // Safety: The buffer only contains valid ASCII because it is
        // initialized with 0 and ThreeBitKmerEncoder::decode_base always
        // returns a char in b"000NACGT"
        unsafe { Kmer::new_unchecked(self.kmer_length, buffer) }
    }

    /// Decode a k-mer. If an erroneous base outside `3..=8` is found, or if the
    /// encoding represents a longer k-mer than the `kmer_length` of the
    /// [`ThreeBitKmerEncoder`], then `None` is returned.
    fn decode_kmer_checked(&self, mut encoded_kmer: Self::EncodedKmer) -> Option<Kmer<MAX_LEN>> {
        let mut buffer = [0; MAX_LEN];
        for i in (0..self.kmer_length).rev() {
            let encoded_base = encoded_kmer.0 & T::from_literal(0b111);
            buffer[i] = Self::decode_base_checked(encoded_base)?;
            encoded_kmer.0 >>= 3;
        }
        if encoded_kmer.0 == T::ZERO {
            // Safety: The buffer only contains valid ASCII because it is
            // initialized with 0 and ThreeBitKmerEncoder::decode_base always
            // returns a char in b"000NACGT"
            Some(unsafe { Kmer::new_unchecked(self.kmer_length, buffer) })
        } else {
            None
        }
    }

    /// Get an iterator over all encoded k-mers that are at most a Hamming
    /// distance of N away from the provided k-mer. This involves replacing each
    /// base with the other bases in `ACGTN`. The original k-mer is included in
    /// the iterator.
    ///
    /// `N` must be a supported number of mismatches. See
    /// [`SupportedMismatchNumber`] for more details.
    #[inline]
    fn get_variants<const N: usize>(
        &self, encoded_kmer: Self::EncodedKmer,
    ) -> <MismatchNumber<N> as SupportedMismatchNumber<MAX_LEN, Self>>::MismatchIter
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, Self>, {
        MismatchNumber::get_iterator(encoded_kmer, self)
    }

    /// Get an iterator over the encoded overlapping k-mers in a sequence, from
    /// left to right. If the sequence is shorter than the `kmer_length` of the
    /// [`ThreeBitKmerEncoder`], then the iterator will be empty.
    #[inline]
    fn iter_from_sequence<'a, S: AsRef<[u8]> + ?Sized>(&self, seq: &'a S) -> Self::SeqIter<'a> {
        let seq = seq.as_ref();
        if seq.len() < self.kmer_length {
            ThreeBitKmerIterator {
                current_kmer: ThreeBitEncodedKmer(T::ZERO),
                remaining:    [].iter(),
                kmer_mask:    self.kmer_mask,
            }
        } else {
            let index = self.kmer_length - 1;
            let pre_first_kmer = self.encode_kmer(&seq[..index]);

            ThreeBitKmerIterator {
                current_kmer: pre_first_kmer,
                remaining:    seq[index..].iter(),
                kmer_mask:    self.kmer_mask,
            }
        }
    }

    /// Get an iterator over the encoded overlapping k-mers in a sequence, from
    /// right to left. If the sequence is shorter than the `kmer_length` of the
    /// [`ThreeBitKmerEncoder`], then the iterator will be empty.
    #[inline]
    fn iter_from_sequence_rev<'a, S: AsRef<[u8]> + ?Sized>(&self, seq: &'a S) -> Self::SeqIterRev<'a> {
        let seq = seq.as_ref();
        if seq.len() < self.kmer_length {
            ThreeBitKmerIteratorRev {
                current_kmer: T::ZERO.into(),
                remaining:    [].iter(),
                kmer_length:  self.kmer_length,
            }
        } else {
            let index = seq.len() + 1 - self.kmer_length;
            let pre_last_kmer = self.encode_kmer(&seq[index..]).0 << 3;

            ThreeBitKmerIteratorRev {
                current_kmer: pre_last_kmer.into(),
                remaining:    seq[..index].iter(),
                kmer_length:  self.kmer_length,
            }
        }
    }
}

/// An iterator over the three-bit encoded overlapping k-mers in a sequence,
/// from left to right.
pub struct ThreeBitKmerIterator<'a, const MAX_LEN: usize>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    current_kmer: ThreeBitEncodedKmer<MAX_LEN>,
    remaining:    std::slice::Iter<'a, u8>,
    kmer_mask:    ThreeBitMaxLenToType<MAX_LEN>,
}

impl<const MAX_LEN: usize> Iterator for ThreeBitKmerIterator<'_, MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    type Item = ThreeBitEncodedKmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Extract the next base to add on the right end of the k-mer
        let base = *self.remaining.next()?;
        let encoded_base = ThreeBitKmerEncoder::encode_base(base);
        // Shift the k-mer to the left to make room for the new base
        let shifted_kmer = self.current_kmer.0 << 3;
        // Add the new base and remove the leftmost base (to preserve the
        // length of the k-mer)
        self.current_kmer = ((shifted_kmer | encoded_base) & self.kmer_mask).into();
        Some(self.current_kmer)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.remaining.size_hint()
    }
}

impl<const MAX_LEN: usize> ExactSizeIterator for ThreeBitKmerIterator<'_, MAX_LEN> where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen
{
}

/// An iterator over the three-bit encoded overlapping k-mers in a sequence,
/// from right to left.
pub struct ThreeBitKmerIteratorRev<'a, const MAX_LEN: usize>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    current_kmer: ThreeBitEncodedKmer<MAX_LEN>,
    remaining:    std::slice::Iter<'a, u8>,
    kmer_length:  usize,
}

impl<const MAX_LEN: usize> Iterator for ThreeBitKmerIteratorRev<'_, MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    type Item = ThreeBitEncodedKmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Extract the next base to add on the left end of the k-mer
        let base = *self.remaining.next_back()?;
        let encoded_base = ThreeBitKmerEncoder::<MAX_LEN>::encode_base(base);
        // Shift the new base to the proper position
        let shifted_encoded_base = encoded_base << (3 * (self.kmer_length - 1));
        // Shift the k-mer to the right to make room for the new base and remove
        // the rightmost base, then add the new base
        self.current_kmer = ((self.current_kmer.0 >> 3) | shifted_encoded_base).into();
        Some(self.current_kmer)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.remaining.size_hint()
    }
}

impl<const MAX_LEN: usize> ExactSizeIterator for ThreeBitKmerIteratorRev<'_, MAX_LEN> where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen
{
}

/// An iterator over all three-bit encoded k-mers that are at most a Hamming
/// distance of 1 away from a provided k-mer. The original k-mer is included in
/// the iterator.
pub struct ThreeBitOneMismatchIter<const MAX_LEN: usize>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    encoded_kmer:       ThreeBitEncodedKmer<MAX_LEN>,
    kmer_length:        usize,
    current_kmer:       ThreeBitEncodedKmer<MAX_LEN>,
    current_index:      usize,
    current_base_num:   usize,
    set_mask_third_bit: ThreeBitMaxLenToType<MAX_LEN>,
    not_finished:       bool,
}

impl<const MAX_LEN: usize, T: Uint> ThreeBitOneMismatchIter<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    /// Create a new [`ThreeBitOneMismatchIter`] from a provided `encoded_kmer`
    /// and `kmer_length`. For most purposes, use
    /// [`ThreeBitKmerEncoder::get_variants_one_mismatch`] to obtain an
    /// iterator.
    #[inline]
    #[must_use]
    pub(crate) fn new(encoded_kmer: ThreeBitEncodedKmer<MAX_LEN>, kmer_encoder: &ThreeBitKmerEncoder<MAX_LEN>) -> Self {
        Self {
            encoded_kmer,
            kmer_length: kmer_encoder.kmer_length(),
            current_kmer: encoded_kmer,
            current_index: 0,
            current_base_num: 0,
            set_mask_third_bit: T::from_literal(0b100),
            not_finished: true,
        }
    }
}

impl<const MAX_LEN: usize> Iterator for ThreeBitOneMismatchIter<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    type Item = ThreeBitEncodedKmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Check for end of iterator
        while self.current_index < self.kmer_length {
            // Mutate 4 times to other bases
            if self.current_base_num < 4 {
                self.current_base_num += 1;
                // When current_kmer = 100, then the second summand will be 0
                // and the new value will be 000 due to the first summand
                // getting masked. When current_kmer = 0**, then the first
                // summand is current_kmer and the second summand is 001.
                self.current_kmer = ((self.current_kmer.0 | self.set_mask_third_bit)
                    - ((self.current_kmer.0 & self.set_mask_third_bit) >> 2))
                    .into();
                return Some(self.current_kmer);
            }

            self.current_kmer = self.encoded_kmer;
            self.current_index += 1;
            self.current_base_num = 0;
            self.set_mask_third_bit <<= 3;
        }

        if self.not_finished {
            self.not_finished = false;
            return Some(self.encoded_kmer);
        }
        None
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = usize::from(self.not_finished) + 4 * (self.kmer_length - self.current_index) - self.current_base_num;
        (size, Some(size))
    }
}

impl<const MAX_LEN: usize> ExactSizeIterator for ThreeBitOneMismatchIter<MAX_LEN> where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen
{
}

/// An iterator over all three-bit encoded k-mers that are at most a Hamming
/// distance of `N` away from a provided k-mer, where `N >= 2`. The original
/// k-mer is included in the iterator.
///
/// ### Acknowledgements
///
/// The code for subset iteration was inspired by [this
/// blog](https://fishi.devtail.io/weblog/2015/05/18/common-bitwise-techniques-subset-iterations/)
pub struct ThreeBitMismatchIter<const MAX_LEN: usize, const N: usize>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen, {
    /// The original encoded k-mer
    encoded_kmer:            ThreeBitMaxLenToType<MAX_LEN>,
    /// The current encoded k-mer, which is mutated as the iterator cycles
    /// through variants
    current_kmer:            ThreeBitMaxLenToType<MAX_LEN>,
    /// The current subset of bases being mutated, encoded as a bitmask where
    /// ith bit from the left is 1 precisely when the ith base in the kmer is
    /// being mutated
    current_subset:          ThreeBitMaxLenToType<MAX_LEN>,
    /// The number of bases being mutated in `current_subset`
    subset_size:             usize,
    /// A buffer to hold the `set_mask_third_bit` values for each index in
    /// `current_subset`, as well as the number of times each base in the subset
    /// was mutated. This value will always be at least 1, since every base in
    /// the subset must be mutated
    masks_and_times_mutated: [(ThreeBitMaxLenToType<MAX_LEN>, usize); N],
    /// The upper bound for `current_subset`. We must have `current_subset`
    /// strictly less than `max_subset`
    max_subset:              ThreeBitMaxLenToType<MAX_LEN>,
}

impl<const MAX_LEN: usize, const N: usize, T: Uint> ThreeBitMismatchIter<MAX_LEN, N>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    #[inline]
    #[must_use]
    pub(crate) fn new(encoded_kmer: ThreeBitEncodedKmer<MAX_LEN>, encoder: &ThreeBitKmerEncoder<MAX_LEN>) -> Self {
        // This implementation assumes N >= 2
        const { assert!(N >= 2) }

        // Due to the value in the initialization of times_mutated being
        // usize::MAX, this will trigger the creation of a new subset (encoded
        // subset 1) on the first call to next. We could modify the
        // initialization here to avoid this, but weirdly it results in a
        // worsening of performance.
        Self {
            encoded_kmer:            encoded_kmer.0,
            current_kmer:            encoded_kmer.0,
            current_subset:          T::ZERO,
            subset_size:             0,
            masks_and_times_mutated: [(T::from_literal(0b100), usize::MAX); N],
            max_subset:              T::ONE << encoder.kmer_length(),
        }
    }
}

impl<const MAX_LEN: usize, const N: usize, T: Uint> Iterator for ThreeBitMismatchIter<MAX_LEN, N>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen<T = T>,
{
    type Item = ThreeBitEncodedKmer<MAX_LEN>;

    fn next(&mut self) -> Option<Self::Item> {
        // Save the output: this iterator yields then updates
        let out = ThreeBitEncodedKmer(self.current_kmer);

        if self.masks_and_times_mutated[0].1 < 4 {
            // This first branch is redundant with the next, but increases
            // efficiency. Mutate the first base in the subset.

            self.current_kmer = (self.current_kmer | self.masks_and_times_mutated[0].0)
                - ((self.current_kmer & self.masks_and_times_mutated[0].0) >> 2);
            self.masks_and_times_mutated[0].1 += 1;
        } else if let Some(i) = self.masks_and_times_mutated[..self.subset_size]
            .iter()
            .position(|(_, base_num)| *base_num < 4)
        {
            // This second branch means that a base was found which hasn't been
            // modified 4 times yet (namely the leftmost).

            // Reset each base before it
            for (set_mask_third_bit, base_num) in &mut self.masks_and_times_mutated[..i] {
                // Modify twice in order to cycle to original base, then first
                // variant base. This may seem inefficient, but it is cheaper
                // than storing a stack of encoded kmers.
                let set_mask_third_bit = *set_mask_third_bit;
                self.current_kmer =
                    (self.current_kmer | set_mask_third_bit) - ((self.current_kmer & set_mask_third_bit) >> 2);
                self.current_kmer =
                    (self.current_kmer | set_mask_third_bit) - ((self.current_kmer & set_mask_third_bit) >> 2);
                *base_num = 1;
            }

            // Mutate the identified base once
            let set_mask_third_bit = self.masks_and_times_mutated[i].0;
            self.current_kmer = (self.current_kmer | set_mask_third_bit) - ((self.current_kmer & set_mask_third_bit) >> 2);
            self.masks_and_times_mutated[i].1 += 1;
        } else {
            // This final branch means we need a new subset. In either case of
            // the below code, the encoded subset (represented as T) will get
            // larger.

            if (self.current_subset.count_ones() as usize) < N {
                // Our subset has not reached maximum size, so we can add 1
                self.current_subset += T::ONE;
            } else {
                // Replace all trailing 0s with 1s, then add 1. This always
                // decreases the number of 1s
                self.current_subset = (self.current_subset | (self.current_subset - T::ONE)) + T::ONE;
            }

            // Check if we are at the end of the iterator
            if self.current_subset >= self.max_subset {
                // The iterator has just reached its end, but we may still need
                // to return Some(out). The subsets are visited in order, so
                // since N >= 1, we will visit self.max_mask (it has only a
                // single element). Hence, when we first arrive at this point in
                // the code, we will have self.current_subset equal to
                // self.max_mask. If this is the case, the below variable
                // becomes Some(out), otherwise it is None.
                let out = (self.current_subset == self.max_subset).then_some(out);

                // Adjust current_subset to make it larger than max_mask, so
                // that we return None from now on. This will not cause
                // overflow, since we are only utilizing one third of the bits
                // in T.
                self.current_subset = self.max_subset + T::ONE;
                return out;
            }

            // Reset current_kmer
            self.current_kmer = self.encoded_kmer;

            // Reprocess the saved state for this new subset
            self.subset_size = 0;
            let mut processed_bases = 0;
            let mut temp_subset = self.current_subset;
            while temp_subset != T::ZERO {
                let num_zeros = temp_subset.trailing_zeros() as usize;
                temp_subset >>= num_zeros + 1;
                processed_bases += num_zeros;
                self.masks_and_times_mutated[self.subset_size].0 = T::from_literal(0b100) << (processed_bases * 3);
                processed_bases += 1;
                self.subset_size += 1;
            }

            // Precompute first variant: must modify each base in the subset
            for (set_mask_third_bit, base_num) in &mut self.masks_and_times_mutated[..self.subset_size] {
                self.current_kmer =
                    (self.current_kmer | *set_mask_third_bit) - ((self.current_kmer & *set_mask_third_bit) >> 2);
                *base_num = 1;
            }
        }

        Some(out)
    }

    // We could implement size_hint and ExactSizeIterator, but it is expensive
    // to compute the size of the iterator after it has already started
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_decoding() {
        assert_eq!(ThreeBitKmerEncoder::<21>::decode_base(0), b'0');
        assert_eq!(ThreeBitKmerEncoder::<21>::decode_base(1), b'0');
        assert_eq!(ThreeBitKmerEncoder::<21>::decode_base(2), b'0');
        assert_eq!(ThreeBitKmerEncoder::<21>::decode_base(3), b'N');
        assert_eq!(ThreeBitKmerEncoder::<21>::decode_base(4), b'A');
        assert_eq!(ThreeBitKmerEncoder::<21>::decode_base(5), b'C');
        assert_eq!(ThreeBitKmerEncoder::<21>::decode_base(6), b'G');
        assert_eq!(ThreeBitKmerEncoder::<21>::decode_base(7), b'T');
    }

    #[test]
    #[should_panic(expected = "index out of bounds: the len is 8 but the index is 8")]
    fn test_invalid_decoding() {
        ThreeBitKmerEncoder::<21>::decode_base(8);
    }
}
