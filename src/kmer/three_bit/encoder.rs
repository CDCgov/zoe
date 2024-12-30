use super::len_mappings::{KmerLen, MaxLenToType, SupportedThreeBitKmerLen};
use crate::data::mappings::THREE_BIT_MAPPING;
use crate::data::types::Int;
use crate::kmer::encoder::Kmer;
use crate::{
    data::types::Uint,
    kmer::{encoder::KmerEncoder, errors::KmerError},
};

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
pub struct ThreeBitEncodedKmer<const MAX_LEN: usize>(MaxLenToType<MAX_LEN>)
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen;

impl<const MAX_LEN: usize, T: Uint> From<T> for ThreeBitEncodedKmer<MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen<T = T>,
{
    #[inline]
    fn from(value: T) -> Self {
        Self(value)
    }
}

impl<const MAX_LEN: usize, T: Uint> std::fmt::Display for ThreeBitEncodedKmer<MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen<T = T>,
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
        // TODO: Potentially make unsafe with sealed traits or ASCII check
        f.write_str(std::str::from_utf8(&buffer[start..]).unwrap())
    }
}

/// A [`KmerEncoder`] using three bits to represent each base. This allows for
/// `A`, `C`, `G`, `T`, and `N` to all be represented. This encoder does not
/// preserve case or the distinction between `T` and `U`. `N` is used as a
/// catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerEncoder<const MAX_LEN: usize>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    kmer_length: usize,
    kmer_mask:   MaxLenToType<MAX_LEN>,
}

impl<const MAX_LEN: usize, T: Uint> KmerEncoder<MAX_LEN> for ThreeBitKmerEncoder<MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen<T = T>,
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
    fn get_kmer_length(&self) -> usize {
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
    /// Any input greater than 8 causes a panic, while less than 3 may cause
    /// unexpected behavior.
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
        Kmer::new(self.kmer_length, buffer)
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
            Some(Kmer::new(self.kmer_length, buffer))
        } else {
            None
        }
    }

    /// Get an iterator over all encoded k-mers that are exactly a Hamming
    /// distance of one away from the provided k-mer. This involves replacing
    /// each base with the other bases in `ACGTN`. The original k-mer is not
    /// included in the iterator.
    #[inline]
    fn get_variants_one_mismatch(&self, encoded_kmer: Self::EncodedKmer) -> Self::OneMismatchIter {
        ThreeBitOneMismatchIter::new(encoded_kmer, self.kmer_length)
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
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    current_kmer: ThreeBitEncodedKmer<MAX_LEN>,
    remaining:    std::slice::Iter<'a, u8>,
    kmer_mask:    MaxLenToType<MAX_LEN>,
}

impl<const MAX_LEN: usize, T: Uint> Iterator for ThreeBitKmerIterator<'_, MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen<T = T>,
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
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen
{
}

/// An iterator over the three-bit encoded overlapping k-mers in a sequence,
/// from right to left.
pub struct ThreeBitKmerIteratorRev<'a, const MAX_LEN: usize>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    current_kmer: ThreeBitEncodedKmer<MAX_LEN>,
    remaining:    std::slice::Iter<'a, u8>,
    kmer_length:  usize,
}

impl<const MAX_LEN: usize, T: Uint> Iterator for ThreeBitKmerIteratorRev<'_, MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen<T = T>,
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
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen
{
}

/// An iterator over all three-bit encoded k-mers that are exactly a Hamming
/// distance of one away from a provided k-mer. The original k-mer is not
/// included in the iterator.
pub struct ThreeBitOneMismatchIter<const MAX_LEN: usize>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    encoded_kmer:       ThreeBitEncodedKmer<MAX_LEN>,
    kmer_length:        usize,
    current_kmer:       ThreeBitEncodedKmer<MAX_LEN>,
    current_index:      usize,
    current_base_num:   usize,
    set_mask_third_bit: MaxLenToType<MAX_LEN>,
}

impl<const MAX_LEN: usize> ThreeBitOneMismatchIter<MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    /// Create a new [`ThreeBitOneMismatchIter`] from a provided `encoded_kmer`
    /// and `kmer_length`. For most purposes, use
    /// [`ThreeBitKmerEncoder::get_variants_one_mismatch`] to obtain an
    /// iterator.
    #[inline]
    #[must_use]
    pub fn new(encoded_kmer: ThreeBitEncodedKmer<MAX_LEN>, kmer_length: usize) -> Self {
        Self {
            encoded_kmer,
            kmer_length,
            current_kmer: encoded_kmer,
            current_index: 0,
            current_base_num: 0,
            set_mask_third_bit: MaxLenToType::<MAX_LEN>::from_literal(0b100),
        }
    }
}

impl<const MAX_LEN: usize> Iterator for ThreeBitOneMismatchIter<MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
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

        None
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = 4 * (self.kmer_length - self.current_index) - self.current_base_num;
        (size, Some(size))
    }
}

impl<const MAX_LEN: usize> ExactSizeIterator for ThreeBitOneMismatchIter<MAX_LEN> where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen
{
}
