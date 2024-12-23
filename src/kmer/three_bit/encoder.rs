use crate::{
    data::DNA_PROFILE_MAP,
    kmer::{encoder::KmerEncoder, errors::KmerError},
};

/// A [`KmerEncoder`] using three bits to represent each base. This allows for
/// `A`, `C`, `G`, `T`, and `N` to all be represented. This encoder does not
/// preserve case or the distinction between `T` and `U`. `N` is used as a
/// catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerEncoder {
    kmer_length: usize,
    kmer_mask:   EncodedKmer,
}

#[derive(PartialEq, Eq, Hash, Clone, Copy)]
#[repr(transparent)]
pub struct EncodedKmer(u64);

impl From<u64> for EncodedKmer {
    #[inline]
    fn from(value: u64) -> Self {
        EncodedKmer(value)
    }
}

impl std::fmt::Display for EncodedKmer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut kmer = self.0;
        let mut buffer = [0; ThreeBitKmerEncoder::MAX_KMER_LENGTH];
        let mut start = ThreeBitKmerEncoder::MAX_KMER_LENGTH;
        while kmer != 0 && start > 0 {
            start -= 1;
            let encoded_base = kmer & 0b111;
            buffer[start] = ThreeBitKmerEncoder::decode_base(encoded_base as u8);
            kmer >>= ThreeBitKmerEncoder::BITS_PER_BASE;
        }
        f.write_str(&String::from_utf8_lossy(&buffer))
    }
}

impl KmerEncoder for ThreeBitKmerEncoder {
    const BITS_PER_BASE: usize = 3;
    const MAX_KMER_LENGTH: usize = 21;
    type Base = u8;
    type Kmer = EncodedKmer;
    type SeqIter<'a> = ThreeBitKmerIterator<'a>;
    type SeqIterRev<'a> = ThreeBitKmerIteratorRev<'a>;
    type OneMismatchIter = ThreeBitOneMismatchIter;

    /// Creates a new [`ThreeBitKmerEncoder`] with the specified k-mer length.
    ///
    /// # Errors
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is 0 or greater
    /// than 21.
    #[inline]
    fn new(kmer_length: usize) -> Result<Self, KmerError> {
        if kmer_length <= Self::MAX_KMER_LENGTH && kmer_length > 0 {
            Ok(Self {
                kmer_length,
                kmer_mask: ((1u64 << (Self::BITS_PER_BASE * kmer_length)) - 1).into(),
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
    fn encode_base(base: u8) -> u8 {
        DNA_PROFILE_MAP[base] + 1
    }

    /// The [`ThreeBitKmerEncoder`] can handle all `u8` inputs, so this function
    /// is equivalent to [`encode_base`]. It will never return None.
    ///
    /// [`encode_base`]: ThreeBitKmerEncoder::encode_base
    #[inline]
    fn encode_base_checked(base: u8) -> Option<u8> {
        Some(Self::encode_base(base))
    }

    /// Decode a single base. The base must be in `0..=4` (an invariant upheld
    /// by [`ThreeBitKmerEncoder`]), and will panic otherwise. Consider
    /// [`decode_base_checked`] when it is not known whether the base will be
    /// valid.
    ///
    /// # Panics
    ///
    /// The [`ThreeBitKmerEncoder`] represents each base by a `u8` in `0..=4`.
    /// Any input of this range will panic.
    ///
    /// [`decode_base_checked`]: ThreeBitKmerEncoder::decode_base_checked
    #[inline]
    fn decode_base(encoded_base: u8) -> u8 {
        b"ACGTN"[(encoded_base - 1) as usize]
    }

    /// Decode a single base. If the base is not in `0..=4`, then `None` is
    /// returned.
    #[inline]
    fn decode_base_checked(encoded_base: u8) -> Option<u8> {
        b"ACGTN".get((encoded_base - 1) as usize).copied()
    }

    /// Encode a k-mer. The length of `kmer` must match the `kmer_length` of the
    /// [`ThreeBitKmerEncoder`], otherwise unexpected results can occur.
    /// Consider [`encode_kmer_checked`] when it is not known whether `kmer` has
    /// the appropriate length.
    ///
    /// [`encode_kmer_checked`]: ThreeBitKmerEncoder::encode_kmer_checked
    #[inline]
    fn encode_kmer(&self, kmer: &[u8]) -> EncodedKmer {
        let mut encoded_kmer = 0;

        for &base in kmer {
            encoded_kmer = (encoded_kmer << Self::BITS_PER_BASE) | u64::from(Self::encode_base(base));
        }

        EncodedKmer(encoded_kmer)
    }

    /// Encode a k-mer. If the length of `kmer` does not match the `kmer_length`
    /// of the [`ThreeBitKmerEncoder`], then `None` is returned.
    #[inline]
    fn encode_kmer_checked(&self, kmer: &[u8]) -> Option<EncodedKmer> {
        if kmer.len() == self.kmer_length {
            Some(self.encode_kmer(kmer))
        } else {
            None
        }
    }

    /// Decode a k-mer. The bases and k-mer length are assumed to be valid for the
    /// given [`ThreeBitKmerEncoder`]. If an invalid base is encountered (a
    /// three-bit value outside `0..=4`), then this function will panic. If the
    /// length of the encoded k-mer does not agree with the `kmer_length` of the
    /// [`ThreeBitKmerEncoder`], then either some bases may be truncated, or
    /// erroneous extra bases may appear in the output. Consider
    /// [`decode_kmer_checked`] when it is not known whether the bases and k-mer
    /// length will be valid.
    ///
    /// # Panics
    ///
    /// If an invalid base is encountered (a three-bit value outside `0..=4`),
    /// then this function will panic.
    ///
    /// [`decode_kmer_checked`]: ThreeBitKmerEncoder::decode_kmer_checked
    fn decode_kmer(&self, mut encoded_kmer: EncodedKmer) -> Vec<u8> {
        let mut kmer = vec![0; self.kmer_length];
        for i in (0..kmer.len()).rev() {
            let encoded_base = encoded_kmer.0 & 0b111;
            kmer[i] = Self::decode_base(encoded_base as u8);
            encoded_kmer.0 >>= Self::BITS_PER_BASE;
        }
        kmer
    }

    /// Decode a k-mer. If an erroneous base outside `0..=4` is found, or if the
    /// encoding represents a longer k-mer than the `kmer_length` of the
    /// [`ThreeBitKmerEncoder`], then `None` is returned.
    fn decode_kmer_checked(&self, mut encoded_kmer: EncodedKmer) -> Option<Vec<u8>> {
        let mut kmer = vec![0; self.kmer_length];
        for i in (0..kmer.len()).rev() {
            let encoded_base = encoded_kmer.0 & 0b111;
            kmer[i] = Self::decode_base_checked(encoded_base as u8)?;
            encoded_kmer.0 >>= Self::BITS_PER_BASE;
        }
        if encoded_kmer.0 == 0 {
            Some(kmer)
        } else {
            None
        }
    }

    /// Get an iterator over all encoded k-mers that are exactly a Hamming
    /// distance of one away from the provided k-mer. This involves replacing
    /// each base with the other bases in `ACGTN`. The original k-mer is not
    /// included in the iterator.
    #[inline]
    fn get_variants_one_mismatch(&self, encoded_kmer: EncodedKmer) -> ThreeBitOneMismatchIter {
        ThreeBitOneMismatchIter::new(encoded_kmer, self.kmer_length)
    }

    /// Get an iterator over the encoded overlapping k-mers in a sequence, from
    /// left to right. If the sequence is shorter than the `kmer_length` of the
    /// [`ThreeBitKmerEncoder`], then the iterator will be empty.
    #[inline]
    fn iter_from_sequence<'a>(&self, seq: &'a [u8]) -> ThreeBitKmerIterator<'a> {
        if seq.len() < self.kmer_length {
            ThreeBitKmerIterator {
                current_kmer: EncodedKmer(0),
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
    fn iter_from_sequence_rev<'a>(&self, seq: &'a [u8]) -> ThreeBitKmerIteratorRev<'a> {
        if seq.len() < self.kmer_length {
            ThreeBitKmerIteratorRev {
                current_kmer: 0.into(),
                remaining:    [].iter(),
                kmer_length:  self.kmer_length,
            }
        } else {
            let index = seq.len() + 1 - self.kmer_length;
            let pre_last_kmer = self.encode_kmer(&seq[index..]).0 << Self::BITS_PER_BASE;

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
pub struct ThreeBitKmerIterator<'a> {
    current_kmer: EncodedKmer,
    remaining:    std::slice::Iter<'a, u8>,
    kmer_mask:    EncodedKmer,
}

impl Iterator for ThreeBitKmerIterator<'_> {
    type Item = EncodedKmer;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Extract the next base to add on the right end of the k-mer
        let base = *self.remaining.next()?;
        let encoded_base = u64::from(ThreeBitKmerEncoder::encode_base(base));
        // Shift the k-mer to the left to make room for the new base
        let shifted_kmer = self.current_kmer.0 << ThreeBitKmerEncoder::BITS_PER_BASE;
        // Add the new base and remove the leftmost base (to preserve the
        // length of the k-mer)
        self.current_kmer = ((shifted_kmer | encoded_base) & self.kmer_mask.0).into();
        Some(self.current_kmer)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.remaining.size_hint()
    }
}

impl ExactSizeIterator for ThreeBitKmerIterator<'_> {}

/// An iterator over the three-bit encoded overlapping k-mers in a sequence,
/// from right to left.
pub struct ThreeBitKmerIteratorRev<'a> {
    current_kmer: EncodedKmer,
    remaining:    std::slice::Iter<'a, u8>,
    kmer_length:  usize,
}

impl Iterator for ThreeBitKmerIteratorRev<'_> {
    type Item = EncodedKmer;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // Extract the next base to add on the left end of the k-mer
        let base = *self.remaining.next_back()?;
        let encoded_base = u64::from(ThreeBitKmerEncoder::encode_base(base));
        // Shift the new base to the proper position
        let shifted_encoded_base = encoded_base << (3 * (self.kmer_length - 1));
        // Shift the k-mer to the right to make room for the new base and remove
        // the rightmost base, then add the new base
        self.current_kmer = ((self.current_kmer.0 >> ThreeBitKmerEncoder::BITS_PER_BASE) | shifted_encoded_base).into();
        Some(self.current_kmer)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.remaining.size_hint()
    }
}

impl ExactSizeIterator for ThreeBitKmerIteratorRev<'_> {}

/// An iterator over all three-bit encoded k-mers that are exactly a Hamming
/// distance of one away from a provided k-mer. The original k-mer is not
/// included in the iterator.
pub struct ThreeBitOneMismatchIter {
    encoded_kmer:       EncodedKmer,
    kmer_length:        usize,
    current_kmer:       EncodedKmer,
    current_index:      usize,
    current_base_num:   usize,
    set_mask_third_bit: u64,
}

impl ThreeBitOneMismatchIter {
    /// Create a new [`ThreeBitOneMismatchIter`] from a provided `encoded_kmer`
    /// and `kmer_length`. For most purposes, use
    /// [`ThreeBitKmerEncoder::get_variants_one_mismatch`] to obtain an
    /// iterator.
    #[inline]
    #[must_use]
    pub fn new(encoded_kmer: EncodedKmer, kmer_length: usize) -> Self {
        Self {
            encoded_kmer,
            kmer_length,
            current_kmer: encoded_kmer,
            current_index: 0,
            current_base_num: 0,
            set_mask_third_bit: 0b100,
        }
    }
}

impl Iterator for ThreeBitOneMismatchIter {
    type Item = EncodedKmer;

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
                self.current_kmer = ((self.current_kmer.0 & !self.set_mask_third_bit)
                    + ((!self.current_kmer.0 & self.set_mask_third_bit) >> 2))
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

impl ExactSizeIterator for ThreeBitOneMismatchIter {}
