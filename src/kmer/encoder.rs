use std::hash::Hash;

use crate::{data::CheckSequence, kmer::errors::KmerError, prelude::Nucleotides};

use super::{KmerEncode, SupportedMismatchNumber};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Kmer<const MAX_LEN: usize> {
    pub(crate) length: u8,
    pub(crate) buffer: [u8; MAX_LEN],
}

impl<const MAX_LEN: usize> Ord for Kmer<MAX_LEN> {
    #[inline]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.as_ref().cmp(other.as_ref())
    }
}

impl<const MAX_LEN: usize> PartialOrd for Kmer<MAX_LEN> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const MAX_LEN: usize> std::fmt::Display for Kmer<MAX_LEN> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Safety: A Kmer instance is guaranteed to be valid ASCII
        f.write_str(unsafe { std::str::from_utf8_unchecked(self.as_ref()) })
    }
}

impl<const MAX_LEN: usize> Kmer<MAX_LEN> {
    /// Create a new `Kmer` instance.
    ///
    /// ## Panics
    ///
    /// Panics if the buffer does not contain valid ASCII, or if the length of
    /// `bases` is less than 2 or longer than `MAX_LEN`.
    #[inline]
    #[must_use]
    pub fn new<S: AsRef<[u8]>>(bases: S) -> Self {
        let bases = bases.as_ref();
        let length = bases.len();

        assert!(bases.is_ascii_simd::<16>());
        assert!(bases.len() >= 2);
        let mut buffer = [0; MAX_LEN];
        buffer[..length].copy_from_slice(bases);

        // Safety: buffer is initialized with 0 (valid ASCII) and is filled with
        // values from bases, which was verified to be valid ASCII
        unsafe { Kmer::new_unchecked(length, buffer) }
    }

    /// Create a new `Kmer` instance from a length and buffer.
    ///
    /// ## Safety
    ///
    /// `buffer` must be valid ASCII, and length must be between 2 and
    /// `MAX_LEN`.
    #[allow(clippy::cast_possible_truncation)]
    #[inline]
    #[must_use]
    pub unsafe fn new_unchecked(length: usize, buffer: [u8; MAX_LEN]) -> Self {
        Self {
            // This will not truncate since we require all k-mers to be of
            // length less than 256
            length: length as u8,
            buffer,
        }
    }

    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.length as usize
    }

    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    #[inline]
    #[must_use]
    pub fn to_vec(&self) -> Vec<u8> {
        self.as_ref().to_vec()
    }

    #[inline]
    #[must_use]
    pub fn to_nucleotides(&self) -> Nucleotides {
        Nucleotides(self.to_vec())
    }
}

impl<const MAX_LEN: usize> AsRef<[u8]> for Kmer<MAX_LEN> {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.buffer[..self.len()]
    }
}

/// [`KmerEncoder`] represents an encoder/decoder used to represent bases and
/// k-mers in a more efficienct form. While individual bases and k-mers can be
/// manually encoded and decoded, the main purpose of a [`KmerEncoder`] is to
/// mesh with k-mer algorithm such as finding a k-mer in a sequence (as is
/// possible with a [`KmerSet`]) or counting the k-mers in a sequence (as is
/// possible with a [`KmerCounter`]).
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
///
/// [`KmerSet`]: super::kmer_set::KmerSet
/// [`KmerCounter`]: super::kmer_counter::KmerCounter
/// [`SupportedKmerLen`]: super::SupportedKmerLen
pub trait KmerEncoder<const MAX_LEN: usize>
where
    Self: Sized, {
    /// The type of an encoded base.
    type EncodedBase;
    /// The type of an encoded k-mer.
    type EncodedKmer: From<Self::EncodedBase> + Eq + Hash + Copy + KmerEncode<MAX_LEN, Self>;
    /// An iterator over the encoded overlapping k-mers in a sequence, from left
    /// to right.
    type SeqIter<'a>: Iterator<Item = Self::EncodedKmer>;
    /// A consuming iterator over the encoded overlapping k-mers in a sequence,
    /// from left to right.
    type SeqIntoIter: Iterator<Item = Self::EncodedKmer>;
    /// An iterator over the encoded overlapping k-mers in a sequence, from
    /// right to left.
    type SeqIterRev<'a>: Iterator<Item = Self::EncodedKmer>;
    /// A consumuing iterator over the encoded overlapping k-mers in a sequence,
    /// from right to left.
    type SeqIntoIterRev: Iterator<Item = Self::EncodedKmer>;
    /// An iterator over all encoded k-mers that are at most a Hamming distance
    /// of one away from a provided k-mer. The original k-mer is included in the
    /// iterator.
    type OneMismatchIter: Iterator<Item = Self::EncodedKmer>;
    /// A zero-size struct used to hold a number of mismatches. For valid values
    /// of `N`, [`SupportedMismatchNumber`] will be implemented, and will
    /// provide the appropriate iterator to use for generating variants. See
    /// [`SupportedMismatchNumber`] for more details.
    type MismatchNumber<const N: usize>;

    /// Creates a new [`KmerEncoder`] with the specified k-mer length.
    ///
    /// ## Errors
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is not a valid
    /// length for the encoding.
    fn new(kmer_length: usize) -> Result<Self, KmerError>;

    /// Retrieves the k-mer length associated with this [`KmerEncoder`].
    fn kmer_length(&self) -> usize;

    /// Encodes a single base. The base is assumed to be valid for the given
    /// [`KmerEncoder`]. Consider [`encode_base_checked`] when it is not known
    /// whether the base will be valid.
    ///
    /// [`encode_base_checked`]: KmerEncoder::encode_base_checked
    fn encode_base(base: u8) -> Self::EncodedBase;

    /// Encodes a single base. If the base is not valid for the given
    /// [`KmerEncoder`], then `None` is returned.
    fn encode_base_checked(base: u8) -> Option<Self::EncodedBase>;

    /// Decodes a single base. The base is assumed to be valid for the given
    /// [`KmerEncoder`]. Consider [`decode_base_checked`] when it is not known
    /// whether the base will be valid.
    ///
    /// [`decode_base_checked`]: KmerEncoder::decode_base_checked
    fn decode_base(encoded_base: Self::EncodedBase) -> u8;

    /// Decodes a single base. If the base is not valid for the given
    /// [`KmerEncoder`], then `None` is returned.
    fn decode_base_checked(encoded_base: Self::EncodedBase) -> Option<u8>;

    /// Encodes a k-mer. The bases and k-mer length are assumed to be valid for
    /// the given [`KmerEncoder`]. Consider [`encode_kmer_checked`] when it is
    /// not known whether the bases and k-mer length will be valid.
    ///
    /// [`encode_kmer_checked`]: KmerEncoder::encode_kmer_checked
    fn encode_kmer<S: AsRef<[u8]>>(&self, kmer: S) -> Self::EncodedKmer;

    /// Encodes a k-mer. If the bases and k-mer length are not valid for the
    /// given [`KmerEncoder`], then `None` is returned.
    fn encode_kmer_checked<S: AsRef<[u8]>>(&self, kmer: S) -> Option<Self::EncodedKmer>;

    /// Decodes a k-mer. The bases and k-mer length are assumed to be valid for
    /// the given [`KmerEncoder`]. Consider [`decode_kmer_checked`] when it is
    /// not known whether the bases and k-mer length will be valid.
    ///
    /// [`decode_kmer_checked`]: KmerEncoder::decode_kmer_checked
    fn decode_kmer(&self, encoded_kmer: Self::EncodedKmer) -> Kmer<MAX_LEN>;

    /// Decodes a k-mer. If the bases and k-mer length are not valid for the
    /// given [`KmerEncoder`], then `None` is returned.
    fn decode_kmer_checked(&self, encoded_kmer: Self::EncodedKmer) -> Option<Kmer<MAX_LEN>>;

    /// Gets an iterator over all encoded k-mers that are at most a Hamming
    /// distance of `N` away from the provided k-mer. The original k-mer is
    /// included in the iterator.
    ///
    /// `N` must be a supported number of mismatches. See
    /// [`SupportedMismatchNumber`] for more details.
    fn get_variants<const N: usize>(
        &self, encoded_kmer: Self::EncodedKmer,
    ) -> <Self::MismatchNumber<N> as SupportedMismatchNumber<MAX_LEN, Self>>::MismatchIter
    where
        Self::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, Self>;

    /// Gets an iterator over the encoded overlapping k-mers in a sequence, from
    /// left to right. If the sequence is shorter than the k-mer length of the
    /// [`KmerEncoder`], then the iterator will be empty. The sequence should
    /// only contain valid bases for the given [`KmerEncoder`].
    fn iter_from_sequence<'a, S: AsRef<[u8]> + ?Sized>(&self, seq: &'a S) -> Self::SeqIter<'a>;

    /// Gets a consuming iterator over the encoded overlapping k-mers in a
    /// sequence, from left to right. If the sequence is shorter than the k-mer
    /// length of the [`KmerEncoder`], then the iterator will be empty. The
    /// sequence should only contain valid bases for the given [`KmerEncoder`].
    ///
    /// This is similar to [`iter_from_sequence`], but the iterator
    /// consumes/stores the sequence. This is useful when attempting to map an
    /// iterator of sequences to an iterator of kmers, in which case
    /// [`iter_from_sequence`] would not work because the sequence would be
    /// dropped too soon.
    ///
    /// [`iter_from_sequence`]: KmerEncoder::iter_from_sequence
    fn iter_consuming_seq<S>(&self, seq: S) -> Self::SeqIntoIter
    where
        S: Into<Vec<u8>>,
        for<'a> &'a S: AsRef<Vec<u8>>;
    // Note: the HRTB is not necessary, but it prevents this function from
    // accepting slices or other types which would incur an expensive
    // conversion. In such cases, `iter_from_sequence` should be used instead.

    /// Gets an iterator over the encoded overlapping k-mers in a sequence, from
    /// right to left. If the sequence is shorter than the k-mer length of the
    /// [`KmerEncoder`], then the iterator will be empty. The sequence should
    /// only contain valid bases for the given [`KmerEncoder`].
    fn iter_from_sequence_rev<'a, S: AsRef<[u8]> + ?Sized>(&self, seq: &'a S) -> Self::SeqIterRev<'a>;

    /// Gets a consuming iterator over the encoded overlapping k-mers in a
    /// sequence, from right to left. If the sequence is shorter than the k-mer
    /// length of the [`KmerEncoder`], then the iterator will be empty. The
    /// sequence should only contain valid bases for the given [`KmerEncoder`].
    ///
    /// This is similar to [`iter_from_sequence_rev`], but the iterator
    /// consumes/stores the sequence. This is useful when attempting to map an
    /// iterator of sequences to an iterator of kmers, in which case
    /// [`iter_from_sequence_rev`] would not work because the sequence would be
    /// dropped too soon.
    ///
    /// [`iter_from_sequence_rev`]: KmerEncoder::iter_from_sequence_rev
    fn iter_consuming_seq_rev<S>(&self, seq: S) -> Self::SeqIntoIterRev
    where
        S: Into<Vec<u8>>,
        for<'a> &'a S: AsRef<Vec<u8>>;
    // Note: the HRTB is not necessary, but it prevents this function from
    // accepting slices or other types which would incur an expensive
    // conversion. In such cases, `iter_from_sequence_rev` should be used
    // instead.
}
