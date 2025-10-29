//! The definition of the [`KmerEncoder`] trait, along with provided encoders.

use crate::kmer::{Kmer, KmerEncode, SupportedMismatchNumber, errors::KmerError};
use std::hash::Hash;

pub mod three_bit;

/// A trait used to define encoders/decoders for representing bases and k-mers
/// in a more efficienct form.
///
/// While individual bases and k-mers can be manually encoded and decoded, the
/// main purpose of a [`KmerEncoder`] is to automatically provide k-mer
/// algorithms, such as finding a k-mer in a sequence (with a [`KmerSet`]) or
/// counting the k-mers in a sequence (with a [`KmerCounter`]).
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
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is not a valid
    /// length for the encoding.
    fn new(kmer_length: usize) -> Result<Self, KmerError>;

    /// Retrieves the k-mer length associated with this [`KmerEncoder`].
    fn kmer_length(&self) -> usize;

    /// Encodes a single base.
    fn encode_base(base: u8) -> Self::EncodedBase;

    /// Decodes a single base.
    ///
    /// The base must have been generated using this [`KmerEncoder`], otherwise
    /// this function may panic or have unexpected behavior.
    fn decode_base(encoded_base: Self::EncodedBase) -> u8;

    /// Encodes a k-mer.
    ///
    /// The k-mer length is assumed to be valid for the given [`KmerEncoder`].
    /// Consider [`encode_kmer_checked`] when it is not known whether the k-mer
    /// length will be valid.
    ///
    /// [`encode_kmer_checked`]: KmerEncoder::encode_kmer_checked
    fn encode_kmer<S: AsRef<[u8]>>(&self, kmer: S) -> Self::EncodedKmer;

    /// Decodes a k-mer.
    ///
    /// The encoding must have been generated using this [`KmerEncoder`],
    /// otherwise this function may panic or have unexpected behavior.
    fn decode_kmer(&self, encoded_kmer: Self::EncodedKmer) -> Kmer<MAX_LEN>;

    /// Returns an iterator over all encoded k-mers that are at most a Hamming
    /// distance of `N` away from the provided k-mer.
    ///
    /// The original k-mer is included in the iterator. `N` must be a supported
    /// number of mismatches. See [`SupportedMismatchNumber`] for more details.
    fn get_variants<const N: usize>(
        &self, encoded_kmer: Self::EncodedKmer,
    ) -> <Self::MismatchNumber<N> as SupportedMismatchNumber<MAX_LEN, Self>>::MismatchIter
    where
        Self::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, Self>;

    /// Returns an iterator over the encoded overlapping k-mers in a sequence,
    /// from left to right.
    ///
    /// If the sequence is shorter than the k-mer length of the [`KmerEncoder`],
    /// then the iterator will be empty.
    fn iter_from_sequence<'a, S: AsRef<[u8]> + ?Sized>(&self, seq: &'a S) -> Self::SeqIter<'a>;

    /// Returns a consuming iterator over the encoded overlapping k-mers in a
    /// sequence, from left to right.
    ///
    /// If the sequence is shorter than the k-mer length of the [`KmerEncoder`],
    /// then the iterator will be empty.
    ///
    /// This is similar to [`iter_from_sequence`], but the iterator
    /// consumes/stores the sequence (which must be able to be converted to a
    /// [`Vec`] with [`AsRef`]). This is useful when attempting to map an
    /// iterator of sequences to an iterator of k-mers, in which case
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

    /// Returns an iterator over the encoded overlapping k-mers in a sequence,
    /// from right to left.
    ///
    /// If the sequence is shorter than the k-mer length of the [`KmerEncoder`],
    /// then the iterator will be empty.
    fn iter_from_sequence_rev<'a, S: AsRef<[u8]> + ?Sized>(&self, seq: &'a S) -> Self::SeqIterRev<'a>;

    /// Returns a consuming iterator over the encoded overlapping k-mers in a
    /// sequence, from right to left.
    ///
    /// If the sequence is shorter than the k-mer length of the [`KmerEncoder`],
    /// then the iterator will be empty
    ///
    /// This is similar to [`iter_from_sequence_rev`], but the iterator
    /// consumes/stores the sequence (which must be able to be converted to a
    /// [`Vec`] with [`AsRef`]). This is useful when attempting to map an
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

    /// Encodes a k-mer, checking to ensure that its length is correct (and
    /// returning [`None`] otherwise).
    #[inline]
    #[must_use]
    fn encode_kmer_checked<S: AsRef<[u8]>>(&self, kmer: S) -> Option<Self::EncodedKmer> {
        if kmer.as_ref().len() == self.kmer_length() {
            Some(self.encode_kmer(kmer))
        } else {
            None
        }
    }

    /// Given an iterator of encoded kmers, returns an iterator over the decoded
    /// k-mers.
    ///
    /// The k-mers must have been encoded with this [`KmerEncoder`].
    #[inline]
    fn decode_iter(&self, iter: impl Iterator<Item = Self::EncodedKmer>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        iter.map(|encoded_kmer| self.decode_kmer(encoded_kmer))
    }
}
