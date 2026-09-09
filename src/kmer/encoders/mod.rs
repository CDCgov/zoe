//! The definition of the [`KmerEncoder`] trait, along with provided encoders.

use crate::kmer::{GetVariants, Kmer, KmerEncode, errors::KmerError};
use std::hash::Hash;

pub mod three_bit;
pub mod two_bit;

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
/// [`KmerSet`]: crate::kmer::KmerSet
/// [`KmerCounter`]: crate::kmer::KmerCounter
/// [`SupportedKmerLen`]: crate::kmer::SupportedKmerLen
pub trait KmerEncoder<const MAX_LEN: usize>
where
    Self: Sized, {
    /// The type of an encoded k-mer.
    type EncodedKmer: Eq + Hash + Copy + KmerEncode<MAX_LEN, Self>;
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

    /// Creates a new [`KmerEncoder`] with the specified k-mer length.
    ///
    /// ## Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is not a valid
    /// length for the encoding.
    fn new(kmer_length: usize) -> Result<Self, KmerError>;

    /// Retrieves the k-mer length associated with this [`KmerEncoder`].
    fn kmer_length(&self) -> usize;

    /// Encodes a k-mer where the length is known to be correct.
    ///
    /// ## Validity
    ///
    /// The k-mer length is assumed to be valid for the given [`KmerEncoder`],
    /// and a check is not made.
    fn encode_kmer_unchecked(&self, kmer: impl AsRef<[u8]>) -> Self::EncodedKmer;

    /// Encodes a k-mer.
    ///
    /// ## Panic
    ///
    /// Panics if the length of `kmer` does not match the expected length set in
    /// the encoder.
    fn encode_kmer(&self, kmer: impl AsRef<[u8]>) -> Self::EncodedKmer {
        assert_eq!(
            kmer.as_ref().len(),
            self.kmer_length(),
            "The length of the k-mer must agree with that of the encoder"
        );

        // Validity: check made above
        self.encode_kmer_unchecked(kmer)
    }

    /// Decodes a k-mer.
    ///
    /// The encoding must have been generated using this [`KmerEncoder`],
    /// otherwise this function may panic or have unexpected behavior.
    fn decode_kmer(&self, encoded_kmer: Self::EncodedKmer) -> Kmer<MAX_LEN>;

    /// Returns an iterator over all encoded k-mers that are at most a Hamming
    /// distance of `N` away from the provided k-mer.
    ///
    /// The original k-mer is included in the iterator. `N` must be a supported
    /// number of mismatches. See [`GetVariants`] for more details.
    fn get_variants<const N: usize>(&self, encoded_kmer: Self::EncodedKmer) -> Self::Iter
    where
        Self: GetVariants<N, MAX_LEN>, {
        <Self as GetVariants<N, MAX_LEN>>::variants(self, encoded_kmer)
    }

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
    /// then the iterator will be empty.
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
    #[deprecated(
        since = "0.0.33",
        note = "consider using encode_kmer or performing a manual check beforehand"
    )]
    fn encode_kmer_checked<S: AsRef<[u8]>>(&self, kmer: S) -> Option<Self::EncodedKmer> {
        if kmer.as_ref().len() == self.kmer_length() {
            Some(self.encode_kmer(kmer))
        } else {
            None
        }
    }

    /// Given an iterator of encoded k-mers, returns an iterator over the
    /// decoded k-mers.
    ///
    /// The k-mers must have been encoded with this [`KmerEncoder`].
    #[inline]
    fn decode_iter(&self, iter: impl Iterator<Item = Self::EncodedKmer>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        iter.map(|encoded_kmer| self.decode_kmer(encoded_kmer))
    }
}
