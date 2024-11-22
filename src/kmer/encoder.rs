use crate::kmer::errors::KmerError;

/// [`KmerEncoder`] represents an encoder/decoder used to represent bases and
/// kmers in a more efficienct form. While individual bases and kmers can be
/// manually encoded and decoded, the main purpose of a [`KmerEncoder`] is to
/// mesh with kmer algorithm such as finding a kmer in a sequence (as is
/// possible with a [`KmerSet`]) or counting the kmers in a sequence (as is
/// possible with a [`KmerCounter`]).
///
/// [`KmerSet`]: super::kmer_set::KmerSet
/// [`KmerCounter`]: super::kmer_counter::KmerCounter
pub trait KmerEncoder
where
    Self: Sized, {
    /// The number of bits used to represent each base.
    const BITS_PER_BASE: usize;
    /// The maximum possible kmer length representable with this
    /// [`KmerEncoder`].
    const MAX_KMER_LENGTH: usize;
    /// The type of an encoded base.
    type Base;
    /// The type of an encoded kmer.
    type Kmer;
    /// An iterator over the encoded overlapping kmers in a sequence, from left
    /// to right.
    type SeqIter<'a>: Iterator<Item = Self::Kmer>;
    /// An iterator over the encoded overlapping kmers in a sequence, from right
    /// to left.
    type SeqIterRev<'a>: Iterator<Item = Self::Kmer>;
    /// An iterator over all encoded kmers that are exactly a Hamming distance
    /// of one away from a provided kmer. The original kmer is not included in
    /// the iterator.
    type OneMismatchIter: Iterator<Item = Self::Kmer>;

    /// Creates a new [`KmerEncoder`] with the specified kmer length.
    ///
    /// # Errors
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is not a valid
    /// length for the encoding.
    fn new(kmer_length: usize) -> Result<Self, KmerError>;

    /// Retrieve the kmer length associated with this [`KmerEncoder`].
    fn get_kmer_length(&self) -> usize;

    /// Encode a single base. The base is assumed to be valid for the given
    /// [`KmerEncoder`]. Consider [`encode_base_checked`] when it is not known
    /// whether the base will be valid.
    ///
    /// [`encode_base_checked`]: KmerEncoder::encode_base_checked
    fn encode_base(base: u8) -> Self::Base;

    /// Encode a single base. If the base is not valid for the given
    /// [`KmerEncoder`], then `None` is returned.
    fn encode_base_checked(base: u8) -> Option<Self::Base>;

    /// Decode a single base. The base is assumed to be valid for the given
    /// [`KmerEncoder`]. Consider [`decode_base_checked`] when it is not known
    /// whether the base will be valid.
    ///
    /// [`decode_base_checked`]: KmerEncoder::decode_base_checked
    fn decode_base(encoded_base: Self::Base) -> u8;

    /// Decode a single base. If the base is not valid for the given
    /// [`KmerEncoder`], then `None` is returned.
    fn decode_base_checked(encoded_base: Self::Base) -> Option<u8>;

    /// Encode a kmer. The bases and kmer length are assumed to be valid for the
    /// given [`KmerEncoder`]. Consider [`encode_kmer_checked`] when it is not
    /// known whether the bases and kmer length will be valid.
    ///
    /// [`encode_kmer_checked`]: KmerEncoder::encode_kmer_checked
    fn encode_kmer(&self, kmer: &[u8]) -> Self::Kmer;

    /// Encode a kmer. If the bases and kmer length are not valid for the given
    /// [`KmerEncoder`], then `None` is returned.
    fn encode_kmer_checked(&self, kmer: &[u8]) -> Option<Self::Kmer>;

    /// Decode a kmer. The bases and kmer length are assumed to be valid for the
    /// given [`KmerEncoder`]. Consider [`decode_kmer_checked`] when it is not
    /// known whether the bases and kmer length will be valid.
    ///
    /// [`decode_kmer_checked`]: KmerEncoder::decode_kmer_checked
    fn decode_kmer(&self, encoded_kmer: Self::Kmer) -> Vec<u8>;

    /// Decode a kmer. If the bases and kmer length are not valid for the given
    /// [`KmerEncoder`], then `None` is returned.
    fn decode_kmer_checked(&self, encoded_kmer: Self::Kmer) -> Option<Vec<u8>>;

    /// Get an iterator over all encoded kmers that are exactly a Hamming
    /// distance of one away from the provided kmer. The original kmer is not
    /// included in the iterator.
    fn get_variants_one_mismatch(&self, encoded_kmer: Self::Kmer) -> Self::OneMismatchIter;

    /// Get an iterator over the encoded overlapping kmers in a sequence, from
    /// left to right. If the sequence is shorter than the kmer length of the
    /// [`KmerEncoder`], then the iterator will be empty. The sequence should
    /// only contain valid bases for the given [`KmerEncoder`].
    fn iter_from_sequence<'a>(&self, seq: &'a [u8]) -> Self::SeqIter<'a>;

    /// Get an iterator over the encoded overlapping kmers in a sequence, from
    /// right to left. If the sequence is shorter than the kmer length of the
    /// [`KmerEncoder`], then the iterator will be empty. The sequence should
    /// only contain valid bases for the given [`KmerEncoder`].
    fn iter_from_sequence_rev<'a>(&self, seq: &'a [u8]) -> Self::SeqIterRev<'a>;
}
