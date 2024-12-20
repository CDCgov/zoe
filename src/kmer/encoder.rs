use crate::kmer::errors::KmerError;

// TODO: Derives
pub struct Kmer<const MAX_LEN: usize> {
    pub(crate) length: u8,
    pub(crate) buffer: [u8; MAX_LEN],
}

impl<const MAX_LEN: usize> Kmer<MAX_LEN> {
    #[allow(clippy::cast_possible_truncation)]
    #[inline]
    #[must_use]
    pub fn new(length: usize, buffer: [u8; MAX_LEN]) -> Self {
        Self {
            // This will not truncate since we require all kmers to be of length
            // less than 256
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
}

impl<const MAX_LEN: usize> std::fmt::Display for Kmer<MAX_LEN> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // TODO: Unless we seal the KmerEncoder trait, we need to avoid unsafe here
        // Or we would need a separate Kmer struct for each encoder
        // Should we avoid unwrap?
        f.write_str(std::str::from_utf8(self.as_ref()).unwrap())
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
/// [`KmerSet`]: super::kmer_set::KmerSet
/// [`KmerCounter`]: super::kmer_counter::KmerCounter
pub trait KmerEncoder<const MAX_LEN: usize>
where
    Self: Sized, {
    /// The type of an encoded base.
    type EncodedBase;
    /// The type of an encoded k-mer.
    type EncodedKmer: From<Self::EncodedBase>;
    /// An iterator over the encoded overlapping k-mers in a sequence, from left
    /// to right.
    type SeqIter<'a>: Iterator<Item = Self::EncodedKmer>;
    /// An iterator over the encoded overlapping k-mers in a sequence, from right
    /// to left.
    type SeqIterRev<'a>: Iterator<Item = Self::EncodedKmer>;
    /// An iterator over all encoded k-mers that are exactly a Hamming distance
    /// of one away from a provided k-mer. The original k-mer is not included in
    /// the iterator.
    type OneMismatchIter: Iterator<Item = Self::EncodedKmer>;

    /// Creates a new [`KmerEncoder`] with the specified k-mer length.
    ///
    /// # Errors
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is not a valid
    /// length for the encoding.
    fn new(kmer_length: usize) -> Result<Self, KmerError>;

    /// Retrieve the k-mer length associated with this [`KmerEncoder`].
    fn get_kmer_length(&self) -> usize;

    /// Encode a single base. The base is assumed to be valid for the given
    /// [`KmerEncoder`]. Consider [`encode_base_checked`] when it is not known
    /// whether the base will be valid.
    ///
    /// [`encode_base_checked`]: KmerEncoder::encode_base_checked
    fn encode_base(base: u8) -> Self::EncodedBase;

    /// Encode a single base. If the base is not valid for the given
    /// [`KmerEncoder`], then `None` is returned.
    fn encode_base_checked(base: u8) -> Option<Self::EncodedBase>;

    /// Decode a single base. The base is assumed to be valid for the given
    /// [`KmerEncoder`]. Consider [`decode_base_checked`] when it is not known
    /// whether the base will be valid.
    ///
    /// [`decode_base_checked`]: KmerEncoder::decode_base_checked
    fn decode_base(encoded_base: Self::EncodedBase) -> u8;

    /// Decode a single base. If the base is not valid for the given
    /// [`KmerEncoder`], then `None` is returned.
    fn decode_base_checked(encoded_base: Self::EncodedBase) -> Option<u8>;

    /// Encode a k-mer. The bases and k-mer length are assumed to be valid for
    /// the given [`KmerEncoder`]. Consider [`encode_kmer_checked`] when it is
    /// not known whether the bases and k-mer length will be valid.
    ///
    /// [`encode_kmer_checked`]: KmerEncoder::encode_kmer_checked
    fn encode_kmer<S: AsRef<[u8]>>(&self, kmer: S) -> Self::EncodedKmer;

    /// Encode a k-mer. If the bases and k-mer length are not valid for the
    /// given [`KmerEncoder`], then `None` is returned.
    fn encode_kmer_checked<S: AsRef<[u8]>>(&self, kmer: S) -> Option<Self::EncodedKmer>;

    /// Decode a k-mer. The bases and k-mer length are assumed to be valid for
    /// the given [`KmerEncoder`]. Consider [`decode_kmer_checked`] when it is
    /// not known whether the bases and k-mer length will be valid.
    ///
    /// [`decode_kmer_checked`]: KmerEncoder::decode_kmer_checked
    fn decode_kmer(&self, encoded_kmer: Self::EncodedKmer) -> Kmer<MAX_LEN>;

    /// Decode a k-mer. If the bases and k-mer length are not valid for the
    /// given [`KmerEncoder`], then `None` is returned.
    fn decode_kmer_checked(&self, encoded_kmer: Self::EncodedKmer) -> Option<Kmer<MAX_LEN>>;

    /// Get an iterator over all encoded k-mers that are exactly a Hamming
    /// distance of one away from the provided k-mer. The original k-mer is not
    /// included in the iterator.
    fn get_variants_one_mismatch(&self, encoded_kmer: Self::EncodedKmer) -> Self::OneMismatchIter;

    /// Get an iterator over the encoded overlapping k-mers in a sequence, from
    /// left to right. If the sequence is shorter than the k-mer length of the
    /// [`KmerEncoder`], then the iterator will be empty. The sequence should
    /// only contain valid bases for the given [`KmerEncoder`].
    fn iter_from_sequence<'a, S: AsRef<[u8]>>(&self, seq: &'a S) -> Self::SeqIter<'a>;

    /// Get an iterator over the encoded overlapping k-mers in a sequence, from
    /// right to left. If the sequence is shorter than the k-mer length of the
    /// [`KmerEncoder`], then the iterator will be empty. The sequence should
    /// only contain valid bases for the given [`KmerEncoder`].
    fn iter_from_sequence_rev<'a, S: AsRef<[u8]>>(&self, seq: &'a S) -> Self::SeqIterRev<'a>;
}
