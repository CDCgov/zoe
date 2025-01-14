use super::{Kmer, KmerError, KmerLen, SupportedKmerLen};
use crate::kmer::encoder::KmerEncoder;
use std::{
    collections::{HashMap, hash_map},
    hash::{BuildHasher, RandomState},
    ops::{Index, Range},
};

pub struct KmerCounter<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S = RandomState> {
    map:         HashMap<E::EncodedKmer, usize, S>,
    pub encoder: E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> KmerCounter<MAX_LEN, E> {
    /// Creates a new [`KmerCounter`] with the specified k-mer length.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn new(kmer_length: usize) -> Result<Self, KmerError> {
        Ok(Self {
            map:     HashMap::default(),
            encoder: E::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> KmerCounter<MAX_LEN, E, S> {
    /// Creates a new [`KmerCounter`] with the specified k-mer length
    /// and hasher.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn with_hasher(kmer_length: usize, hasher: S) -> Result<Self, KmerError> {
        Ok(Self {
            map:     HashMap::with_hasher(hasher),
            encoder: E::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> Index<E::EncodedKmer> for KmerCounter<MAX_LEN, E, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    type Output = usize;

    #[inline]
    fn index(&self, index: E::EncodedKmer) -> &Self::Output {
        &self.map[&index]
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerCounter<MAX_LEN, E, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    /// Get the length of the k-mers being stored in the counter.
    #[inline]
    pub fn get_kmer_length(&self) -> usize {
        self.encoder.get_kmer_length()
    }

    /// If the already encoded k-mer is present in this counter, then increment
    /// its count. Otherwise, add it to the counter with a count of 1. The
    /// encoded k-mer must have been generated using the [`KmerEncoder`]
    /// associated with this [`KmerCounter`].
    #[inline]
    pub fn insert_encoded_kmer(&mut self, encoded_kmer: E::EncodedKmer) {
        *self.map.entry(encoded_kmer).or_default() += 1;
    }

    /// If the k-mer is present in this counter, then increment its count.
    /// Otherwise, add it to the counter with a count of 1. The bases and k-mer
    /// length are assumed to be valid for the [`KmerEncoder`] associated with
    /// this [`KmerCounter`]. Consider [`insert_kmer_checked`] when it is not known
    /// whether the bases and k-mer length will be valid.
    ///
    /// [`insert_kmer_checked`]: KmerCounter::insert_kmer_checked
    #[inline]
    pub fn insert_kmer<Q: AsRef<[u8]>>(&mut self, kmer: Q) {
        self.insert_encoded_kmer(self.encoder.encode_kmer(kmer));
    }

    /// If the k-mer is present in this counter, then increment its count.
    /// Otherwise, add it to the counter with a count of 1. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerCounter`], then `false` is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_checked<Q: AsRef<[u8]>>(&mut self, kmer: Q) -> bool {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer(encoded_kmer);
        true
    }

    /// Return whether an already encoded k-mer is present in this counter
    /// (i.e., has a count of at least 1). The encoded k-mer must have been
    /// generated using the [`KmerEncoder`] associated with this
    /// [`KmerCounter`].
    #[inline]
    pub fn contains_encoded(&self, encoded_kmer: E::EncodedKmer) -> bool {
        self.map.contains_key(&encoded_kmer)
    }

    /// Return whether a k-mer is present in this counter (i.e., has a count of
    /// at least 1). The bases and k-mer length are assumed to be valid for the
    /// [`KmerEncoder`] associated with this [`KmerCounter`]. Consider
    /// [`contains_checked`] when it is not known whether the bases and k-mer
    /// length will be valid.
    ///
    /// [`contains_checked`]: KmerCounter::contains_checked
    #[inline]
    pub fn contains<Q: AsRef<[u8]>>(&self, kmer: Q) -> bool {
        self.contains_encoded(self.encoder.encode_kmer(kmer))
    }

    /// Return whether a k-mer is present in this counter (i.e., has a count of
    /// at least 1). If the bases and k-mer length are not valid for the
    /// [`KmerEncoder`] associated with this [`KmerCounter`], then `None` is
    /// returned.
    #[inline]
    pub fn contains_checked<Q: AsRef<[u8]>>(&self, kmer: Q) -> Option<bool> {
        Some(self.contains_encoded(self.encoder.encode_kmer_checked(kmer)?))
    }

    /// Get the count of an already encoded k-mer. If the k-mer is not present
    /// in the counter, then `0` is returned. The encoded k-mer must have been
    /// generated using the [`KmerEncoder`] associated with this
    /// [`KmerCounter`].
    #[inline]
    pub fn get_encoded(&self, encoded_kmer: E::EncodedKmer) -> usize {
        self.map.get(&encoded_kmer).copied().unwrap_or_default()
    }

    /// Get the count of a k-mer. If the k-mer is not present in the counter,
    /// then `0` is returned. The bases and k-mer length are assumed to be valid
    /// for the [`KmerEncoder`] associated with this [`KmerCounter`]. Consider
    /// [`get_checked`] when it is not known whether the bases and k-mer length
    /// will be valid.
    ///
    /// [`get_checked`]: KmerCounter::get_checked
    #[inline]
    pub fn get<Q: AsRef<[u8]>>(&self, kmer: Q) -> usize {
        self.get_encoded(self.encoder.encode_kmer(kmer))
    }

    /// Get the count of a k-mer. If the k-mer is not present in the counter,
    /// then `0` is returned. If the bases and k-mer length are not valid for
    /// the [`KmerEncoder`] associated with this [`KmerCounter`], then `None` is
    /// returned.
    #[inline]
    pub fn get_checked<Q: AsRef<[u8]>>(&self, kmer: Q) -> Option<usize> {
        Some(self.get_encoded(self.encoder.encode_kmer_checked(kmer)?))
    }

    /// Iterate over the decoded k-mers and counts in the counter.
    #[inline]
    pub fn iter_decoded(&self) -> KmerCounterDecodedIter<'_, MAX_LEN, E> {
        KmerCounterDecodedIter {
            map_into_iter: self.map.iter(),
            encoder:       &self.encoder,
        }
    }

    /// Iterate over the encoded k-mers and counts in the counter.
    #[inline]
    pub fn iter_encoded(&self) -> hash_map::Iter<'_, E::EncodedKmer, usize> {
        self.map.iter()
    }

    /// Insert all k-mers into the counter with at most one mismatch compared to
    /// the provided, already encoded k-mer. The original k-mer is also tallied.
    /// The encoded k-mer must have been generated using the [`KmerEncoder`]
    /// associated with this [`KmerCounter`].
    #[inline]
    pub fn insert_encoded_kmer_one_mismatch(&mut self, encoded_kmer: E::EncodedKmer) {
        self.insert_encoded_kmer(encoded_kmer);
        for variant in self.encoder.get_variants_one_mismatch(encoded_kmer) {
            self.insert_encoded_kmer(variant);
        }
    }

    /// Insert all k-mers into the counter with at most one mismatch compared to
    /// the provided k-mer. The original k-mer is also tallied. The bases and
    /// k-mer length are assumed to be valid for the [`KmerEncoder`] associated
    /// with this [`KmerCounter`]. Consider [`insert_kmer_one_mismatch_checked`]
    /// when it is not known whether the bases and k-mer length will be valid.
    ///
    /// [`insert_kmer_one_mismatch_checked`]:
    ///     KmerCounter::insert_kmer_one_mismatch_checked
    #[inline]
    pub fn insert_kmer_one_mismatch<Q: AsRef<[u8]>>(&mut self, kmer: Q) {
        self.insert_encoded_kmer_one_mismatch(self.encoder.encode_kmer(kmer));
    }

    /// Insert all k-mers into the counter with at most one mismatch compared to
    /// the provided k-mer. The original k-mer is also tallied. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerCounter`], then `false` is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_one_mismatch_checked<Q: AsRef<[u8]>>(&mut self, kmer: Q) -> bool {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer_one_mismatch(encoded_kmer);
        true
    }

    /// Insert all k-mers from a sequence into the [`KmerCounter`]. The bases in
    /// the sequence must be valid for the [`KmerEncoder`] associated with this
    /// [`KmerCounter`].
    #[inline]
    pub fn insert_from_sequence<Q: AsRef<[u8]>>(&mut self, seq: Q) {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer(encoded_kmer);
        }
    }

    /// Insert all k-mers from a sequence into the [`KmerCounter`], in addition
    /// to all k-mers with up to one mismatch from those in the sequence. The
    /// bases in the sequence must be valid for the [`KmerEncoder`] associated
    /// with this [`KmerCounter`].
    #[inline]
    pub fn insert_from_sequence_one_mismatch<Q: AsRef<[u8]>>(&mut self, seq: Q) {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer_one_mismatch(encoded_kmer);
        }
    }

    /// Return the indices of the leftmost occurrence of any of the k-mers in
    /// this [`KmerCounter`] within a provided sequence. The bases in the
    /// sequence must
    /// be valid for the [`KmerEncoder`] associated with this [`KmerCounter`].
    /// If no occurrence is found, then `None` is returned.
    pub fn find_kmers<Q: AsRef<[u8]>>(&self, seq: Q) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder.iter_from_sequence(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                return Some(i..i + self.get_kmer_length());
            }
        }
        None
    }

    /// Return the indices of the rightmost occurrence of any of the k-mers in
    /// this [`KmerCounter`] within a provided sequence. The bases in the
    /// sequence must be valid for the [`KmerEncoder`] associated with this
    /// [`KmerCounter`]. If no occurrence is found, then `None` is returned.
    pub fn find_kmers_rev<Q: AsRef<[u8]>>(&self, seq: Q) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder.iter_from_sequence_rev(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                let end = seq.as_ref().len() - i;
                return Some(end - self.get_kmer_length()..end);
            }
        }
        None
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> IntoIterator for KmerCounter<MAX_LEN, E, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    type Item = (Kmer<MAX_LEN>, usize);
    type IntoIter = KmerCounterDecodedIntoIter<MAX_LEN, E, S>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            map_into_iter: self.map.into_iter(),
            encoder:       self.encoder,
        }
    }
}

/// An iterator over a [`KmerCounter`] yielding decoded k-mers and their counts.
/// The iterator consumes the original counter.
pub struct KmerCounterDecodedIntoIter<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen, {
    pub(crate) map_into_iter: <HashMap<E::EncodedKmer, usize, S> as IntoIterator>::IntoIter,
    pub(crate) encoder:       E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> Iterator for KmerCounterDecodedIntoIter<MAX_LEN, E, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    type Item = (Kmer<MAX_LEN>, usize);

    #[inline]
    fn next(&mut self) -> Option<(Kmer<MAX_LEN>, usize)> {
        self.map_into_iter.next().map(|(x, c)| (self.encoder.decode_kmer(x), c))
    }
}

/// An iterator over a [`KmerCounter`] yielding decoded k-mers and their
/// counts. The iterator takes the original counter by reference.
pub struct KmerCounterDecodedIter<'a, const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen, {
    pub(crate) map_into_iter: hash_map::Iter<'a, E::EncodedKmer, usize>,
    pub(crate) encoder:       &'a E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> Iterator for KmerCounterDecodedIter<'_, MAX_LEN, E>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    type Item = (Kmer<MAX_LEN>, usize);

    #[inline]
    fn next(&mut self) -> Option<(Kmer<MAX_LEN>, usize)> {
        self.map_into_iter.next().map(|(x, &c)| (self.encoder.decode_kmer(*x), c))
    }
}
