use super::{
    EncodedKmerCollection, Kmer, KmerCollectionContains, KmerError, KmerLen, MismatchNumber, SupportedKmerLen,
    SupportedMismatchNumber,
};
use crate::{kmer::encoder::KmerEncoder, prelude::Len};
use std::{
    collections::{HashMap, hash_map},
    hash::{BuildHasher, RandomState},
    ops::Index,
};

/// A [`KmerCounter`] stores counts of encoded k-mers, or it can be considered
/// as a multiset. [`KmerCounter`] has many of the same methods as [`KmerSet`].
///
/// Consider using the alias [`ThreeBitKmerCounter`], unless you are using a
/// custom [`KmerEncoder`].
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
///
/// [`KmerSet`]: super::KmerSet
/// [`ThreeBitKmerCounter`]: super::ThreeBitKmerCounter
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct KmerCounter<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S = RandomState>
where
    S: BuildHasher, {
    map:     HashMap<E::EncodedKmer, usize, S>,
    encoder: E,
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

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerCounter<MAX_LEN, E, S> {
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
    pub fn insert_kmer(&mut self, kmer: impl AsRef<[u8]>) {
        self.insert_encoded_kmer(self.encoder.encode_kmer(kmer));
    }

    /// If the k-mer is present in this counter, then increment its count.
    /// Otherwise, add it to the counter with a count of 1. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerCounter`], then `false` is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_checked(&mut self, kmer: impl AsRef<[u8]>) -> bool {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer(encoded_kmer);
        true
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
    pub fn get(&self, kmer: impl AsRef<[u8]>) -> usize {
        self.get_encoded(self.encoder.encode_kmer(kmer))
    }

    /// Get the count of a k-mer. If the k-mer is not present in the counter,
    /// then `0` is returned. If the bases and k-mer length are not valid for
    /// the [`KmerEncoder`] associated with this [`KmerCounter`], then `None` is
    /// returned.
    #[inline]
    pub fn get_checked(&self, kmer: impl AsRef<[u8]>) -> Option<usize> {
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

    /// Insert all k-mers into the counter with at most `N` mismatches compared
    /// to the provided, already encoded k-mer. The original k-mer is also
    /// tallied. The encoded k-mer must have been generated using the
    /// [`KmerEncoder`] associated with this [`KmerCounter`].
    ///
    /// `N` must be a supported number of mismatches. See
    /// [`SupportedMismatchNumber`] for more details.
    #[inline]
    pub fn insert_encoded_kmer_with_variants<const N: usize>(&mut self, encoded_kmer: E::EncodedKmer)
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        for variant in self.encoder.get_variants::<N>(encoded_kmer) {
            self.insert_encoded_kmer(variant);
        }
    }

    /// Insert all k-mers into the counter with at most N mismatches compared to
    /// the provided k-mer. The original k-mer is also tallied. The bases and
    /// k-mer length are assumed to be valid for the [`KmerEncoder`] associated
    /// with this [`KmerCounter`]. Consider
    /// [`insert_kmer_with_variants_checked`] when it is not known whether the
    /// bases and k-mer length will be valid.
    ///
    /// [`insert_kmer_with_variants_checked`]:
    ///     KmerCounter::insert_kmer_with_variants_checked
    #[inline]
    pub fn insert_kmer_with_variants<const N: usize>(&mut self, kmer: impl AsRef<[u8]>)
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        self.insert_encoded_kmer_with_variants::<N>(self.encoder.encode_kmer(kmer));
    }

    /// Insert all k-mers into the counter with at most N mismatches compared to
    /// the provided k-mer. The original k-mer is also tallied. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerCounter`], then `false` is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_with_variants_checked<const N: usize>(&mut self, kmer: impl AsRef<[u8]>) -> bool
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer_with_variants::<N>(encoded_kmer);
        true
    }

    /// Insert all k-mers from a sequence into the [`KmerCounter`]. The bases in
    /// the sequence must be valid for the [`KmerEncoder`] associated with this
    /// [`KmerCounter`].
    #[inline]
    pub fn insert_from_sequence(&mut self, seq: impl AsRef<[u8]>) {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer(encoded_kmer);
        }
    }

    /// Insert all k-mers from a sequence into the [`KmerCounter`], in addition
    /// to all k-mers with up to N mismatches from those in the sequence. The
    /// bases in the sequence must be valid for the [`KmerEncoder`] associated
    /// with this [`KmerCounter`].
    ///
    /// ```
    /// # use zoe::kmer::ThreeBitKmerCounter;
    /// let mut counter = ThreeBitKmerCounter::<8>::new(8).unwrap();
    /// let seq = b"GATAGGGGATTGT";
    /// counter.insert_from_sequence_with_variants::<2>(seq);
    /// ```
    #[inline]
    pub fn insert_from_sequence_with_variants<const N: usize>(&mut self, seq: impl AsRef<[u8]>)
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer_with_variants::<N>(encoded_kmer);
        }
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> EncodedKmerCollection<MAX_LEN>
    for KmerCounter<MAX_LEN, E, S>
{
    type Encoder = E;
    type EncodedKmer = E::EncodedKmer;

    #[inline]
    fn encoder(&self) -> &Self::Encoder {
        &self.encoder
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerCollectionContains<MAX_LEN>
    for KmerCounter<MAX_LEN, E, S>
{
    #[inline]
    fn contains_encoded(&self, kmer: Self::EncodedKmer) -> bool {
        self.map.contains_key(&kmer)
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> IntoIterator for KmerCounter<MAX_LEN, E, S>
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

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> Len for KmerCounter<MAX_LEN, E> {
    #[inline]
    fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.map.len()
    }
}
