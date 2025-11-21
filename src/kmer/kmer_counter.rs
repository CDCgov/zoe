//! Defines a collection storing encoded k-mers and their counts.

use crate::{
    kmer::{
        EncodedKmerCollection, FindKmersInSeq, Kmer, KmerEncode, KmerError, SupportedMismatchNumber, encoders::KmerEncoder,
    },
    prelude::Len,
};
use std::{
    collections::{HashMap, hash_map},
    hash::{BuildHasher, RandomState},
    ops::Index,
};

/// A [`HashMap`] of k-mers and their counts, with k-mers stored in an encoded
/// format.
///
/// K-mers can be tallied into a [`KmerCounter`] multiple ways:
///
/// 1. A single k-mer can be tallied with [`tally_kmer`]
/// 2. A k-mer as well as similar k-mers (up to `N` mismatches) can be tallied
///    with [`tally_kmer_with_variants`]
/// 3. Multiple k-mers from an iterator can be tallied with [`tally_from_iter`]
/// 4. Overlapping k-mers from a sequence can be tallied with
///    [`tally_from_sequence`]
/// 5. Overlapping k-mers from a sequence with mismatches can be tallied with
///    [`tally_from_sequence_with_variants`]
///
/// After a k-mer counter is populated, it can be used in multiple ways:
///
/// 1. Check for k-mer with [`contains`], or get its count with [`get`]
/// 2. Iteration: [`iter_encoded`] and [`iter_decoded`] provide the k-mers and
///    counts, while [`keys_encoded`] and [`keys_decoded`] provide just the
///    k-mers
/// 3. Search for the k-mers within a sequence using [`FindKmersInSeq`] (or the
///    related trait [`FindKmers`])
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
/// [`ThreeBitKmerCounter`]:
///     crate::kmer::encoders::three_bit::ThreeBitKmerCounter
/// [`SupportedKmerLen`]: super::SupportedKmerLen
/// [`tally_kmer`]: KmerCounter::tally_kmer
/// [`tally_kmer_with_variants`]: KmerCounter::tally_kmer_with_variants
/// [`tally_from_iter`]: KmerCounter::tally_from_iter
/// [`tally_from_sequence`]: KmerCounter::tally_from_sequence
/// [`tally_from_sequence_with_variants`]:
///     KmerCounter::tally_from_sequence_with_variants
/// [`contains`]: KmerCounter::contains
/// [`get`]: KmerCounter::get
/// [`iter_encoded`]: KmerCounter::iter_encoded
/// [`iter_decoded`]: KmerCounter::iter_decoded
/// [`keys_encoded`]: KmerCounter::keys_encoded
/// [`keys_decoded`]: KmerCounter::keys_decoded
/// [`FindKmers`]: super::FindKmers
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct KmerCounter<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S = RandomState>
where
    S: BuildHasher, {
    /// The hashmap storing the encoded k-mers and their counts.
    map:     HashMap<E::EncodedKmer, usize, S>,
    /// The encoder used to encode the k-mers.
    encoder: E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> KmerCounter<MAX_LEN, E> {
    /// Creates a new [`KmerCounter`] with the specified k-mer length.
    ///
    /// ## Errors
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
    /// ## Errors
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

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> Index<E::EncodedKmer> for KmerCounter<MAX_LEN, E, S> {
    type Output = usize;

    #[inline]
    fn index(&self, index: E::EncodedKmer) -> &Self::Output {
        &self.map[&index]
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerCounter<MAX_LEN, E, S> {
    /// Tallies the k-mer in the [`KmerCounter`].
    ///
    /// If the k-mer is not present in the counter, this inserts it with a count
    /// of 1. Otherwise, this increments its count.
    ///
    /// The k-mer can be either encoded or decoded (in which case it is encoded
    /// before tallying). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerCounter`]. If it is
    /// decoded, it must be of length `self.kmer_length()`.
    pub fn tally_kmer<K>(&mut self, kmer: &K)
    where
        K: KmerEncode<MAX_LEN, E>, {
        *self.map.entry(kmer.encode_kmer(&self.encoder)).or_default() += 1;
    }

    /// Tallies all k-mers into the [`KmerCounter`] with at most `N` mismatches
    /// compared to the provided k-mer.
    ///
    /// The original k-mer is also tallied. The original k-mer can be either
    /// encoded or decoded (in which case it is encoded before tallying and
    /// variant generation). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerCounter`]. If it is decoded,
    /// it must be of length `self.kmer_length()`.
    #[inline]
    pub fn tally_kmer_with_variants<const N: usize>(&mut self, kmer: &impl KmerEncode<MAX_LEN, E>)
    where
        E::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        self.encoder
            .get_variants::<N>(kmer.encode_kmer(&self.encoder))
            .for_each(|variant| self.tally_kmer(&variant));
    }

    /// Tallies k-mers from an iterator in the [`KmerCounter`].
    ///
    /// The k-mers can be either encoded or decoded (in which case it is encoded
    /// before insertion). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerCounter`]. If it is
    /// decoded, it must be of length `self.kmer_length()`.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// When there is a choice, it is more efficient to use an iterator over
    /// encoded k-mers rather than decoded ones.
    ///
    /// </div>
    #[inline]
    pub fn tally_from_iter<I: IntoIterator<Item: KmerEncode<MAX_LEN, E>>>(&mut self, iter: I) {
        iter.into_iter().for_each(|kmer| self.tally_kmer(&kmer));
    }

    /// Tallies all overlapping k-mers from a sequence in the [`KmerCounter`].
    #[inline]
    pub fn tally_from_sequence(&mut self, seq: impl AsRef<[u8]>) {
        self.encoder.iter_from_sequence(&seq).for_each(|kmer| self.tally_kmer(&kmer));
    }

    /// Tallies all k-mers from a sequence in the [`KmerCounter`], in addition
    /// to all k-mers with up to N mismatches from those in the sequence.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::kmer::encoders::three_bit::ThreeBitKmerCounter;
    /// let mut counter = ThreeBitKmerCounter::<8>::new(8).unwrap();
    /// let seq = b"GATAGGGGATTGT";
    /// counter.tally_from_sequence_with_variants::<2>(seq);
    /// ```
    #[inline]
    pub fn tally_from_sequence_with_variants<const N: usize>(&mut self, seq: impl AsRef<[u8]>)
    where
        E::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.tally_kmer_with_variants::<N>(&encoded_kmer);
        }
    }

    /// Checks whether the [`KmerCounter`] contains a k-mer with a nonzero
    /// count.
    ///
    /// The k-mer can be either encoded or decoded (in which case it is encoded
    /// before checking). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerCounter`]. If it is
    /// decoded, it must be of length `self.kmer_length()`.
    #[inline]
    #[must_use]
    pub fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.map.contains_key(&kmer.encode_kmer(&self.encoder))
    }

    /// Gets the count of a k-mer.
    ///
    /// If the k-mer is not present in the counter, then `0` is returned. The
    /// k-mer can be either encoded or decoded (in which case it is encoded
    /// before checking). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerCounter`]. If it is
    /// decoded, it must be of length `self.kmer_length()`.
    #[inline]
    pub fn get<K>(&self, kmer: &K) -> usize
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.map.get(&kmer.encode_kmer(&self.encoder)).copied().unwrap_or_default()
    }

    /// Returns an iterator over the encoded k-mers and their counts.
    #[inline]
    pub fn iter_encoded(&self) -> hash_map::Iter<'_, E::EncodedKmer, usize> {
        self.map.iter()
    }

    /// Returns an iterator over the decoded k-mers and their counts.
    #[inline]
    pub fn iter_decoded(&self) -> impl Iterator<Item = (Kmer<MAX_LEN>, &usize)> {
        self.map.iter().map(|(k, c)| (self.encoder.decode_kmer(*k), c))
    }

    /// Returns an iterator over the decoded k-mers in the counter without
    /// duplicates.
    #[inline]
    pub fn keys_decoded(&self) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.map.keys().map(|encoded_kmer| self.encoder().decode_kmer(*encoded_kmer))
    }

    /// Returns an iterator over the encoded k-mers in the counter without
    /// duplicates.
    #[inline]
    pub fn keys_encoded(&self) -> impl Iterator<Item = &E::EncodedKmer> {
        self.map.keys()
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

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> FindKmersInSeq<MAX_LEN> for KmerCounter<MAX_LEN, E, S> {
    #[inline]
    fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, Self::Encoder>, {
        self.contains(kmer)
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> IntoIterator for KmerCounter<MAX_LEN, E, S> {
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
pub struct KmerCounterDecodedIntoIter<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> {
    pub(crate) map_into_iter: <HashMap<E::EncodedKmer, usize, S> as IntoIterator>::IntoIter,
    pub(crate) encoder:       E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> Iterator for KmerCounterDecodedIntoIter<MAX_LEN, E, S> {
    type Item = (Kmer<MAX_LEN>, usize);

    #[inline]
    fn next(&mut self) -> Option<(Kmer<MAX_LEN>, usize)> {
        self.map_into_iter.next().map(|(x, c)| (self.encoder.decode_kmer(x), c))
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> Len for KmerCounter<MAX_LEN, E, S>
where
    S: BuildHasher,
{
    #[inline]
    fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.map.len()
    }
}
