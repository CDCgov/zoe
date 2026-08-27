//! Defines [`KmerSet`], which holds a set of encoded k-mers using a hash set.

use crate::{
    kmer::{EncodedKmerCollection, FindKmersInSeq, GetVariants, Kmer, KmerEncode, KmerEncoder, KmerError},
    prelude::Len,
};
use std::{
    collections::{HashSet, hash_set},
    hash::{BuildHasher, RandomState},
    iter::Copied,
};

/// A [`HashSet`] of k-mers, stored in an encoded format.
///
/// K-mers can be inserted into a [`KmerSet`] multiple ways:
///
/// - A single k-mer can be inserted with [`insert_kmer`]
/// - A k-mer as well as similar k-mers (up to `N` mismatches) can be inserted
///   with [`insert_kmer_with_variants`]
/// - Multiple k-mers from an iterator can be inserted with [`insert_from_iter`]
/// - Overlapping k-mers from a sequence can be inserted with
///   [`insert_from_sequence`]
/// - Overlapping k-mers from a sequence with mismatches can be inserted with
///   [`insert_from_sequence_with_variants`]
///
/// After a k-mer set is populated, it can be used in multiple ways:
///
/// - Check for k-mer with [`contains`]
/// - Iteration: [`iter_encoded`] and [`iter_decoded`] provide the k-mers in the
///   set without duplicates
/// - Set operations: encoded and decoded iterators for set operations between
///   two k-mer sets are implemented, including difference, intersection,
///   symmetric difference, and union
/// - Search for the k-mers within a sequence using [`FindKmersInSeq`] (or the
///   related trait [`FindKmers`])
///
/// <div class="warning tip">
///
/// **Tip**
///
/// For guidance on picking the appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// </div>
///
/// [`find_in_seq`]: KmerSet::find_in_seq
/// [`find_in_seq_rev`]: KmerSet::find_in_seq_rev
/// [`insert_from_sequence`]: KmerSet::insert_from_sequence
/// [`insert_from_sequence_with_variants`]:
///     KmerSet::insert_from_sequence_with_variants
/// [`SupportedKmerLen`]: crate::kmer::SupportedKmerLen
/// [`insert_kmer`]: KmerSet::insert_kmer
/// [`FindKmers`]: crate::kmer::FindKmers
/// [`contains`]: KmerSet::contains
/// [`iter_encoded`]: KmerSet::iter_encoded
/// [`iter_decoded`]: KmerSet::iter_decoded
/// [`insert_from_iter`]: KmerSet::insert_from_iter
/// [`insert_kmer_with_variants`]: KmerSet::insert_kmer_with_variants
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct KmerSet<const MAX_LEN: usize, E, S = RandomState>
where
    E: KmerEncoder<MAX_LEN>,
    S: BuildHasher, {
    /// The hashset storing the encoded k-mers.
    set:     HashSet<E::EncodedKmer, S>,
    /// The encoder used to encode the k-mers.
    encoder: E,
}

impl<const MAX_LEN: usize, E> KmerSet<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN>,
{
    /// Creates a new [`KmerSet`] with the specified k-mer length.
    ///
    /// ## Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn new(kmer_length: usize) -> Result<Self, KmerError> {
        Ok(Self {
            set:     HashSet::default(),
            encoder: E::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, E, S> KmerSet<MAX_LEN, E, S>
where
    E: KmerEncoder<MAX_LEN>,
    S: BuildHasher,
{
    /// Creates a new [`KmerSet`] with the specified k-mer length and hasher.
    ///
    /// ## Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn with_hasher(kmer_length: usize, hasher: S) -> Result<Self, KmerError> {
        Ok(Self {
            set:     HashSet::with_hasher(hasher),
            encoder: E::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, E, S> KmerSet<MAX_LEN, E, S>
where
    E: KmerEncoder<MAX_LEN>,
    S: BuildHasher,
{
    /// Inserts a k-mer into the [`KmerSet`].
    ///
    /// The k-mer can be either encoded or decoded (in which case it is encoded
    /// before insertion). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerSet`]. If it is decoded,
    /// it must be of length [`Self::kmer_length`].
    #[inline]
    pub fn insert_kmer<K>(&mut self, kmer: &K)
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.set.insert(kmer.encode_kmer(&self.encoder));
    }

    /// Inserts all k-mers into the [`KmerSet`] with at most `N` mismatches
    /// compared to the provided k-mer.
    ///
    /// The original k-mer is also inserted. The original k-mer can be either
    /// encoded or decoded (in which case it is encoded before insertion and
    /// variant generation). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerSet`]. If it is decoded,
    /// it must be of length [`Self::kmer_length`].
    #[inline]
    pub fn insert_kmer_with_variants<const N: usize>(&mut self, kmer: &impl KmerEncode<MAX_LEN, E>)
    where
        E: GetVariants<N, MAX_LEN>, {
        self.encoder
            .get_variants::<N>(kmer.encode_kmer(&self.encoder))
            .for_each(|variant| self.insert_kmer(&variant));
    }

    /// Inserts k-mers from an iterator into the [`KmerSet`].
    ///
    /// The k-mers can be either encoded or decoded (in which case it is encoded
    /// before insertion). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerSet`]. If it is decoded,
    /// it must be of length [`Self::kmer_length`].
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
    pub fn insert_from_iter<I: IntoIterator<Item: KmerEncode<MAX_LEN, E>>>(&mut self, iter: I) {
        iter.into_iter().for_each(|kmer| self.insert_kmer(&kmer));
    }

    /// Inserts all overlapping k-mers from a sequence into the [`KmerSet`].
    #[inline]
    pub fn insert_from_sequence(&mut self, seq: impl AsRef<[u8]>) {
        self.encoder.iter_from_sequence(&seq).for_each(|kmer| self.insert_kmer(&kmer));
    }

    /// Insert all k-mers from a sequence into the [`KmerSet`], in addition to
    /// all k-mers with up to `N` mismatches from those in the sequence.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::kmer::encoders::three_bit::ThreeBitKmerSet;
    /// let mut set = ThreeBitKmerSet::<8>::new(8).unwrap();
    /// let seq = b"GATAGGGGATTGT";
    /// set.insert_from_sequence_with_variants::<2>(seq);
    /// ```
    #[inline]
    pub fn insert_from_sequence_with_variants<const N: usize>(&mut self, seq: impl AsRef<[u8]>)
    where
        E: GetVariants<N, MAX_LEN>, {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_kmer_with_variants::<N>(&encoded_kmer);
        }
    }

    /// Checks whether the [`KmerSet`] contains a k-mer.
    ///
    /// The k-mers can be either encoded or decoded (in which case it is encoded
    /// before checking). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerSet`]. If it is decoded,
    /// it must be of length [`Self::kmer_length`].
    #[inline]
    #[must_use]
    pub fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.set.contains(&kmer.encode_kmer(&self.encoder))
    }

    /// Returns an iterator over the encoded k-mers in the set.
    #[inline]
    pub fn iter_encoded(&self) -> Copied<hash_set::Iter<'_, E::EncodedKmer>> {
        self.set.iter().copied()
    }

    /// Returns an iterator over the decoded k-mers in the set.
    #[inline]
    pub fn iter_decoded(&self) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.iter_encoded())
    }

    /// Returns an iterator over the encoded k-mers representing a set
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` but not in
    /// `other`.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`KmerSet`]s are not
    /// equal.
    #[inline]
    pub fn difference_encoded<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> Copied<hash_set::Difference<'a, E::EncodedKmer, S>> {
        assert_eq!(self.kmer_length(), other.kmer_length());

        self.set.difference(&other.set).copied()
    }

    /// Returns an iterator over the decoded k-mers representing a set
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` but not in
    /// `other`.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`KmerSet`]s are not equal.
    #[inline]
    pub fn difference_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.difference_encoded(other))
    }

    /// Returns an iterator over the encoded k-mers representing a set
    /// intersection.
    ///
    /// In other words, this returns the k-mers that are in both `self` and
    /// `other`.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`KmerSet`]s are not equal.
    #[inline]
    pub fn intersection_encoded<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> Copied<hash_set::Intersection<'a, E::EncodedKmer, S>> {
        assert_eq!(self.kmer_length(), other.kmer_length());

        self.set.intersection(&other.set).copied()
    }

    /// Returns an iterator over the decoded k-mers representing a set
    /// intersection.
    ///
    /// In other words, this returns the k-mers that are in both `self` and
    /// `other`.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`KmerSet`]s are not equal.
    #[inline]
    pub fn intersection_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.intersection_encoded(other))
    }

    /// Returns an iterator over the encoded k-mers representing a set symmetric
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` or in `other`
    /// but not in both.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`KmerSet`]s are not equal.
    #[inline]
    pub fn symmetric_difference_encoded<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> Copied<hash_set::SymmetricDifference<'a, E::EncodedKmer, S>> {
        assert_eq!(self.kmer_length(), other.kmer_length());

        self.set.symmetric_difference(&other.set).copied()
    }

    /// Returns an iterator over the decoded k-mers representing a set symmetric
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` or in `other`
    /// but not in both.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`KmerSet`]s are not equal.
    #[inline]
    pub fn symmetric_difference_decoded<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.symmetric_difference_encoded(other))
    }

    /// Returns an iterator over the encoded k-mers representing a set union.
    ///
    /// In other words, this returns the k-mers that are in `self` or `other`,
    /// without duplicates.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`KmerSet`]s are not equal.
    #[inline]
    pub fn union_encoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> Copied<hash_set::Union<'a, E::EncodedKmer, S>> {
        assert_eq!(self.kmer_length(), other.kmer_length());

        self.set.union(&other.set).copied()
    }

    /// Returns an iterator over the decoded k-mers representing a set union.
    ///
    /// In other words, this returns the k-mers that are in `self` or `other`,
    /// without duplicates.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`KmerSet`]s are not equal.
    #[inline]
    pub fn union_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.union_encoded(other))
    }
}

impl<const MAX_LEN: usize, E, S> EncodedKmerCollection<MAX_LEN> for KmerSet<MAX_LEN, E, S>
where
    E: KmerEncoder<MAX_LEN>,
    S: BuildHasher,
{
    type Encoder = E;
    type EncodedKmer = E::EncodedKmer;

    #[inline]
    fn encoder(&self) -> &Self::Encoder {
        &self.encoder
    }
}

impl<const MAX_LEN: usize, E, S> FindKmersInSeq<MAX_LEN> for KmerSet<MAX_LEN, E, S>
where
    E: KmerEncoder<MAX_LEN>,
    S: BuildHasher,
{
    #[inline]
    fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, Self::Encoder>, {
        self.contains(kmer)
    }
}

impl<const MAX_LEN: usize, E, S> IntoIterator for KmerSet<MAX_LEN, E, S>
where
    E: KmerEncoder<MAX_LEN>,
    S: BuildHasher,
{
    type Item = Kmer<MAX_LEN>;
    type IntoIter = KmerSetDecodedIntoIter<MAX_LEN, E, S>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            set_into_iter: self.set.into_iter(),
            encoder:       self.encoder,
        }
    }
}

/// An iterator over a [`KmerSet`] yielding decoded k-mers. The iterator
/// consumes the original set.
pub struct KmerSetDecodedIntoIter<const MAX_LEN: usize, E, S>
where
    E: KmerEncoder<MAX_LEN>, {
    set_into_iter: <HashSet<E::EncodedKmer, S> as IntoIterator>::IntoIter,
    encoder:       E,
}

impl<const MAX_LEN: usize, E, S> Iterator for KmerSetDecodedIntoIter<MAX_LEN, E, S>
where
    E: KmerEncoder<MAX_LEN>,
{
    type Item = Kmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Kmer<MAX_LEN>> {
        self.set_into_iter.next().map(|x| self.encoder.decode_kmer(x))
    }
}

impl<const MAX_LEN: usize, E, S> Len for KmerSet<MAX_LEN, E, S>
where
    E: KmerEncoder<MAX_LEN>,
    S: BuildHasher,
{
    #[inline]
    fn is_empty(&self) -> bool {
        self.set.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.set.len()
    }
}
