//! Defines a collection representing a set of encoded k-mers.

use crate::{
    kmer::{EncodedKmerCollection, FindKmersInSeq, Kmer, KmerEncode, KmerEncoder, KmerError, SupportedMismatchNumber},
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
/// 1. A single k-mer can be inserted with [`insert_kmer`]
/// 2. A k-mer as well as similar k-mers (up to `N` mismatches) can be inserted
///    with [`insert_kmer_with_variants`]
/// 3. Multiple k-mers from an iterator can be inserted with
///    [`insert_from_iter`]
/// 4. Overlapping k-mers from a sequence can be inserted with
///    [`insert_from_sequence`]
/// 5. Overlapping k-mers from a sequence with mismatches can be inserted with
///    [`insert_from_sequence_with_variants`]
///
/// After a k-mer set is populated, it can be used in multiple ways:
///
/// 1. Check for k-mer with [`contains`]
/// 2. Iteration: [`iter_encoded`] and [`iter_decoded`] provide the k-mers
///    without duplicates
/// 3. Set operations: encoded and decoded iterators for set operations are
///    implemented, including difference, intersection, symmetric difference,
///    and union
/// 4. Search for the k-mers within a sequence using [`FindKmersInSeq`] (or the
///    related trait [`FindKmers`]).
///
/// Consider using the alias [`ThreeBitKmerSet`], unless you are using a custom
/// [`KmerEncoder`].
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
/// [`ThreeBitKmerSet`]: crate::kmer::encoders::three_bit::ThreeBitKmerSet
/// [`SupportedKmerLen`]: super::SupportedKmerLen
/// [`insert_kmer`]: KmerSet::insert_kmer
/// [`FindKmers`]: super::FindKmers
/// [`contains`]: KmerSet::contains
/// [`iter_encoded`]: KmerSet::iter_encoded
/// [`iter_decoded`]: KmerSet::iter_decoded
/// [`insert_from_iter`]: KmerSet::insert_from_iter
/// [`insert_kmer_with_variants`]: KmerSet::insert_kmer_with_variants
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct KmerSet<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S = RandomState>
where
    S: BuildHasher, {
    set:     HashSet<E::EncodedKmer, S>,
    encoder: E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> KmerSet<MAX_LEN, E> {
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

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerSet<MAX_LEN, E, S> {
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

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerSet<MAX_LEN, E, S> {
    /// Inserts a k-mer into the [`KmerSet`].
    ///
    /// The k-mer can be either encoded or decoded (in which case it is encoded
    /// before insertion). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerSet`]. If it is decoded,
    /// it must be of length `self.kmer_length()`.
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
    /// it must be of length `self.kmer_length()`.
    #[inline]
    pub fn insert_kmer_with_variants<const N: usize>(&mut self, kmer: &impl KmerEncode<MAX_LEN, E>)
    where
        E::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        self.encoder
            .get_variants::<N>(kmer.encode_kmer(&self.encoder))
            .for_each(|variant| self.insert_kmer(&variant));
    }

    /// Inserts k-mers from an iterator into the [`KmerSet`].
    ///
    /// The k-mers can be either encoded or decoded (in which case it is encoded
    /// before insertion). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerSet`]. If it is decoded,
    /// it must be of length `self.kmer_length()`.
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
        E::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_kmer_with_variants::<N>(&encoded_kmer);
        }
    }

    /// Checks whether the [`KmerSet`] contains a k-mer.
    ///
    /// The k-mers can be either encoded or decoded (in which case it is encoded
    /// before checking). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerSet`]. If it is decoded,
    /// it must be of length `self.kmer_length()`.
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
    /// `other`. The two sets must have the same k-mer length.
    #[inline]
    pub fn difference_encoded<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> Copied<hash_set::Difference<'a, E::EncodedKmer, S>> {
        self.set.difference(&other.set).copied()
    }

    /// Returns an iterator over the decoded k-mers representing a set
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` but not in
    /// `other`. The two sets must have the same k-mer length.
    #[inline]
    pub fn difference_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.difference_encoded(other))
    }

    /// Returns an iterator over the encoded k-mers representing a set
    /// intersection.
    ///
    /// In other words, this returns the k-mers that are in both `self` and
    /// `other`. The two sets must have the same k-mer length.
    #[inline]
    pub fn intersection<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> Copied<hash_set::Intersection<'a, E::EncodedKmer, S>> {
        self.set.intersection(&other.set).copied()
    }

    /// Returns an iterator over the decoded k-mers representing a set
    /// intersection.
    ///
    /// In other words, this returns the k-mers that are in both `self` and
    /// `other`. The two sets must have the same k-mer length.
    #[inline]
    pub fn intersection_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.intersection(other))
    }

    /// Returns an iterator over the encoded k-mers representing a set symmetric
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` or in `other`
    /// but not in both. The two sets must have the same k-mer length.
    #[inline]
    pub fn symmetric_difference<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> Copied<hash_set::SymmetricDifference<'a, E::EncodedKmer, S>> {
        self.set.symmetric_difference(&other.set).copied()
    }

    /// Returns an iterator over the decoded k-mers representing a set symmetric
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` or in `other`
    /// but not in both. The two sets must have the same k-mer length.
    #[inline]
    pub fn symmetric_difference_decoded<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.symmetric_difference(other))
    }

    /// Returns an iterator over the encoded k-mers representing a set union.
    ///
    /// In other words, this returns the k-mers that are in `self` or `other`,
    /// without duplicates. The two sets must have the same k-mer length.
    #[inline]
    pub fn union<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> Copied<hash_set::Union<'a, E::EncodedKmer, S>> {
        self.set.union(&other.set).copied()
    }

    /// Returns an iterator over the decoded k-mers representing a set union.
    ///
    /// In other words, this returns the k-mers that are in `self` or `other`,
    /// without duplicates. The two sets must have the same k-mer length.
    #[inline]
    pub fn union_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.union(other))
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> EncodedKmerCollection<MAX_LEN>
    for KmerSet<MAX_LEN, E, S>
{
    type Encoder = E;
    type EncodedKmer = E::EncodedKmer;

    #[inline]
    fn encoder(&self) -> &Self::Encoder {
        &self.encoder
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> FindKmersInSeq<MAX_LEN> for KmerSet<MAX_LEN, E, S> {
    #[inline]
    fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, Self::Encoder>, {
        self.contains(kmer)
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> IntoIterator for KmerSet<MAX_LEN, E, S> {
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
pub struct KmerSetDecodedIntoIter<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> {
    pub(crate) set_into_iter: <HashSet<E::EncodedKmer, S> as IntoIterator>::IntoIter,
    pub(crate) encoder:       E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> Iterator for KmerSetDecodedIntoIter<MAX_LEN, E, S> {
    type Item = Kmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Kmer<MAX_LEN>> {
        self.set_into_iter.next().map(|x| self.encoder.decode_kmer(x))
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> Len for KmerSet<MAX_LEN, E, S>
where
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
