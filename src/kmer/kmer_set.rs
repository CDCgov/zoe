use super::{EncodedKmerCollection, KmerCollectionContains, KmerEncode, SupportedMismatchNumber};
use crate::{
    kmer::{Kmer, KmerEncoder, KmerError},
    prelude::Len,
};
use std::{
    collections::HashSet,
    hash::{BuildHasher, RandomState},
};

/// A [`KmerSet`] holds a set of encoded k-mers and provides methods for
/// efficiently using them. For instance, a [`KmerSet`] can be used to find the
/// leftmost or rightmost occurrence of k-mers in a sequence, using
/// [`find_in_seq`] or [`find_in_seq_rev`]. A [`KmerSet`] can also be
/// constructed using [`insert_from_sequence`] or
/// [`insert_from_sequence_with_variants`].
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
/// [`ThreeBitKmerSet`]: super::ThreeBitKmerSet
/// [`SupportedKmerLen`]: super::SupportedKmerLen
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
    /// Insert an already encoded k-mer into the set. The encoded k-mer must
    /// have been generated using the [`KmerEncoder`] associated with this
    /// [`KmerSet`].
    #[inline]
    pub fn insert_encoded_kmer(&mut self, encoded_kmer: E::EncodedKmer) {
        self.set.insert(encoded_kmer);
    }

    /// Insert a k-mer into the set. The bases and k-mer length are assumed to be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. Consider
    /// [`insert_kmer_checked`] when it is not known whether the bases and k-mer
    /// length will be valid.
    ///
    /// [`insert_kmer_checked`]: KmerSet::insert_kmer_checked
    #[inline]
    pub fn insert_kmer(&mut self, kmer: impl AsRef<[u8]>) {
        self.insert_encoded_kmer(self.encoder.encode_kmer(kmer));
    }

    /// Insert a k-mer into the set. If the bases and k-mer length are not valid
    /// for the [`KmerEncoder`] associated with this [`KmerSet`], then `false`
    /// is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_checked(&mut self, kmer: impl AsRef<[u8]>) -> bool {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer(encoded_kmer);
        true
    }

    /// Insert k-mers from an iterator into the [`KmerSet`]. The k-mers can be
    /// either encoded or decoded. Encoded k-mers must have been generated using
    /// the [`KmerEncoder`] associated with this [`KmerSet`]. Decoded k-mers
    /// must have bases and k-mer lengths which are valid for the
    /// [`KmerEncoder`].
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
        for kmer in iter {
            self.insert_encoded_kmer(kmer.encode_kmer(&self.encoder));
        }
    }

    /// Visits the decoded k-mers in an iterator of encoded k-mers
    #[inline]
    fn decode_iter<'a>(&self, iter: impl Iterator<Item = &'a E::EncodedKmer>) -> impl Iterator<Item = Kmer<MAX_LEN>>
    where
        <E as KmerEncoder<MAX_LEN>>::EncodedKmer: 'a, {
        iter.map(|encoded_kmer| self.encoder.decode_kmer(*encoded_kmer))
    }

    /// Visits the encoded k-mers in the set.
    #[inline]
    pub fn iter_encoded(&self) -> impl Iterator<Item = &E::EncodedKmer> {
        self.set.iter()
    }

    /// Visits the decoded k-mers in the set.
    #[inline]
    pub fn iter_decoded(&self) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.decode_iter(self.iter_encoded())
    }

    /// Visits the encoded k-mers representing the difference, i.e., the k-mers
    /// that are in `self` but not in `other`.
    #[inline]
    pub fn difference_encoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = &'a E::EncodedKmer> {
        self.set.difference(&other.set)
    }

    /// Visits the decoded k-mers representing the difference, i.e., the k-mers
    /// that are in `self` but not in `other`.
    #[inline]
    pub fn difference_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.decode_iter(self.set.difference(&other.set))
    }

    /// Visits the encoded k-mers representing the intersection, i.e., the
    /// k-mers that are both in `self` and `other`.
    #[inline]
    pub fn intersection<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = &'a E::EncodedKmer> {
        self.set.intersection(&other.set)
    }

    /// Visits the decoded k-mers representing the intersection, i.e., the
    /// k-mers that are both in `self` and `other`.
    #[inline]
    pub fn intersection_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.decode_iter(self.intersection(other))
    }

    /// Visits the encoded k-mers representing the symmetric difference, i.e.,
    /// the k-mers that are in `self` or in `other` but not in both.
    #[inline]
    pub fn symmetric_difference<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> impl Iterator<Item = &'a E::EncodedKmer> {
        self.set.symmetric_difference(&other.set)
    }

    /// Visits the decoded k-mers representing the symmetric difference, i.e.,
    /// the k-mers that are in `self` or in `other` but not in both.
    #[inline]
    pub fn symmetric_difference_decoded<'a>(
        &'a self, other: &'a KmerSet<MAX_LEN, E, S>,
    ) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.decode_iter(self.symmetric_difference(other))
    }

    /// Visits the encoded k-mers representing the union, i.e., the k-mers that
    /// are in `self` or `other`, without duplicates.
    #[inline]
    pub fn union<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = &'a E::EncodedKmer> {
        self.set.union(&other.set)
    }

    /// Visits the decoded k-mers representing the union, i.e., the k-mers that
    /// are in `self` or `other`, without duplicates.
    #[inline]
    pub fn union_decoded<'a>(&'a self, other: &'a KmerSet<MAX_LEN, E, S>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.decode_iter(self.union(other))
    }

    /// Insert all k-mers into the set with at most `N` mismatches compared to the
    /// provided, already encoded k-mer. The original k-mer is also inserted.
    /// The encoded k-mer must have been generated using the encoder associated
    /// with this [`KmerSet`].
    ///
    /// `N` must be a supported number of mismatches. See
    /// [`SupportedMismatchNumber`] for more details.
    #[inline]
    pub fn insert_encoded_kmer_with_variants<const N: usize>(&mut self, encoded_kmer: E::EncodedKmer)
    where
        E::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        for variant in self.encoder.get_variants::<N>(encoded_kmer) {
            self.insert_encoded_kmer(variant);
        }
    }

    /// Insert all k-mers into the set with at most N mismatches compared to the
    /// provided k-mer. The original k-mer is also inserted. The bases and k-mer
    /// length are assumed to be valid for the [`KmerEncoder`] associated with
    /// this [`KmerSet`]. Consider [`insert_kmer_with_variants_checked`] when it
    /// is not known whether the bases and k-mer length will be valid.
    ///
    /// [`insert_kmer_with_variants_checked`]:
    ///     KmerSet::insert_kmer_with_variants_checked
    #[inline]
    pub fn insert_kmer_with_variants<const N: usize>(&mut self, kmer: impl AsRef<[u8]>)
    where
        E::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        self.insert_encoded_kmer_with_variants::<N>(self.encoder.encode_kmer(kmer));
    }

    /// Insert all k-mers into the set with at most N mismatches compared to the
    /// provided k-mer. The original k-mer is also inserted. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`], then `false` is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_with_variants_checked<const N: usize>(&mut self, kmer: impl AsRef<[u8]>) -> bool
    where
        E::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer_with_variants::<N>(encoded_kmer);
        true
    }

    /// Insert all k-mers from a sequence into the [`KmerSet`]. The bases in the
    /// sequence must be valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`].
    #[inline]
    pub fn insert_from_sequence(&mut self, seq: impl AsRef<[u8]>) {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer(encoded_kmer);
        }
    }

    /// Insert all k-mers from a sequence into the [`KmerSet`], in addition to
    /// all k-mers with up to N mismatches from those in the sequence. The bases
    /// in the sequence must be valid for the [`KmerEncoder`] associated with
    /// this [`KmerSet`].
    ///
    /// ```
    /// # use zoe::kmer::ThreeBitKmerSet;
    /// let mut set = ThreeBitKmerSet::<8>::new(8).unwrap();
    /// let seq = b"GATAGGGGATTGT";
    /// set.insert_from_sequence_with_variants::<2>(seq);
    /// ```
    #[inline]
    pub fn insert_from_sequence_with_variants<const N: usize>(&mut self, seq: impl AsRef<[u8]>)
    where
        E::MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer_with_variants::<N>(encoded_kmer);
        }
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

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerCollectionContains<MAX_LEN>
    for KmerSet<MAX_LEN, E, S>
{
    #[inline]
    fn contains_encoded(&self, kmer: Self::EncodedKmer) -> bool {
        self.set.contains(&kmer)
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
