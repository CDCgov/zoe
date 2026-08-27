//! Defines [`IndexedKmerSet`], which holds a set of encoded k-mers using
//! indexing into a `Vec`.

use crate::{
    data::views::Len,
    kmer::{
        EncodedKmerCollection, FindKmersInSeq, GetVariants, Kmer, KmerEncode, KmerEncoder, KmerError, KmerIndex, KmerLen,
        SupportedKmerLen,
    },
};
use std::iter::Enumerate;

/// A set of k-mers, stored as a boolean vector indexed by the encoded k-mers,
/// rather than using hashing like [`KmerSet`].
///
/// For example, there are 64 possible 2-bit 3-mers, which can be encoded as
/// `u8`s in `[0, 63]`. If a 3-mer of `TCG` is encountered, it can be encoded as
/// `54`, which can then be used to index the vector and set index 54 to be
/// `true`.
///
/// K-mers can be inserted into an [`IndexedKmerSet`] multiple ways:
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
/// After an indexed k-mer set is populated, it can be used in multiple ways:
///
/// - Check for a k-mer with [`contains`]
/// - Iteration: [`iter_encoded`] and [`iter_decoded`] provide the k-mers in the
///   set without duplicates
/// - Set operations: encoded and decoded iterators for set operations between
///   two indexed k-mer sets are implemented, including difference,
///   intersection, symmetric difference, and union
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
/// [`insert_kmer`]: IndexedKmerSet::insert_kmer
/// [`HashSet`]: std::collections::HashSet
/// [`KmerSet`]: crate::kmer::collections::KmerSet
/// [`insert_kmer_with_variants`]: IndexedKmerSet::insert_kmer_with_variants
/// [`insert_from_iter`]: IndexedKmerSet::insert_from_iter
/// [`insert_from_sequence`]: IndexedKmerSet::insert_from_sequence
/// [`insert_from_sequence_with_variants`]:
///     IndexedKmerSet::insert_from_sequence_with_variants
/// [`contains`]: IndexedKmerSet::contains
/// [`iter_encoded`]: IndexedKmerSet::iter_encoded
/// [`iter_decoded`]: IndexedKmerSet::iter_decoded
/// [`FindKmers`]: crate::kmer::FindKmers
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct IndexedKmerSet<const MAX_LEN: usize, E>
where
    E: KmerEncoder<MAX_LEN>, {
    /// The vec indexed by the encoded kmers.
    vec:     Box<[bool]>,
    /// The encoder used to encode the k-mers.
    encoder: E,
}

impl<const MAX_LEN: usize, E> IndexedKmerSet<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    /// Creates a new [`IndexedKmerSet`] with the specified k-mer length.
    ///
    /// ## Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn new(kmer_length: usize) -> Result<Self, KmerError> {
        Ok(Self {
            vec:     vec![false; E::EncodedKmer::max_index_for_length(kmer_length)].into_boxed_slice(),
            encoder: E::new(kmer_length)?,
        })
    }

    /// Inserts a k-mer into the [`IndexedKmerSet`].
    ///
    /// The k-mer can be either encoded or decoded (in which case it is encoded
    /// before insertion). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`IndexedKmerSet`]. If it is
    /// decoded, it must be of length [`Self::kmer_length`].
    #[inline]
    pub fn insert_kmer<K>(&mut self, kmer: &K)
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.vec[kmer.encode_kmer(&self.encoder).as_usize()] = true;
    }

    /// Inserts all k-mers into the [`IndexedKmerSet`] with at most `N`
    /// mismatches compared to the provided k-mer.
    ///
    /// The original k-mer is also inserted. The original k-mer can be either
    /// encoded or decoded (in which case it is encoded before insertion and
    /// variant generation). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`IndexedKmerSet`]. If it is
    /// decoded, it must be of length [`Self::kmer_length`].
    #[inline]
    pub fn insert_kmer_with_variants<const N: usize>(&mut self, kmer: &impl KmerEncode<MAX_LEN, E>)
    where
        E: GetVariants<N, MAX_LEN>, {
        self.encoder
            .get_variants::<N>(kmer.encode_kmer(&self.encoder))
            .for_each(|variant| self.insert_kmer(&variant));
    }

    /// Inserts k-mers from an iterator into the [`IndexedKmerSet`].
    ///
    /// The k-mers can be either encoded or decoded (in which case it is encoded
    /// before insertion). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`IndexedKmerSet`]. If it is
    /// decoded, it must be of length [`Self::kmer_length`].
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

    /// Inserts all overlapping k-mers from a sequence into the
    /// [`IndexedKmerSet`].
    #[inline]
    pub fn insert_from_sequence(&mut self, seq: impl AsRef<[u8]>) {
        self.encoder.iter_from_sequence(&seq).for_each(|kmer| self.insert_kmer(&kmer));
    }

    /// Insert all k-mers from a sequence into the [`IndexedKmerSet`], in
    /// addition to all k-mers with up to `N` mismatches from those in the
    /// sequence.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::kmer::{IndexedKmerSet, encoders::three_bit::ThreeBitKmerEncoder};
    /// let mut set = IndexedKmerSet::<4, ThreeBitKmerEncoder<4>>::new(4).unwrap();
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

    /// Checks whether the [`IndexedKmerSet`] contains a k-mer.
    ///
    /// The k-mers can be either encoded or decoded (in which case it is encoded
    /// before checking). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this [`IndexedKmerSet`]. If it is
    /// decoded, it must be of length [`Self::kmer_length`].
    #[inline]
    #[must_use]
    pub fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.vec[kmer.encode_kmer(&self.encoder).as_usize()]
    }

    /// Returns an iterator over the encoded k-mers in the indexed set.
    #[inline]
    pub fn iter_encoded(&self) -> impl Iterator<Item = E::EncodedKmer> {
        self.vec
            .iter()
            .enumerate()
            .filter(|(_, seen)| **seen)
            .map(|(index, _)| E::EncodedKmer::from_usize(index))
    }

    /// Returns an iterator over the decoded k-mers in the indexed set.
    #[inline]
    pub fn iter_decoded(&self) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.iter_encoded().map(|encoded_kmer| self.encoder.decode_kmer(encoded_kmer))
    }

    /// Returns an iterator over the encoded k-mers representing a set
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` but not in
    /// `other`.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`IndexedKmerSet`]s are not
    /// equal.
    #[inline]
    pub fn difference_encoded<'a>(
        &'a self, other: &'a IndexedKmerSet<MAX_LEN, E>,
    ) -> impl Iterator<Item = E::EncodedKmer> + 'a {
        assert_eq!(self.kmer_length(), other.kmer_length());

        self.vec
            .iter()
            .zip(&other.vec)
            .enumerate()
            .filter(|(_, (selfs, others))| **selfs && !**others)
            .map(|(index, _)| E::EncodedKmer::from_usize(index))
    }

    /// Returns an iterator over the decoded k-mers representing a set
    /// difference
    ///
    /// In other words, this returns the k-mers that are in `self` but not in
    /// `other`.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`IndexedKmerSet`]s are not
    /// equal.
    #[inline]
    pub fn difference_decoded<'a>(&'a self, other: &'a IndexedKmerSet<MAX_LEN, E>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
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
    /// Will panic if the `kmer_length` of the two [`IndexedKmerSet`]s are not
    /// equal.
    #[inline]
    pub fn intersection_encoded<'a>(
        &'a self, other: &'a IndexedKmerSet<MAX_LEN, E>,
    ) -> impl Iterator<Item = E::EncodedKmer> + 'a {
        assert_eq!(self.kmer_length(), other.kmer_length());

        self.vec
            .iter()
            .zip(&other.vec)
            .enumerate()
            .filter(|(_, (selfs, others))| **selfs && **others)
            .map(|(index, _)| E::EncodedKmer::from_usize(index))
    }

    /// Returns an iterator over the decoded k-mers representing a set
    /// intersection.
    ///
    /// In other words, this returns the k-mers that are in both `self` and
    /// `other`.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`IndexedKmerSet`]s are not
    /// equal.
    #[inline]
    pub fn intersection_decoded<'a>(&'a self, other: &'a IndexedKmerSet<MAX_LEN, E>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
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
    /// Will panic if the `kmer_length` of the two [`IndexedKmerSet`]s are not
    /// equal.
    #[inline]
    pub fn symmetric_difference_encoded<'a>(
        &'a self, other: &'a IndexedKmerSet<MAX_LEN, E>,
    ) -> impl Iterator<Item = E::EncodedKmer> + 'a {
        assert_eq!(self.kmer_length(), other.kmer_length());

        self.vec
            .iter()
            .zip(&other.vec)
            .enumerate()
            .filter(|(_, (selfs, others))| **selfs ^ **others)
            .map(|(index, _)| E::EncodedKmer::from_usize(index))
    }

    /// Returns an iterator over the decoded k-mers representing a set symmetric
    /// difference.
    ///
    /// In other words, this returns the k-mers that are in `self` or in `other`
    /// but not in both.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`IndexedKmerSet`]s are not
    /// equal.
    #[inline]
    pub fn symmetric_difference_decoded<'a>(
        &'a self, other: &'a IndexedKmerSet<MAX_LEN, E>,
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
    /// Will panic if the `kmer_length` of the two [`IndexedKmerSet`]s are not
    /// equal.
    #[inline]
    pub fn union_encoded<'a>(&'a self, other: &'a IndexedKmerSet<MAX_LEN, E>) -> impl Iterator<Item = E::EncodedKmer> + 'a {
        assert_eq!(self.kmer_length(), other.kmer_length());

        self.vec
            .iter()
            .zip(&other.vec)
            .enumerate()
            .filter(|(_, (selfs, others))| **selfs || **others)
            .map(|(index, _)| E::EncodedKmer::from_usize(index))
    }

    /// Returns an iterator over the decoded k-mers representing a set union.
    ///
    /// In other words, this returns the k-mers that are in `self` or `other`,
    /// without duplicates.
    ///
    /// ## Panics
    ///
    /// Will panic if the `kmer_length` of the two [`IndexedKmerSet`]s are not
    /// equal.
    #[inline]
    pub fn union_decoded<'a>(&'a self, other: &'a IndexedKmerSet<MAX_LEN, E>) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.union_encoded(other))
    }
}

impl<const MAX_LEN: usize, E> EncodedKmerCollection<MAX_LEN> for IndexedKmerSet<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    type Encoder = E;
    type EncodedKmer = E::EncodedKmer;

    #[inline]
    fn encoder(&self) -> &Self::Encoder {
        &self.encoder
    }
}

impl<const MAX_LEN: usize, E> FindKmersInSeq<MAX_LEN> for IndexedKmerSet<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    #[inline]
    fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, Self::Encoder>, {
        self.contains(kmer)
    }
}

impl<const MAX_LEN: usize, E> IntoIterator for IndexedKmerSet<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
{
    type Item = Kmer<MAX_LEN>;
    type IntoIter = IndexedKmerSetDecodedIntoIter<MAX_LEN, E>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            vec_into_iter: self.vec.into_iter().enumerate(),
            encoder:       self.encoder,
        }
    }
}

/// An iterator over an [`IndexedKmerSet`] yielding decoded k-mers. The iterator
/// consumes the original set.
pub struct IndexedKmerSetDecodedIntoIter<const MAX_LEN: usize, E>
where
    E: KmerEncoder<MAX_LEN>, {
    vec_into_iter: Enumerate<std::vec::IntoIter<bool>>,
    encoder:       E,
}

impl<const MAX_LEN: usize, E> Iterator for IndexedKmerSetDecodedIntoIter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
{
    type Item = Kmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Kmer<MAX_LEN>> {
        self.vec_into_iter
            .find_map(|(index, seen)| seen.then(|| self.encoder.decode_kmer(E::EncodedKmer::from_usize(index))))
    }
}

impl<const MAX_LEN: usize, E> Len for IndexedKmerSet<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN>,
{
    #[inline]
    fn is_empty(&self) -> bool {
        !self.vec.iter().any(|&seen| seen)
    }

    #[inline]
    fn len(&self) -> usize {
        self.vec.iter().filter(|&&seen| seen).count()
    }
}

#[cfg(test)]
mod tests {
    use crate::kmer::{FindKmersInSeq, IndexedKmerSet, encoders::two_bit::TwoBitKmerEncoder};

    #[test]
    fn finds_kmers_in_sequence() {
        let mut set = IndexedKmerSet::<3, TwoBitKmerEncoder<3>>::new(3).unwrap();
        set.insert_kmer(&b"ACG");

        assert_eq!(set.find_in_seq(b"TTACGAA"), Some(2..5));
        assert_eq!(set.find_in_seq_rev(b"ACGTTACG"), Some(5..8));
        assert_eq!(set.find_all_in_seq(b"ACGACG").collect::<Vec<_>>(), vec![0..3, 3..6]);
    }
}
