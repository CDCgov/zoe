//! Defines [`IndexedKmerCounter`], which holds encoded k-mers and their counts,
//! using the encoded k-mers to index into a `Vec`.

use crate::{
    data::views::Len,
    kmer::{EncodedKmerCollection, FindKmersInSeq, GetVariants, Kmer, KmerEncode, KmerEncoder, KmerError, KmerIndex},
};
use std::{
    iter::Enumerate,
    ops::{Index, IndexMut},
};

/// A collection of k-mers and their counts, with encoded k-mers used to index
/// the collection and access the counts. This works similarly to
/// [`KmerCounter`], which wraps a [`HashMap`] of k-mers.
///
/// K-mers can be tallied into an [`IndexedKmerCounter`] multiple ways:
///
/// - A single k-mer can be tallied with [`tally_kmer`]
/// - A k-mer as well as similar k-mers (up to `N` mismatches) can be tallied
///   with [`tally_kmer_with_variants`]
/// - Multiple k-mers from an iterator can be tallied with [`tally_from_iter`]
/// - Overlapping k-mers from a sequence can be tallied with
///   [`tally_from_sequence`]
/// - Overlapping k-mers with mismatches from a sequence can be tallied with
///   [`tally_from_sequence_with_variants`]
///
/// After a k-mer counter is populated, it can be used in multiple ways:
///
/// - Check for k-mer with [`contains`], or get its count with [`get`]
/// - Iteration: [`iter_encoded`] and [`iter_decoded`] provide the k-mers and
///   counts, while [`keys_encoded`] and [`keys_decoded`] provide just the
///   k-mers
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
/// [`KmerCounter`]: crate::kmer::KmerCounter
/// [`HashMap`]: std::collections::HashMap
/// [`tally_kmer`]: IndexedKmerCounter::tally_kmer
/// [`tally_kmer_with_variants`]: IndexedKmerCounter::tally_kmer_with_variants
/// [`tally_from_iter`]: IndexedKmerCounter::tally_from_iter
/// [`tally_from_sequence`]: IndexedKmerCounter::tally_from_sequence
/// [`tally_from_sequence_with_variants`]:
///     IndexedKmerCounter::tally_from_sequence_with_variants
/// [`contains`]: IndexedKmerCounter::contains
/// [`get`]: IndexedKmerCounter::get
/// [`iter_encoded`]: IndexedKmerCounter::iter_encoded
/// [`iter_decoded`]: IndexedKmerCounter::iter_decoded
/// [`keys_encoded`]: IndexedKmerCounter::keys_encoded
/// [`keys_decoded`]: IndexedKmerCounter::keys_decoded
/// [`FindKmers`]: crate::kmer::FindKmers
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct IndexedKmerCounter<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> {
    /// The vec, indexed by encoded k-mers, storing their respective counts.
    vec:     Box<[usize]>,
    /// The encoder used to encode the k-mers.
    encoder: E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> IndexedKmerCounter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
{
    /// Creates a new [`IndexedKmerCounter`] with the specified k-mer length.
    ///
    /// ## Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn new(kmer_length: usize) -> Result<Self, KmerError> {
        Ok(Self {
            vec:     vec![0; E::EncodedKmer::max_index_for_length(kmer_length)].into_boxed_slice(),
            encoder: E::new(kmer_length)?,
        })
    }

    /// Tallies the k-mer in the [`IndexedKmerCounter`].
    ///
    /// The k-mer can be either encoded or decoded (in which case it gets
    /// automatically encoded). If it is encoded, it must have been generated
    /// using the [`KmerEncoder`] associated with this [`IndexedKmerCounter`].
    /// If it is decoded, it must be of length [`Self::kmer_length`].
    #[inline]
    pub fn tally_kmer<K>(&mut self, kmer: &K)
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.vec[kmer.encode_kmer(&self.encoder).as_usize()] += 1;
    }

    /// Tallies all k-mers into the [`IndexedKmerCounter`] with at most `N`
    /// mismatches compared to the provided k-mer.
    ///
    /// The original k-mer is also tallied. The k-mer can be either encoded or
    /// decoded (in which case it gets automatically encoded). If it is encoded,
    /// it must have been generated using the [`KmerEncoder`] associated with
    /// this [`IndexedKmerCounter`]. If it is decoded, it must be of length
    /// [`Self::kmer_length`].
    #[inline]
    pub fn tally_kmer_with_variants<const N: usize>(&mut self, kmer: &impl KmerEncode<MAX_LEN, E>)
    where
        E: GetVariants<N, MAX_LEN>, {
        self.encoder
            .get_variants::<N>(kmer.encode_kmer(&self.encoder))
            .for_each(|variant| self.tally_kmer(&variant));
    }

    /// Tallies k-mers from an iterator into the [`IndexedKmerCounter`].
    ///
    /// The k-mer can be either encoded or decoded (in which case it is encoded
    /// for indexing). If it is encoded, it must have been generated using the
    /// [`KmerEncoder`] associated with this [`IndexedKmerCounter`]. If it is
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
    pub fn tally_from_iter<I: IntoIterator<Item: KmerEncode<MAX_LEN, E>>>(&mut self, iter: I) {
        iter.into_iter().for_each(|kmer| self.tally_kmer(&kmer));
    }

    /// Tallies all overlapping k-mers from a sequence in the
    /// [`IndexedKmerCounter`].
    #[inline]
    pub fn tally_from_sequence(&mut self, seq: impl AsRef<[u8]>) {
        self.encoder.iter_from_sequence(&seq).for_each(|kmer| self.tally_kmer(&kmer));
    }

    /// Tallies all k-mers from a sequence into the [`IndexedKmerCounter`], in
    /// addition to all k-mers with up to `N` mismatches from those in the
    /// sequence.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::kmer::{IndexedKmerCounter, encoders::three_bit::ThreeBitKmerEncoder};
    /// let mut counter = IndexedKmerCounter::<4, ThreeBitKmerEncoder<4>>::new(4).unwrap();
    /// let seq = b"GATAGGGGATTGT";
    /// counter.tally_from_sequence_with_variants::<2>(seq);
    /// ```
    #[inline]
    pub fn tally_from_sequence_with_variants<const N: usize>(&mut self, seq: impl AsRef<[u8]>)
    where
        E: GetVariants<N, MAX_LEN>, {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.tally_kmer_with_variants::<N>(&encoded_kmer);
        }
    }

    /// Checks whether the [`IndexedKmerCounter`] contains a k-mer (with a
    /// nonzero count).
    ///
    /// The k-mer can be either encoded or decoded (in which case it gets
    /// automatically encoded). If it is encoded, it must have been generated
    /// using the [`KmerEncoder`] associated with this [`IndexedKmerCounter`].
    /// If it is decoded, it must be of length [`Self::kmer_length`].
    #[inline]
    #[must_use]
    pub fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.vec[kmer.encode_kmer(&self.encoder).as_usize()] > 0
    }

    /// Gets the count of a k-mer.
    ///
    /// If the k-mer is not present in the counter, then `0` is returned. The
    /// k-mer can be either encoded or decoded (in which case it gets
    /// automatically encoded). If it is encoded, it must have been generated
    /// using the [`KmerEncoder`] associated with this [`IndexedKmerCounter`].
    /// If it is decoded, it must be of length [`Self::kmer_length`].
    #[inline]
    #[must_use]
    pub fn get<K>(&self, kmer: &K) -> usize
    where
        K: KmerEncode<MAX_LEN, E>, {
        self.vec[kmer.encode_kmer(&self.encoder).as_usize()]
    }

    /// Returns an iterator over the encoded k-mers and their counts.
    ///
    /// K-mers must have a nonzero count to be included.
    #[inline]
    pub fn iter_encoded(&self) -> impl Iterator<Item = (E::EncodedKmer, &usize)> {
        self.vec
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .map(|(index, count)| (E::EncodedKmer::from_usize(index), count))
    }

    /// Returns an iterator over the decoded k-mers and their counts, if the
    /// count is nonzero.
    ///
    /// K-mers must have a nonzero count to be included.
    #[inline]
    pub fn iter_decoded(&self) -> impl Iterator<Item = (Kmer<MAX_LEN>, &usize)> {
        self.iter_encoded()
            .map(|(encoded_kmer, count)| (self.encoder.decode_kmer(encoded_kmer), count))
    }

    /// Returns an iterator of decoded k-mers for k-mers in the
    /// [`IndexedKmerCounter`] for k-mers with nonzero counts.
    #[inline]
    pub fn keys_encoded(&self) -> impl Iterator<Item = E::EncodedKmer> {
        self.vec
            .iter()
            .enumerate()
            .filter(|(_index, count)| **count > 0)
            .map(|(index, _count)| E::EncodedKmer::from_usize(index))
    }

    /// Returns an iterator of decoded k-mers for k-mers in the
    /// [`IndexedKmerCounter`] for k-mers with nonzero counts.
    #[inline]
    pub fn keys_decoded(&self) -> impl Iterator<Item = Kmer<MAX_LEN>> {
        self.encoder.decode_iter(self.keys_encoded())
    }
}

impl<const MAX_LEN: usize, E> Index<E::EncodedKmer> for IndexedKmerCounter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
{
    type Output = usize;

    #[inline]
    fn index(&self, index: E::EncodedKmer) -> &Self::Output {
        &self.vec[index.as_usize()]
    }
}

// This implementation is okay for IndexedKmerCounter since a count of zero
// means it isn't present in the counter, but this same implementation for
// KmerCounter is not allowed and could lead to invalid state
impl<const MAX_LEN: usize, E> IndexMut<E::EncodedKmer> for IndexedKmerCounter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
{
    #[inline]
    fn index_mut(&mut self, index: E::EncodedKmer) -> &mut Self::Output {
        &mut self.vec[index.as_usize()]
    }
}

impl<const MAX_LEN: usize, E> EncodedKmerCollection<MAX_LEN> for IndexedKmerCounter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN>,
{
    type Encoder = E;
    type EncodedKmer = E::EncodedKmer;

    #[inline]
    fn encoder(&self) -> &Self::Encoder {
        &self.encoder
    }
}

impl<const MAX_LEN: usize, E> FindKmersInSeq<MAX_LEN> for IndexedKmerCounter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
{
    #[inline]
    fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, Self::Encoder>, {
        self.contains(kmer)
    }
}

/// An iterator over a [`IndexedKmerCounter`] yielding decoded k-mers and their
/// counts.
///
/// The iterator consumes the original counter.
pub struct IndexedKmerCounterDecodedIntoIter<const MAX_LEN: usize, E> {
    vec_into_iter: Enumerate<std::vec::IntoIter<usize>>,
    encoder:       E,
}

impl<const MAX_LEN: usize, E> Iterator for IndexedKmerCounterDecodedIntoIter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
{
    type Item = (Kmer<MAX_LEN>, usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.vec_into_iter.find_map(|(index, count)| {
            (count > 0).then(|| (self.encoder.decode_kmer(E::EncodedKmer::from_usize(index)), count))
        })
    }
}

impl<const MAX_LEN: usize, E> IntoIterator for IndexedKmerCounter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN, EncodedKmer: KmerIndex>,
{
    type Item = (Kmer<MAX_LEN>, usize);
    type IntoIter = IndexedKmerCounterDecodedIntoIter<MAX_LEN, E>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            vec_into_iter: self.vec.into_iter().enumerate(),
            encoder:       self.encoder,
        }
    }
}

impl<const MAX_LEN: usize, E> Len for IndexedKmerCounter<MAX_LEN, E>
where
    E: KmerEncoder<MAX_LEN>,
{
    #[inline]
    fn is_empty(&self) -> bool {
        !self.vec.iter().any(|&count| count > 0)
    }

    #[inline]
    fn len(&self) -> usize {
        self.vec.iter().filter(|&&seen| seen > 0).count()
    }
}

#[cfg(test)]
mod tests {
    use crate::kmer::{FindKmersInSeq, IndexedKmerCounter, encoders::two_bit::TwoBitKmerEncoder};

    #[test]
    fn test_kmer_tally() {
        let mut counter = IndexedKmerCounter::<3, TwoBitKmerEncoder<3>>::new(3).unwrap();
        counter.tally_from_sequence(b"ACGACG");

        let counts = counter
            .iter_decoded()
            .filter(|(_, count)| **count > 0)
            .map(|(kmer, count)| (kmer.as_ref().to_vec(), *count))
            .collect::<Vec<_>>();

        assert_eq!(counts, vec![(b"ACG".to_vec(), 2), (b"CGA".to_vec(), 1), (b"GAC".to_vec(), 1)]);
    }

    #[test]
    fn iterators_include_only_counted_kmers() {
        let mut counter = IndexedKmerCounter::<3, TwoBitKmerEncoder<3>>::new(3).unwrap();
        counter.tally_from_sequence(b"ACGACG");

        let keys = counter.keys_decoded().map(|kmer| kmer.as_ref().to_vec()).collect::<Vec<_>>();
        assert_eq!(keys, vec![b"ACG", b"CGA", b"GAC"]);

        let counts = counter
            .into_iter()
            .map(|(kmer, count)| (kmer.as_ref().to_vec(), count))
            .collect::<Vec<_>>();
        assert_eq!(counts, vec![(b"ACG".to_vec(), 2), (b"CGA".to_vec(), 1), (b"GAC".to_vec(), 1)]);
    }

    #[test]
    fn finds_counted_kmers_in_a_sequence() {
        let mut counter = IndexedKmerCounter::<3, TwoBitKmerEncoder<3>>::new(3).unwrap();
        counter.tally_kmer(&b"ACG");

        assert_eq!(counter.find_all_in_seq(b"ACGACG").collect::<Vec<_>>(), vec![0..3, 3..6]);
    }
}
