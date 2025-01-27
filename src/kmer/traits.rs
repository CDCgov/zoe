use crate::search::RangeSearch;

use super::{Kmer, KmerEncoder, SupportedKmerLen, ThreeBitEncodedKmer, ThreeBitKmerLen};
use std::ops::Range;

// TODO: Maybe should consume self
/// Trait to support converting an iterator of k-mers into a collection. This
/// trait provides a way to convert an encoded or decoded k-mer into an encoded
/// k-mer. Hence, it is implemented on [`Kmer`], as well as each potential type
/// used as an encoded k-mer.
pub trait KmerEncode<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> {
    fn encode_kmer(&self, encoder: &E) -> E::EncodedKmer;
}

impl<const MAX_LEN: usize, E> KmerEncode<MAX_LEN, E> for ThreeBitEncodedKmer<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
    E: KmerEncoder<MAX_LEN, EncodedKmer = Self>,
{
    fn encode_kmer(&self, _encoder: &E) -> E::EncodedKmer {
        *self
    }
}

impl<const MAX_LEN: usize, E> KmerEncode<MAX_LEN, E> for &ThreeBitEncodedKmer<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
    E: KmerEncoder<MAX_LEN, EncodedKmer = ThreeBitEncodedKmer<MAX_LEN>>,
{
    fn encode_kmer(&self, _encoder: &E) -> E::EncodedKmer {
        **self
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> KmerEncode<MAX_LEN, E> for Kmer<MAX_LEN> {
    fn encode_kmer(&self, encoder: &E) -> E::EncodedKmer {
        encoder.encode_kmer(self)
    }
}

/// Trait for encoded k-mer collections.
pub trait EncodedKmerCollection<const MAX_LEN: usize> {
    type Encoder: KmerEncoder<MAX_LEN, EncodedKmer = Self::EncodedKmer>;
    type EncodedKmer;

    /// Gets the encoder associated with this collection.
    #[must_use]
    fn encoder(&self) -> &Self::Encoder;

    /// Gets the length of the k-mers being stored in this collection.
    #[inline]
    #[must_use]
    fn kmer_length(&self) -> usize {
        self.encoder().kmer_length()
    }
}

/// Trait for k-mer collections with a `contains` method (i.e., that can act as
/// a set).
pub trait KmerCollectionContains<const MAX_LEN: usize>: EncodedKmerCollection<MAX_LEN> {
    /// Return whether an already encoded k-mer is present in the collection.
    /// The encoded k-mer must have been generated using the [`KmerEncoder`]
    /// associated with this collection.
    #[must_use]
    fn contains_encoded(&self, kmer: Self::EncodedKmer) -> bool;

    /// Return whether a k-mer is present in the collection. The bases and k-mer
    /// length are assumed to be valid for the [`KmerEncoder`] associated with
    /// this collection. Consider [`contains_checked`] when it is not known
    /// whether the bases and k-mer length will be valid.
    ///
    /// [`contains_checked`]: KmerCollectionContains::contains_checked
    #[inline]
    #[must_use]
    fn contains(&self, kmer: impl AsRef<[u8]>) -> bool {
        self.contains_encoded(self.encoder().encode_kmer(kmer))
    }

    /// Return whether a k-mer is present in the collection. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// collection, then `None` is returned.
    #[inline]
    #[must_use]
    fn contains_checked(&self, kmer: impl AsRef<[u8]>) -> Option<bool> {
        Some(self.contains_encoded(self.encoder().encode_kmer_checked(kmer)?))
    }

    /// Return the indices of the leftmost occurrence of any of the k-mers in
    /// this collection within a provided sequence. The bases in the sequence
    /// must be valid for the [`KmerEncoder`] associated with this collection.
    /// If no occurrence is found, then `None` is returned.
    #[must_use]
    fn find_in_seq(&self, seq: impl AsRef<[u8]>) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder().iter_from_sequence(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                return Some(i..i + self.kmer_length());
            }
        }
        None
    }

    /// Return the indices of the rightmost occurrence of any of the k-mers in
    /// this collection within a provided sequence. The bases in the sequence
    /// must be valid for the [`KmerEncoder`] associated with this collection.
    /// If no occurrence is found, then `None` is returned.
    #[must_use]
    fn find_in_seq_rev(&self, seq: impl AsRef<[u8]>) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder().iter_from_sequence_rev(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                let end = seq.as_ref().len() - i;
                return Some(end - self.kmer_length()..end);
            }
        }
        None
    }
}

/// Provides methods for searching for a kmer in a collection.
///
/// This trait is also compatible with [`RangeSearch`]. See [Restricting the
/// search range](crate::search#restricting-the-search-range) for more details.
pub trait FindKmers<const MAX_LEN: usize> {
    /// Return the indices of the leftmost occurrence of any of the k-mers in a
    /// collection. The bases in the sequence must be valid for the
    /// [`KmerEncoder`] associated with the collection. If no occurrence is
    /// found, then `None` is returned.
    #[must_use]
    fn find_kmers<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>>;

    /// Return the indices of the rightmost occurrence of any of the k-mers in a
    /// collection. The bases in the sequence must be valid for the
    /// [`KmerEncoder`] associated with the collection. If no occurrence is
    /// found, then `None` is returned.
    #[must_use]
    fn find_kmers_rev<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>>;
}

impl<const MAX_LEN: usize, Q: AsRef<[u8]>> FindKmers<MAX_LEN> for Q {
    #[inline]
    fn find_kmers<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        kmers.find_in_seq(self.as_ref())
    }

    #[inline]
    fn find_kmers_rev<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        kmers.find_in_seq_rev(self.as_ref())
    }
}

impl<const MAX_LEN: usize> FindKmers<MAX_LEN> for RangeSearch<'_> {
    #[inline]
    fn find_kmers<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        self.slice.find_kmers(kmers).map(|r| self.adjust_to_context(&r))
    }

    #[inline]
    fn find_kmers_rev<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        self.slice.find_kmers_rev(kmers).map(|r| self.adjust_to_context(&r))
    }
}
