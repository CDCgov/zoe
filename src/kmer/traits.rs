use crate::private::Sealed;

use super::{Kmer, KmerEncoder, SupportedKmerLen, ThreeBitEncodedKmer, ThreeBitKmerLen};
use std::{iter::Enumerate, ops::Range};

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
pub trait EncodedKmerCollection<const MAX_LEN: usize>: Sealed {
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
pub trait KmerCollectionContains<const MAX_LEN: usize>: EncodedKmerCollection<MAX_LEN> + Sealed {
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
    ///
    /// To find all matches rather than just the first, use [`find_all_in_seq`].
    ///
    /// [`find_all_in_seq`]: KmerCollectionContains::find_all_in_seq
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
    ///
    /// To find all matches rather than just the last, use [`find_all_in_seq`].
    ///
    /// [`find_all_in_seq`]: KmerCollectionContains::find_all_in_seq
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

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in a sequence and those in this collection.
    ///
    /// Consider using [`find_in_seq`] if only a single match is needed.
    ///
    /// [`find_in_seq`]: KmerCollectionContains::find_in_seq
    #[inline]
    #[must_use]
    fn find_all_in_seq<'a>(&'a self, seq: &'a (impl AsRef<[u8]> + ?Sized)) -> KmerMatchesInSeq<'a, MAX_LEN, Self> {
        let seq = seq.as_ref();
        KmerMatchesInSeq::new(self, seq)
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in a sequence and those in this collection, in reverse order.
    ///
    /// Consider using [`find_in_seq_rev`] if only a single match is needed.
    ///
    /// [`find_in_seq_rev`]: KmerCollectionContains::find_in_seq_rev
    #[inline]
    #[must_use]
    fn find_all_in_seq_rev<'a>(&'a self, seq: &'a (impl AsRef<[u8]> + ?Sized)) -> KmerMatchesInSeqRev<'a, MAX_LEN, Self> {
        let seq = seq.as_ref();
        KmerMatchesInSeqRev::new(self, seq)
    }
}

/// An iterator over the matches between the k-mers in a sequence and those in a
/// collection. See [`find_all_in_seq`] for more details.
///
/// [`find_all_in_seq`]: KmerCollectionContains::find_all_in_seq
pub struct KmerMatchesInSeq<'a, const MAX_LEN: usize, T: KmerCollectionContains<MAX_LEN> + ?Sized> {
    collection:   &'a T,
    seq_iterator: Enumerate<<<T as EncodedKmerCollection<MAX_LEN>>::Encoder as KmerEncoder<MAX_LEN>>::SeqIter<'a>>,
}

impl<'a, const MAX_LEN: usize, T: KmerCollectionContains<MAX_LEN> + ?Sized> KmerMatchesInSeq<'a, MAX_LEN, T> {
    /// Create a new instance of the iterator.
    #[inline]
    fn new(collection: &'a T, seq: &'a [u8]) -> Self {
        Self {
            collection,
            seq_iterator: collection.encoder().iter_from_sequence(seq).enumerate(),
        }
    }
}

impl<const MAX_LEN: usize, T> Iterator for KmerMatchesInSeq<'_, MAX_LEN, T>
where
    T: KmerCollectionContains<MAX_LEN>,
{
    type Item = Range<usize>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        for (i, kmer) in &mut self.seq_iterator {
            if self.collection.contains_encoded(kmer) {
                return Some(i..i + self.collection.kmer_length());
            }
        }
        None
    }
}

/// An iterator over the matches between the k-mers in a sequence and those in a
/// collection, in reverse order. See [`find_all_in_seq_rev`] for more details.
///
/// [`find_all_in_seq_rev`]: KmerCollectionContains::find_all_in_seq_rev
pub struct KmerMatchesInSeqRev<'a, const MAX_LEN: usize, T: KmerCollectionContains<MAX_LEN> + ?Sized> {
    collection:   &'a T,
    seq_iterator: Enumerate<<<T as EncodedKmerCollection<MAX_LEN>>::Encoder as KmerEncoder<MAX_LEN>>::SeqIterRev<'a>>,
    seq_len:      usize,
}

impl<'a, const MAX_LEN: usize, T: KmerCollectionContains<MAX_LEN> + ?Sized> KmerMatchesInSeqRev<'a, MAX_LEN, T> {
    /// Create a new instance of the iterator.
    #[inline]
    fn new(collection: &'a T, seq: &'a [u8]) -> Self {
        Self {
            collection,
            seq_iterator: collection.encoder().iter_from_sequence_rev(seq).enumerate(),
            seq_len: seq.len(),
        }
    }
}

impl<const MAX_LEN: usize, T> Iterator for KmerMatchesInSeqRev<'_, MAX_LEN, T>
where
    T: KmerCollectionContains<MAX_LEN>,
{
    type Item = Range<usize>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        for (i, kmer) in &mut self.seq_iterator {
            if self.collection.contains_encoded(kmer) {
                let end = self.seq_len - i;
                return Some(end - self.collection.kmer_length()..end);
            }
        }
        None
    }
}

/// Provides methods for searching for a kmer in a collection.
///
/// Corresponding methods exist for [`RangeSearch`]. See [Restricting the search
/// range](crate::search#restricting-the-search-range) for more details.
///
/// [`RangeSearch`]: crate::search::RangeSearch
pub trait FindKmers<const MAX_LEN: usize>: AsRef<[u8]> {
    /// Return the indices of the leftmost occurrence of any of the k-mers in a
    /// collection. The bases in the sequence must be valid for the
    /// [`KmerEncoder`] associated with the collection. If no occurrence is
    /// found, then `None` is returned.
    ///
    /// To find all matches rather than just the first, use [`find_all_kmers`].
    ///
    /// [`find_all_kmers`]: FindKmers::find_all_kmers
    #[inline]
    #[must_use]
    fn find_kmers<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        kmers.find_in_seq(self)
    }

    /// Return the indices of the rightmost occurrence of any of the k-mers in a
    /// collection. The bases in the sequence must be valid for the
    /// [`KmerEncoder`] associated with the collection. If no occurrence is
    /// found, then `None` is returned.
    ///
    /// To find all matches rather than just the last, use [`find_all_kmers`].
    ///
    /// [`find_all_kmers`]: FindKmers::find_all_kmers
    #[inline]
    #[must_use]
    fn find_kmers_rev<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        kmers.find_in_seq_rev(self)
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in the sequence and those in a collection.
    ///
    /// Consider using [`find_kmers`] if only a single match is needed.
    ///
    /// [`find_kmers`]: FindKmers::find_kmers
    #[inline]
    #[must_use]
    fn find_all_kmers<'a, T: KmerCollectionContains<MAX_LEN>>(
        &'a self, kmers: &'a T,
    ) -> impl Iterator<Item = Range<usize>> + 'a {
        kmers.find_all_in_seq(self)
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in the sequence and those in a collection, in reverse order.
    ///
    /// Consider using [`find_kmers_rev`] if only a single match is needed.
    ///
    /// [`find_kmers_rev`]: FindKmers::find_kmers_rev
    #[inline]
    #[must_use]
    fn find_all_kmers_rev<'a, T: KmerCollectionContains<MAX_LEN>>(
        &'a self, kmers: &'a T,
    ) -> impl Iterator<Item = Range<usize>> + 'a {
        kmers.find_all_in_seq_rev(self)
    }
}

impl<const MAX_LEN: usize, Q: AsRef<[u8]>> FindKmers<MAX_LEN> for Q {}
