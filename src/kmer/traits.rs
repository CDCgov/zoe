//! Traits providing generic k-mer functionality to collections

use crate::{
    kmer::{
        KmerEncoder, SupportedKmerLen,
        encoders::three_bit::{ThreeBitEncodedKmer, ThreeBitKmerLen},
    },
    private::Sealed,
};
use std::{iter::Enumerate, ops::Range};

/// A helper trait for unifying encoded and decoded k-mers, so that functions
/// can accept either.
///
/// The trait provides a way to convert an encoded or decoded k-mer into an
/// encoded k-mer. It is implemented on anything implementing `AsRef<[u8]>` (a
/// decoded k-mer), as well as each potential type used as an encoded k-mer.
///
/// ## Parameters
///
/// - `MAX_LEN`: The maximum length that the [`KmerEncoder`] can support.
/// - `E`: The [`KmerEncoder`] used to encode the k-mer.
pub trait KmerEncode<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> {
    /// Ensures the provided k-mer is encoded. If the k-mer is already encoded,
    /// returns the input value.
    ///
    /// If the k-mer is decoded, then it must be of the correct length for the
    /// `encoder`.
    #[must_use]
    fn encode_kmer(&self, encoder: &E) -> E::EncodedKmer;
}

impl<const MAX_LEN: usize, E> KmerEncode<MAX_LEN, E> for ThreeBitEncodedKmer<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
    E: KmerEncoder<MAX_LEN, EncodedKmer = Self>,
{
    #[inline]
    fn encode_kmer(&self, _encoder: &E) -> E::EncodedKmer {
        *self
    }
}

impl<const MAX_LEN: usize, E> KmerEncode<MAX_LEN, E> for &ThreeBitEncodedKmer<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
    E: KmerEncoder<MAX_LEN, EncodedKmer = ThreeBitEncodedKmer<MAX_LEN>>,
{
    #[inline]
    fn encode_kmer(&self, _encoder: &E) -> E::EncodedKmer {
        **self
    }
}

impl<T: AsRef<[u8]>, const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> KmerEncode<MAX_LEN, E> for T {
    #[inline]
    fn encode_kmer(&self, encoder: &E) -> E::EncodedKmer {
        encoder.encode_kmer(self)
    }
}

/// A trait unifying encoded k-mer collections.
pub trait EncodedKmerCollection<const MAX_LEN: usize>: Sealed {
    /// The encoder used to encode the k-mers.
    type Encoder: KmerEncoder<MAX_LEN, EncodedKmer = Self::EncodedKmer>;
    /// The type used to store the encoded k-mer.
    type EncodedKmer: KmerEncode<MAX_LEN, Self::Encoder>;

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

/// A trait providing methods for searching a sequence for any of the k-mers in
/// a collection.
///
/// These methods are called on the collection. An extension trait for sequences
/// with similar methods is provided by [`FindKmers`], in which case the
/// collection is passed as the argument instead.
///
/// ## Parameters
///
/// `MAX_LEN`: The maximum k-mer length supported by the collection
pub trait FindKmersInSeq<const MAX_LEN: usize>: EncodedKmerCollection<MAX_LEN> + Sealed {
    /// Checks whether the collection contains a k-mer.
    ///
    /// The k-mers can be either encoded or decoded (in which case it is encoded
    /// before checking). If it is encoded, it must have been generated using
    /// the [`KmerEncoder`] associated with this collection. If it is decoded,
    /// it must be of length `self.kmer_length()`.
    #[must_use]
    fn contains<K>(&self, kmer: &K) -> bool
    where
        K: KmerEncode<MAX_LEN, Self::Encoder>;

    /// Returns the indices of the leftmost occurrence of any of the k-mers in
    /// this collection within the provided sequence.
    ///
    /// If no occurrence is found, then `None` is returned. To find all matches
    /// rather than just the first, use [`find_all_in_seq`].
    ///
    /// [`find_all_in_seq`]: FindKmersInSeq::find_all_in_seq
    #[must_use]
    fn find_in_seq(&self, seq: impl AsRef<[u8]>) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder().iter_from_sequence(&seq).enumerate() {
            if self.contains(&kmer) {
                return Some(i..i + self.kmer_length());
            }
        }
        None
    }

    /// Returns the indices of the rightmost occurrence of any of the k-mers in
    /// this collection within the provided sequence.
    ///
    /// If no occurrence is found, then `None` is returned. To find all matches
    /// rather than just the last, use [`find_all_in_seq_rev`].
    ///
    /// [`find_all_in_seq_rev`]: FindKmersInSeq::find_all_in_seq_rev
    #[must_use]
    fn find_in_seq_rev(&self, seq: impl AsRef<[u8]>) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder().iter_from_sequence_rev(&seq).enumerate() {
            if self.contains(&kmer) {
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
    /// [`find_in_seq`]: FindKmersInSeq::find_in_seq
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
    /// [`find_in_seq_rev`]: FindKmersInSeq::find_in_seq_rev
    #[inline]
    #[must_use]
    fn find_all_in_seq_rev<'a>(&'a self, seq: &'a (impl AsRef<[u8]> + ?Sized)) -> KmerMatchesInSeqRev<'a, MAX_LEN, Self> {
        let seq = seq.as_ref();
        KmerMatchesInSeqRev::new(self, seq)
    }
}

/// An iterator over the matches between the k-mers in a sequence and those in a
/// collection.
///
/// See [`find_all_in_seq`] for more details.
///
/// [`find_all_in_seq`]: FindKmersInSeq::find_all_in_seq
pub struct KmerMatchesInSeq<'a, const MAX_LEN: usize, T: FindKmersInSeq<MAX_LEN> + ?Sized> {
    collection:   &'a T,
    seq_iterator: Enumerate<<<T as EncodedKmerCollection<MAX_LEN>>::Encoder as KmerEncoder<MAX_LEN>>::SeqIter<'a>>,
}

impl<'a, const MAX_LEN: usize, T: FindKmersInSeq<MAX_LEN> + ?Sized> KmerMatchesInSeq<'a, MAX_LEN, T> {
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
    T: FindKmersInSeq<MAX_LEN>,
{
    type Item = Range<usize>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        for (i, kmer) in &mut self.seq_iterator {
            if self.collection.contains(&kmer) {
                return Some(i..i + self.collection.kmer_length());
            }
        }
        None
    }
}

/// An iterator over the matches between the k-mers in a sequence and those in a
/// collection, in reverse order.
///
/// See [`find_all_in_seq_rev`] for more details.
///
/// [`find_all_in_seq_rev`]: FindKmersInSeq::find_all_in_seq_rev
pub struct KmerMatchesInSeqRev<'a, const MAX_LEN: usize, T: FindKmersInSeq<MAX_LEN> + ?Sized> {
    collection:   &'a T,
    seq_iterator: Enumerate<<<T as EncodedKmerCollection<MAX_LEN>>::Encoder as KmerEncoder<MAX_LEN>>::SeqIterRev<'a>>,
    seq_len:      usize,
}

impl<'a, const MAX_LEN: usize, T: FindKmersInSeq<MAX_LEN> + ?Sized> KmerMatchesInSeqRev<'a, MAX_LEN, T> {
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
    T: FindKmersInSeq<MAX_LEN>,
{
    type Item = Range<usize>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        for (i, kmer) in &mut self.seq_iterator {
            if self.collection.contains(&kmer) {
                let end = self.seq_len - i;
                return Some(end - self.collection.kmer_length()..end);
            }
        }
        None
    }
}

/// Provides methods for searching for a k-mer in a collection.
///
/// Corresponding methods exist for [`RangeSearch`]. See [Restricting the search
/// range](crate::search#restricting-the-search-range) for more details.
///
/// These methods are called on the sequence. To instead call them on the
/// collection (and pass the sequence as an argument), see [`FindKmersInSeq`].
///
/// [`RangeSearch`]: crate::search::RangeSearch
pub trait FindKmers<const MAX_LEN: usize>: AsRef<[u8]> {
    /// Returns the indices of the leftmost occurrence of any of the k-mers in
    /// the passed collection.
    ///
    /// If no occurrence is found, then `None` is returned. To find all matches
    /// rather than just the first, use [`find_all_kmers`].
    ///
    /// [`find_all_kmers`]: FindKmers::find_all_kmers
    #[inline]
    #[must_use]
    fn find_kmers<T: FindKmersInSeq<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        kmers.find_in_seq(self)
    }

    /// Returns the indices of the rightmost occurrence of any of the k-mers in
    /// the passed collection.
    ///
    /// If no occurrence is found, then `None` is returned. To find all matches
    /// rather than just the last, use [`find_all_kmers`].
    ///
    /// [`find_all_kmers`]: FindKmers::find_all_kmers
    #[inline]
    #[must_use]
    fn find_kmers_rev<T: FindKmersInSeq<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        kmers.find_in_seq_rev(self)
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in the sequence and those in the passed collection.
    ///
    /// Consider using [`find_kmers`] if only a single match is needed.
    ///
    /// [`find_kmers`]: FindKmers::find_kmers
    #[inline]
    #[must_use]
    fn find_all_kmers<'a, T: FindKmersInSeq<MAX_LEN>>(&'a self, kmers: &'a T) -> impl Iterator<Item = Range<usize>> + 'a {
        kmers.find_all_in_seq(self)
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in the sequence and those in the passed collection, in reverse order.
    ///
    /// Consider using [`find_kmers_rev`] if only a single match is needed.
    ///
    /// [`find_kmers_rev`]: FindKmers::find_kmers_rev
    #[inline]
    #[must_use]
    fn find_all_kmers_rev<'a, T: FindKmersInSeq<MAX_LEN>>(
        &'a self, kmers: &'a T,
    ) -> impl Iterator<Item = Range<usize>> + 'a {
        kmers.find_all_in_seq_rev(self)
    }
}

impl<const MAX_LEN: usize, Q: AsRef<[u8]>> FindKmers<MAX_LEN> for Q {}
