use crate::{
    data::view_traits::{IndexAdjustable, SliceRange},
    kmer::{FindKmers, KmerCollectionContains},
};
use std::ops::{Bound, Index, Range};

/// A subsequence along with its starting index. This is used for restricting
/// the search range when performing a search. See [Restricting the search
/// range](crate::search#restricting-the-search-range) for more details.
///
/// It is compatible with:
/// * Performing string search via the [`ByteSubstring`] trait
/// * Searching for k-mers using methods similar to those in [`FindKmers`]
///
/// <div class="warning note">
///
/// **Note**
///
/// All search methods called on a [`RangeSearch`] struct return indices with
/// respect to the original sequence, not the subsequence.
///
/// </div>
///
/// [`ByteSubstring`]: crate::search::substring::ByteSubstring
/// [`FindKmers`]: crate::kmer::FindKmers
pub struct RangeSearch<'a> {
    pub(crate) slice:       &'a [u8],
    pub(crate) starting_at: usize,
}

impl<'a> RangeSearch<'a> {
    /// Create a new [`RangeSearch`] from a sequence and a range.
    pub(crate) fn new<Q: AsRef<[u8]> + 'a + ?Sized, R: SliceRange>(sequence: &'a Q, range: R) -> Self {
        let slice = sequence.as_ref().index(range.clone());
        let starting_at = match range.start_bound() {
            Bound::Included(&start) => start,
            Bound::Unbounded => 0,
            Bound::Excluded(&start) => start.saturating_add(1),
        };
        Self { slice, starting_at }
    }

    /// Given an index/range in the frame of reference of `slice`, adjust it to
    /// be in the original frame of reference.
    #[inline]
    pub(crate) fn adjust_to_context<I: IndexAdjustable>(&self, index: &I) -> I {
        index.add(self.starting_at)
    }

    // FindKmers methods. These need to be distinct from the FindKmers trait
    // because the methods need to consume the RangeSearch struct, rather than
    // take them by reference. Otherwise, a bunch of useless let bindings are
    // required.

    /// Return the indices of the leftmost occurrence of any of the k-mers in a
    /// collection. The bases in the sequence must be valid for the
    /// [`KmerEncoder`] associated with the collection. If no occurrence is
    /// found, then `None` is returned.
    ///
    /// To find all matches rather than just the first, use [`find_all_kmers`].
    /// A similar unrestricted version of this method exists in [`FindKmers`].
    ///
    /// [`find_all_kmers`]: FindKmers::find_all_kmers
    /// [`KmerEncoder`]: crate::kmer::KmerEncoder
    #[inline]
    pub fn find_kmers<const MAX_LEN: usize, T>(self, kmers: &T) -> Option<Range<usize>>
    where
        T: KmerCollectionContains<MAX_LEN>, {
        self.slice.find_kmers(kmers).map(|r| self.adjust_to_context(&r))
    }

    /// Return the indices of the rightmost occurrence of any of the k-mers in a
    /// collection. The bases in the sequence must be valid for the
    /// [`KmerEncoder`] associated with the collection. If no occurrence is
    /// found, then `None` is returned.
    ///
    /// To find all matches rather than just the last, use [`find_all_kmers`]. A
    /// similar unrestricted version of this method exists in [`FindKmers`].
    ///
    /// [`find_all_kmers`]: FindKmers::find_all_kmers
    /// [`KmerEncoder`]: crate::kmer::KmerEncoder
    #[inline]
    pub fn find_kmers_rev<const MAX_LEN: usize, T>(self, kmers: &T) -> Option<Range<usize>>
    where
        T: KmerCollectionContains<MAX_LEN>, {
        self.slice.find_kmers_rev(kmers).map(|r| self.adjust_to_context(&r))
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in the sequence and those in a collection.
    ///
    /// Consider using [`find_kmers`] if only a single match is needed. A
    /// similar unrestricted version of this method exists in [`FindKmers`].
    ///
    /// [`find_kmers`]: FindKmers::find_kmers
    #[inline]
    pub fn find_all_kmers<const MAX_LEN: usize, T>(self, kmers: &'a T) -> impl Iterator<Item = Range<usize>> + 'a
    where
        T: KmerCollectionContains<MAX_LEN>, {
        kmers.find_all_in_seq(self.slice).map(move |r| self.adjust_to_context(&r))
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in the sequence and those in a collection, in reverse order.
    ///
    /// Consider using [`find_kmers_rev`] if only a single match is needed. A
    /// similar unrestricted version of this method exists in [`FindKmers`].
    ///
    /// [`find_kmers_rev`]: FindKmers::find_kmers_rev
    #[inline]
    pub fn find_all_kmers_rev<const MAX_LEN: usize, T>(self, kmers: &'a T) -> impl Iterator<Item = Range<usize>> + 'a
    where
        T: KmerCollectionContains<MAX_LEN>, {
        kmers.find_all_in_seq_rev(self.slice).map(move |r| self.adjust_to_context(&r))
    }

    /// Takes a closure and returns the index of the first byte found within
    /// the sequence that satisfies the closure. Wraps [`Iterator::position`].
    ///
    /// [`Iterator::position`]: std::iter::Iterator::position
    #[inline]
    pub fn position<P>(&self, mut predicate: P) -> Option<usize>
    where
        P: FnMut(u8) -> bool, {
        self.slice
            .iter()
            .position(|&item| predicate(item))
            .map(|index| self.adjust_to_context(&index))
    }

    /// Takes a closure and returns the index of the last byte found within the
    /// sequence that satisfies the closure. Wraps [`Iterator::rposition`].
    ///
    /// [`Iterator::rposition`]: std::iter::Iterator::rposition
    #[inline]
    pub fn rposition<P>(&self, mut predicate: P) -> Option<usize>
    where
        P: FnMut(u8) -> bool, {
        self.slice
            .iter()
            .rposition(|&item| predicate(item))
            .map(|index| self.adjust_to_context(&index))
    }
}

/// Trait for performing restricted string searches on byte substrings. In
/// particular, it provides methods for generating a [`RangeSearch`] struct from
/// a byte string. See [Restricting the search
/// range](crate::search#restricting-the-search-range) for more details.
pub trait ToRangeSearch: AsRef<[u8]> {
    /// Restrict the search to be in `range`.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// All string search methods called on the resulting struct return indices
    /// with respect to the original sequence, not the subsequence.
    ///
    /// </div>
    #[inline]
    fn search_in<R: SliceRange>(&self, range: R) -> RangeSearch {
        RangeSearch::new(self, range)
    }

    /// Restrict the search to be in only the first `n` bytes. If the sequence
    /// is less than `n` bytes long, then the full sequence is searched.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// All string search methods called on the resulting struct return indices
    /// with respect to the original sequence, not the subsequence.
    ///
    /// </div>
    #[inline]
    fn search_in_first(&self, n: usize) -> RangeSearch {
        let seq = self.as_ref();
        RangeSearch::new(seq, ..n.min(seq.len()))
    }

    /// Restrict the search to be in only the last `n` bytes. If the sequence
    /// is less than `n` bytes long, then the full sequence is searched.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// All string search methods called on the resulting struct return indices
    /// with respect to the original sequence, not the subsequence.
    ///
    /// </div>
    #[inline]
    fn search_in_last(&self, n: usize) -> RangeSearch {
        let seq = self.as_ref();
        RangeSearch::new(seq, seq.len().saturating_sub(n)..)
    }
}

impl<T: AsRef<[u8]>> ToRangeSearch for T {}
