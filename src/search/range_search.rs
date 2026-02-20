use crate::{
    data::views::{IndexAdjustable, SliceRange},
    kmer::{FindKmers, FindKmersInSeq},
    prelude::{NucleotidesView, Translate},
    private::Sealed,
};
use std::{
    marker::PhantomData,
    ops::{Bound, Range},
};

/// A subsequence along with its starting index. This is used for restricting
/// the search range when performing a search. See [Restricting the search
/// range](crate::search#restricting-the-search-range) for more details.
///
/// It is compatible with:
///
/// - Performing string search via the [`ByteSubstring`] trait
/// - Searching for k-mers using methods similar to those in [`FindKmers`]
/// - Searching for amino acids in a nucleotides sequence similar to methods in
///   [`Translate`]. This requires the original type to implement [`Translate`]
///   also
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
///
/// ## Parameters
///
/// `Q` tracks the type of the original sequence type. This controls whether
/// [`find_next_aa`] and [`find_next_aa_in_frame`] are available.
///
/// [`find_next_aa`]: RangeSearch::find_next_aa
/// [`find_next_aa_in_frame`]: RangeSearch::find_next_aa_in_frame
pub struct RangeSearch<'a, Q: ?Sized> {
    pub(crate) slice:       &'a [u8],
    pub(crate) starting_at: usize,
    /// A marker for the type of input data
    phantom:                PhantomData<Q>,
}

impl<'a, Q: ?Sized> RangeSearch<'a, Q> {
    /// Creates a new [`RangeSearch`] from a byte sequence and a range.
    fn new<R>(sequence: &'a [u8], range: R) -> Self
    where
        R: SliceRange, {
        let slice = &sequence[range.clone()];
        let starting_at = match range.start_bound() {
            Bound::Included(&start) => start,
            Bound::Unbounded => 0,
            Bound::Excluded(&start) => start.saturating_add(1),
        };
        Self {
            slice,
            starting_at,
            phantom: PhantomData,
        }
    }
}

impl<Q> RangeSearch<'_, Q>
where
    Q: Translate,
{
    /// Slides base by base over the sequence until the next translated `aa` is
    /// found and returns that index (or `None` otherwise).
    #[inline]
    #[must_use]
    pub fn find_next_aa(self, aa: u8) -> Option<usize> {
        NucleotidesView::from(self.slice)
            .find_next_aa(aa)
            .map(|index| self.adjust_to_context(&index))
    }

    /// Slides codon by codon over the sequence (reading frame starting from 0
    /// in the full sequence) until the next translated `aa` is found and
    /// returns that index (or `None` otherwise).
    ///
    /// This is different than calling [`Translate::find_next_aa_in_frame`] on a
    /// slice, since it uses the reading frame of the full sequence, irrelevant
    /// of the starting position of the search.
    #[inline]
    #[must_use]
    pub fn find_next_aa_in_frame(mut self, aa: u8) -> Option<usize> {
        let starting_at_in_frame = self.starting_at.next_multiple_of(3);
        let offset = starting_at_in_frame - self.starting_at;
        self.move_forward_by(offset);
        NucleotidesView::from(self.slice)
            .find_next_aa_in_frame(aa)
            .map(|index| self.adjust_to_context(&index))
    }
}

impl<'a, Q: ?Sized> RangeSearch<'a, Q> {
    /// Moves the starting position of the [`RangeSearch`] forward by `amt`.
    #[inline]
    pub(crate) fn move_forward_by(&mut self, amt: usize) {
        self.starting_at += amt;
        self.slice = &self.slice[amt..];
    }

    /// Given an index/range in the frame of reference of `slice`, adjust it to
    /// be in the original frame of reference.
    #[inline]
    pub(crate) fn adjust_to_context<I: IndexAdjustable>(&self, index: &I) -> I {
        index.add(self.starting_at)
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

    // FindKmers methods. These need to be distinct from the FindKmers trait
    // because the methods need to consume the RangeSearch struct, rather than
    // take them by reference. Otherwise, a bunch of useless let bindings are
    // required.

    /// Returns the indices of the leftmost occurrence of any of the k-mers in
    /// the passed collection.
    ///
    /// If no occurrence is found, then `None` is returned. To find all matches
    /// rather than just the first, use [`find_all_kmers`]. A similar
    /// unrestricted version of this method exists in [`FindKmers`].
    ///
    /// [`find_all_kmers`]: FindKmers::find_all_kmers
    /// [`KmerEncoder`]: crate::kmer::KmerEncoder
    #[inline]
    pub fn find_kmers<const MAX_LEN: usize, T>(self, kmers: &T) -> Option<Range<usize>>
    where
        T: FindKmersInSeq<MAX_LEN>, {
        self.slice.find_kmers(kmers).map(|r| self.adjust_to_context(&r))
    }

    /// Returns the indices of the rightmost occurrence of any of the k-mers in
    /// the passed collection.
    ///
    /// If no occurrence is found, then `None` is returned. To find all matches
    /// rather than just the last, use [`find_all_kmers`]. A similar
    /// unrestricted version of this method exists in [`FindKmers`].
    ///
    /// [`find_all_kmers`]: FindKmers::find_all_kmers
    /// [`KmerEncoder`]: crate::kmer::KmerEncoder
    #[inline]
    pub fn find_kmers_rev<const MAX_LEN: usize, T>(self, kmers: &T) -> Option<Range<usize>>
    where
        T: FindKmersInSeq<MAX_LEN>, {
        self.slice.find_kmers_rev(kmers).map(|r| self.adjust_to_context(&r))
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in the sequence and those in the passed collection.
    ///
    /// Consider using [`find_kmers`] if only a single match is needed. A
    /// similar unrestricted version of this method exists in [`FindKmers`].
    ///
    /// [`find_kmers`]: FindKmers::find_kmers
    #[inline]
    pub fn find_all_kmers<const MAX_LEN: usize, T>(self, kmers: &'a T) -> impl Iterator<Item = Range<usize>> + 'a
    where
        Q: 'a,
        T: FindKmersInSeq<MAX_LEN>, {
        kmers
            .find_all_in_seq(self.slice.as_ref())
            .map(move |r| self.adjust_to_context(&r))
    }

    /// Returns an iterator over the indices of all matches between the k-mers
    /// in the sequence and those in the passed collection, in reverse order.
    ///
    /// Consider using [`find_kmers_rev`] if only a single match is needed. A
    /// similar unrestricted version of this method exists in [`FindKmers`].
    ///
    /// [`find_kmers_rev`]: FindKmers::find_kmers_rev
    #[inline]
    pub fn find_all_kmers_rev<const MAX_LEN: usize, T>(self, kmers: &'a T) -> impl Iterator<Item = Range<usize>> + 'a
    where
        Q: 'a,
        T: FindKmersInSeq<MAX_LEN>, {
        kmers.find_all_in_seq_rev(self.slice).map(move |r| self.adjust_to_context(&r))
    }
}

/// Trait for performing restricted string searches on byte substrings. In
/// particular, it provides methods for generating a [`RangeSearch`] struct from
/// a byte string. See [Restricting the search
/// range](crate::search#restricting-the-search-range) for more details.
pub trait ToRangeSearch: AsRef<[u8]> + Sealed {
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
    fn search_in<R: SliceRange>(&self, range: R) -> RangeSearch<'_, Self> {
        RangeSearch::new(self.as_ref(), range)
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
    fn search_in_first(&self, n: usize) -> RangeSearch<'_, Self> {
        Self::search_in(self, ..n.min(self.as_ref().len()))
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
    fn search_in_last(&self, n: usize) -> RangeSearch<'_, Self> {
        Self::search_in(self, self.as_ref().len().saturating_sub(n)..)
    }
}

impl<T: AsRef<[u8]> + Sealed> ToRangeSearch for T {}

#[cfg(test)]
mod test {
    use super::*;
    use crate::search::ByteSubstring;

    #[test]
    fn search_nucleotides() {
        let seq = NucleotidesView::from(b"CCCATGCCCCATGCCATG");
        let idx = seq.search_in(8..).find_next_aa_in_frame(b'M');
        assert_eq!(idx, Some(15));
    }

    #[test]
    fn search_regular() {
        let seq = b"GGGAAGCATCACGTATCGA";
        let contains = seq.search_in(2..).contains_substring(b"GGG");
        assert!(!contains);
    }
}
