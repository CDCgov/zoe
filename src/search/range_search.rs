use crate::data::view_traits::{IndexAdjustable, SliceRange};
use std::ops::{Bound, Index};

/// A subsequence along with its starting index. This is used for restricting
/// the search range when performing a string search. See [Restricting the
/// search range](crate::search#restricting-the-search-range) for more details.
///
/// <div class="warning">
/// All string search methods called on this struct return indices with respect
/// to the original sequence, not the subsequence.
/// </div>
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
}

/// Trait for performing restricted string searches on byte substrings. In
/// particular, it provides methods for generating a [`RangeSearch`] struct from
/// a byte string. See [Restricting the search
/// range](crate::search#restricting-the-search-range) for more details.
pub trait ToRangeSearch: AsRef<[u8]> {
    /// Restrict the search to be in `range`.
    ///
    /// <div class="warning">
    /// All string search methods called on the resulting struct return indices
    /// with respect to the original sequence, not the subsequence.
    /// </div>
    #[inline]
    fn search_in<R: SliceRange>(&self, range: R) -> RangeSearch {
        RangeSearch::new(self, range)
    }

    /// Restrict the search to be in only the first `n` bytes. If the sequence
    /// is less than `n` bytes long, then the full sequence is searched.
    ///
    /// <div class="warning">
    /// All string search methods called on the resulting struct return indices
    /// with respect to the original sequence, not the subsequence.
    /// </div>
    #[inline]
    fn search_in_first(&self, n: usize) -> RangeSearch {
        let seq = self.as_ref();
        RangeSearch::new(seq, ..n.min(seq.len()))
    }

    /// Restrict the search to be in only the last `n` bytes. If the sequence
    /// is less than `n` bytes long, then the full sequence is searched.
    ///
    /// <div class="warning">
    /// All string search methods called on the resulting struct return indices
    /// with respect to the original sequence, not the subsequence.
    /// </div>
    #[inline]
    fn search_in_last(&self, n: usize) -> RangeSearch {
        let seq = self.as_ref();
        RangeSearch::new(seq, seq.len().saturating_sub(n)..)
    }
}

impl<T: AsRef<[u8]>> ToRangeSearch for T {}
