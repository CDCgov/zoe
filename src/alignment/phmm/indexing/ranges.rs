use crate::alignment::phmm::indexing::{AlnIndex, AlnIndexable, Begin, DpIndex, End, FirstResidue, SeqIndex};
use std::ops::{Bound, Range, RangeBounds, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive};

/// A trait similar to [`RangeBounds`] but for ranges of [`AlnIndex`] values.
pub trait AlnIndexRange {
    /// The type of the starting index.
    type Start: AlnIndex;
    /// The type of the ending index.
    type End: AlnIndex;

    /// Returns the starting index in the range as a [`Bound`].
    fn start_bound(&self) -> Bound<Self::Start>;

    /// Returns the ending index in the range as a [`Bound`].
    fn end_bound(&self) -> Bound<Self::End>;

    /// Returns the range as dynamic programming indices.
    #[inline]
    #[must_use]
    fn to_dp_range<Q>(&self, seq: &Q) -> Range<DpIndex>
    where
        Q: AlnIndexable + ?Sized, {
        let start = match self.start_bound() {
            Bound::Included(start) => start.to_dp_index(seq),
            // +1 due to converting Excluded to Included
            Bound::Excluded(start) => start.to_dp_index(seq) + 1,
            Bound::Unbounded => Begin.to_dp_index(),
        };
        let end = match self.end_bound() {
            // +1 due to converting Included to Excluded
            Bound::Included(end) => end.to_dp_index(seq) + 1,
            Bound::Excluded(end) => end.to_dp_index(seq),
            Bound::Unbounded => End.to_dp_index(seq),
        };

        start..end
    }

    /// Returns the range as sequence indices, saturating at the ends. If either
    /// value is [`Begin`], this will replace it with [`FirstResidue`]. The
    /// returned range will be shortened so as not to exceed an exclusive end of
    /// [`End`].
    #[inline]
    #[must_use]
    fn saturating_to_seq_range(&self, seq: &impl AlnIndexable) -> Range<SeqIndex> {
        let Range { start, end } = self.to_dp_range(seq);
        let start = start.to_seq_index().unwrap_or(FirstResidue.to_seq_index());
        let end = end
            .to_seq_index()
            .unwrap_or(FirstResidue.to_seq_index())
            .min(End.to_seq_index(seq));
        let start = start.min(end);
        start..end
    }

    /// Returns an iterator over the [`DpIndex`] values in the range.
    fn iter_dp_index(&self, seq: &impl AlnIndexable) -> impl Iterator<Item = DpIndex> {
        self.to_dp_range(seq).into_inner().into_iter().map(DpIndex)
    }
}

impl<I: AlnIndex, J: AlnIndex> AlnIndexRange for (Bound<I>, Bound<J>) {
    type Start = I;
    type End = J;

    #[inline]
    fn start_bound(&self) -> Bound<Self::Start> {
        self.0
    }

    #[inline]
    fn end_bound(&self) -> Bound<Self::End> {
        self.1
    }
}

impl<I: AlnIndex> AlnIndexRange for Range<I> {
    type Start = I;
    type End = I;

    fn start_bound(&self) -> Bound<Self::Start> {
        <Self as RangeBounds<I>>::start_bound(self).map(|x| *x)
    }

    fn end_bound(&self) -> Bound<Self::End> {
        <Self as RangeBounds<I>>::end_bound(self).map(|x| *x)
    }
}

impl<I: AlnIndex> AlnIndexRange for RangeInclusive<I> {
    type Start = I;
    type End = I;

    fn start_bound(&self) -> Bound<Self::Start> {
        <Self as RangeBounds<I>>::start_bound(self).map(|x| *x)
    }

    fn end_bound(&self) -> Bound<Self::End> {
        <Self as RangeBounds<I>>::end_bound(self).map(|x| *x)
    }
}

impl<I: AlnIndex> AlnIndexRange for RangeFrom<I> {
    type Start = I;
    type End = I;

    fn start_bound(&self) -> Bound<Self::Start> {
        <Self as RangeBounds<I>>::start_bound(self).map(|x| *x)
    }

    fn end_bound(&self) -> Bound<Self::End> {
        <Self as RangeBounds<I>>::end_bound(self).map(|x| *x)
    }
}

impl<I: AlnIndex> AlnIndexRange for RangeTo<I> {
    type Start = I;
    type End = I;

    fn start_bound(&self) -> Bound<Self::Start> {
        <Self as RangeBounds<I>>::start_bound(self).map(|x| *x)
    }

    fn end_bound(&self) -> Bound<Self::End> {
        <Self as RangeBounds<I>>::end_bound(self).map(|x| *x)
    }
}

impl<I: AlnIndex> AlnIndexRange for RangeToInclusive<I> {
    type Start = I;
    type End = I;

    fn start_bound(&self) -> Bound<Self::Start> {
        <Self as RangeBounds<I>>::start_bound(self).map(|x| *x)
    }

    fn end_bound(&self) -> Bound<Self::End> {
        <Self as RangeBounds<I>>::end_bound(self).map(|x| *x)
    }
}

/// An extension trait for `Range<DpIndex>` and `Range<SeqIndex>` to convert
/// them into `Range<usize>`.
pub trait IndexRangeInner {
    /// Extracts the contained `usize` indices within the range.
    fn into_inner(self) -> Range<usize>;
}

impl IndexRangeInner for Range<DpIndex> {
    #[inline]
    fn into_inner(self) -> Range<usize> {
        self.start.0..self.end.0
    }
}

impl IndexRangeInner for Range<SeqIndex> {
    #[inline]
    fn into_inner(self) -> Range<usize> {
        self.start.0..self.end.0
    }
}
