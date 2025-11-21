use std::{
    ops::{Range, RangeBounds, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive},
    slice::SliceIndex,
};

/// A trait alias for range types (with `usize` bounds, used to index into
/// `&[u8]`).
pub trait SliceRange: SliceIndex<[u8], Output = [u8]> + Clone + RangeBounds<usize> {}
impl<T: SliceIndex<[u8], Output = [u8]> + Clone + RangeBounds<usize>> SliceRange for T {}

/// Helper trait for adding a constant value to an index or range.
pub(crate) trait IndexAdjustable {
    /// Adds a constant value to the range or index.
    #[must_use]
    fn add(&self, n: usize) -> Self;

    /// Subtracts a constant value to the range or index.
    #[must_use]
    #[allow(dead_code)]
    fn sub(&self, n: usize) -> Self;

    /// Subtracts a constant value to the range or index, saturating at bounds.
    #[must_use]
    #[allow(dead_code)]
    fn saturating_sub(&self, n: usize) -> Self;
}

impl IndexAdjustable for usize {
    #[inline]
    fn add(&self, n: usize) -> Self {
        self + n
    }

    #[inline]
    fn sub(&self, n: usize) -> Self {
        self - n
    }

    #[inline]
    fn saturating_sub(&self, n: usize) -> Self {
        usize::saturating_sub(*self, n)
    }
}

impl IndexAdjustable for Range<usize> {
    #[inline]
    fn add(&self, n: usize) -> Self {
        self.start + n..self.end + n
    }

    #[inline]
    fn sub(&self, n: usize) -> Self {
        self.start - n..self.end - n
    }

    #[inline]
    fn saturating_sub(&self, n: usize) -> Self {
        self.start.saturating_sub(n)..self.end.saturating_sub(n)
    }
}

impl IndexAdjustable for RangeInclusive<usize> {
    #[inline]
    fn add(&self, n: usize) -> Self {
        self.start() + n..=self.end() + n
    }

    #[inline]
    fn sub(&self, n: usize) -> Self {
        self.start() - n..=self.end() - n
    }

    #[inline]
    fn saturating_sub(&self, n: usize) -> Self {
        self.start().saturating_sub(n)..=self.end().saturating_sub(n)
    }
}

impl IndexAdjustable for RangeFrom<usize> {
    #[inline]
    fn add(&self, n: usize) -> Self {
        self.start + n..
    }

    #[inline]
    fn sub(&self, n: usize) -> Self {
        self.start - n..
    }

    #[inline]
    fn saturating_sub(&self, n: usize) -> Self {
        self.start.saturating_sub(n)..
    }
}

impl IndexAdjustable for RangeTo<usize> {
    #[inline]
    fn add(&self, n: usize) -> Self {
        ..self.end + n
    }

    #[inline]
    fn sub(&self, n: usize) -> Self {
        ..self.end - n
    }

    #[inline]
    fn saturating_sub(&self, n: usize) -> Self {
        ..self.end.saturating_sub(n)
    }
}

impl IndexAdjustable for RangeToInclusive<usize> {
    #[inline]
    fn add(&self, n: usize) -> Self {
        ..=self.end + n
    }

    #[inline]
    fn sub(&self, n: usize) -> Self {
        ..=self.end - n
    }

    #[inline]
    fn saturating_sub(&self, n: usize) -> Self {
        ..=self.end.saturating_sub(n)
    }
}
