use std::{
    cmp::Ordering,
    ops::{Add, AddAssign},
};

use crate::alignment::phmm::indexing::AlnIndexable;

/// A trait representing different ways to index into sequence data associated
/// with a dynamic programming alignment algorithm.
pub trait AlnIndex: Copy {
    /// Returns the index as a [`DpIndex`].
    ///
    /// It is not checked whether the index is past the end of `seq`.
    #[must_use]
    fn to_dp_index<Q>(self, seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized;

    /// Returns the index as a [`SeqIndex`].
    ///
    /// If the index corresponds to [`Begin`], then `None` is returned since
    /// this does not correspond to a position in the sequence. It is not
    /// checked whether the index is past the end of `seq`.
    #[inline]
    #[must_use]
    fn to_seq_index<Q>(&self, seq: &Q) -> Option<SeqIndex>
    where
        Q: AlnIndexable + ?Sized, {
        self.to_dp_index(seq).0.checked_sub(1).map(SeqIndex)
    }

    /// Gets the index before the current one, as a [`DpIndex`].
    ///
    /// If the index is equivalent to [`Begin`], then `None` is returned.
    #[inline]
    #[must_use]
    fn prev_index<Q>(self, seq: &Q) -> Option<DpIndex>
    where
        Q: AlnIndexable + ?Sized, {
        self.to_dp_index(seq).0.checked_sub(1).map(DpIndex)
    }

    /// Gets the index after the current one, as a [`DpIndex`].
    ///
    /// No check is performed for whether this index is in bounds for the given
    /// `seq`.
    #[inline]
    #[must_use]
    fn next_index<Q>(self, seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        DpIndex(self.to_dp_index(seq).0 + 1)
    }

    /// Gets the minimum of two [`AlnIndex`] structs as a [`DpIndex`] (the
    /// leftmost in `seq`).
    ///
    /// No bounds checking is performed.
    #[inline]
    #[must_use]
    fn min_index<Q>(self, other: impl AlnIndex, seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        self.to_dp_index(seq).min(other.to_dp_index(seq))
    }

    /// Gets the maximum of two [`AlnIndex`] structs as a [`DpIndex`] (the
    /// rightmost in `seq`).
    ///
    /// No bounds checking is performed.
    #[inline]
    #[must_use]
    fn max_index<Q>(self, other: impl AlnIndex, seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        self.to_dp_index(seq).max(other.to_dp_index(seq))
    }

    /// Tests two indices for equality by converting them both to [`DpIndex`].
    ///
    /// No bounds checking is performed.
    #[inline]
    #[must_use]
    fn eq_index<Q>(self, other: impl AlnIndex, seq: &Q) -> bool
    where
        Q: AlnIndexable + ?Sized, {
        self.to_dp_index(seq) == other.to_dp_index(seq)
    }
}

/// An [`AlnIndex`] with respect to a dynamic programming table.
///
/// For sequences, 0 represents no residues aligned and 1 represents the first
/// residue. For pHMMs, 0 represents the BEGIN state and 1 represents the first
/// match state with emissions.
#[repr(transparent)]
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct DpIndex(pub usize);

/// An [`AlnIndex`] with respect to the sequence coordinates.
///
/// For sequences, 0 represents the first residue. For pHMMs, 0 represents the
/// first match state with emissions (the first reference coordinate position).
#[repr(transparent)]
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct SeqIndex(pub usize);

/// An [`AlnIndex`] representing no residues aligned or the BEGIN state of a
/// pHMM, equivalent to `DpIndex(0)`.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct Begin;

/// An [`AlnIndex`] representing the first residue or the first match state with
/// emissions, equivalent to `DpIndex(1)` or `SeqIndex(0)`.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct FirstResidue;

/// An [`AlnIndex`] representing the last residue or the last match state with
/// emissions, whose value depends on the length of the sequence/pHMM.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct LastResidue;

/// An [`AlnIndex`] after [`LastResidue`]. This represents the END state of a
/// pHMM, and also occurs as an exclusive end bound on ranges.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct End;

impl AlnIndex for DpIndex {
    fn to_dp_index<Q>(self, _seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        self
    }
}

impl AlnIndex for SeqIndex {
    fn to_dp_index<Q>(self, _seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        DpIndex(self.0 + 1)
    }
}

impl AlnIndex for Begin {
    #[inline]
    fn to_dp_index<Q>(self, _seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        DpIndex(0)
    }
}

impl AlnIndex for FirstResidue {
    #[inline]
    fn to_dp_index<Q>(self, _seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        DpIndex(1)
    }
}

impl AlnIndex for LastResidue {
    #[inline]
    fn to_dp_index<Q>(self, seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        DpIndex(seq.seq_len())
    }
}

impl AlnIndex for End {
    #[inline]
    fn to_dp_index<Q>(self, seq: &Q) -> DpIndex
    where
        Q: AlnIndexable + ?Sized, {
        DpIndex(seq.seq_len() + 1)
    }
}

impl Begin {
    #[must_use]
    pub fn to_dp_index(self) -> DpIndex {
        DpIndex(0)
    }
}

impl FirstResidue {
    #[must_use]
    pub fn to_dp_index(self) -> DpIndex {
        DpIndex(1)
    }

    #[must_use]
    pub fn to_seq_index(self) -> SeqIndex {
        SeqIndex(0)
    }
}

impl SeqIndex {
    #[must_use]
    pub fn to_dp_index(self) -> DpIndex {
        DpIndex(self.0 + 1)
    }
}

impl DpIndex {
    #[must_use]
    pub fn to_seq_index(self) -> Option<SeqIndex> {
        self.0.checked_sub(1).map(SeqIndex)
    }
}

impl End {
    /// An inherent method override for [`AlnIndex::to_seq_index`] that is
    /// infallible.
    #[must_use]
    pub fn to_seq_index<Q>(self, seq: &Q) -> SeqIndex
    where
        Q: AlnIndexable + ?Sized, {
        SeqIndex(seq.seq_len())
    }
}

impl PartialEq<Begin> for SeqIndex {
    fn eq(&self, _other: &Begin) -> bool {
        // Begin has DpIndex 0, and SeqIndex has DpIndex >= 1
        false
    }
}

impl PartialEq<SeqIndex> for Begin {
    fn eq(&self, _other: &SeqIndex) -> bool {
        // Begin is DpIndex(0), and SeqIndex has DpIndex >= 1
        false
    }
}

impl PartialOrd<Begin> for SeqIndex {
    fn partial_cmp(&self, _other: &Begin) -> Option<Ordering> {
        // Begin is DpIndex(0), and SeqIndex has DpIndex >= 1
        Some(Ordering::Greater)
    }
}

impl PartialOrd<SeqIndex> for Begin {
    fn partial_cmp(&self, _other: &SeqIndex) -> Option<Ordering> {
        // Begin is DpIndex(0), and SeqIndex has DpIndex >= 1
        Some(Ordering::Less)
    }
}

impl PartialEq<FirstResidue> for SeqIndex {
    fn eq(&self, other: &FirstResidue) -> bool {
        *self == (*other).to_seq_index()
    }
}

impl PartialEq<SeqIndex> for FirstResidue {
    fn eq(&self, other: &SeqIndex) -> bool {
        (*self).to_seq_index() == *other
    }
}

impl PartialOrd<FirstResidue> for SeqIndex {
    fn partial_cmp(&self, other: &FirstResidue) -> Option<Ordering> {
        self.partial_cmp(&(*other).to_seq_index())
    }
}

impl PartialOrd<SeqIndex> for FirstResidue {
    fn partial_cmp(&self, other: &SeqIndex) -> Option<Ordering> {
        (*self).to_seq_index().partial_cmp(other)
    }
}

impl PartialEq<Begin> for DpIndex {
    fn eq(&self, other: &Begin) -> bool {
        *self == other.to_dp_index()
    }
}

impl PartialEq<DpIndex> for Begin {
    fn eq(&self, other: &DpIndex) -> bool {
        self.to_dp_index() == *other
    }
}

impl PartialOrd<Begin> for DpIndex {
    fn partial_cmp(&self, other: &Begin) -> Option<Ordering> {
        self.partial_cmp(&other.to_dp_index())
    }
}

impl PartialOrd<DpIndex> for Begin {
    fn partial_cmp(&self, other: &DpIndex) -> Option<Ordering> {
        self.to_dp_index().partial_cmp(other)
    }
}

impl PartialEq<FirstResidue> for DpIndex {
    fn eq(&self, other: &FirstResidue) -> bool {
        *self == other.to_dp_index()
    }
}

impl PartialEq<DpIndex> for FirstResidue {
    fn eq(&self, other: &DpIndex) -> bool {
        self.to_dp_index() == *other
    }
}

impl PartialOrd<FirstResidue> for DpIndex {
    fn partial_cmp(&self, other: &FirstResidue) -> Option<Ordering> {
        self.partial_cmp(&other.to_dp_index())
    }
}

impl PartialOrd<DpIndex> for FirstResidue {
    fn partial_cmp(&self, other: &DpIndex) -> Option<Ordering> {
        self.to_dp_index().partial_cmp(other)
    }
}

impl PartialEq<DpIndex> for SeqIndex {
    fn eq(&self, other: &DpIndex) -> bool {
        self.to_dp_index() == *other
    }
}

impl PartialEq<SeqIndex> for DpIndex {
    fn eq(&self, other: &SeqIndex) -> bool {
        *self == other.to_dp_index()
    }
}

impl PartialOrd<DpIndex> for SeqIndex {
    fn partial_cmp(&self, other: &DpIndex) -> Option<Ordering> {
        self.to_dp_index().partial_cmp(other)
    }
}

impl PartialOrd<SeqIndex> for DpIndex {
    fn partial_cmp(&self, other: &SeqIndex) -> Option<Ordering> {
        self.partial_cmp(&other.to_dp_index())
    }
}

impl Add<usize> for DpIndex {
    type Output = DpIndex;

    fn add(mut self, rhs: usize) -> Self::Output {
        self.0 += rhs;
        self
    }
}

impl Add<usize> for SeqIndex {
    type Output = SeqIndex;

    fn add(mut self, rhs: usize) -> Self::Output {
        self.0 += rhs;
        self
    }
}

impl AddAssign<usize> for DpIndex {
    fn add_assign(&mut self, rhs: usize) {
        self.0 += rhs;
    }
}

impl AddAssign<usize> for SeqIndex {
    fn add_assign(&mut self, rhs: usize) {
        self.0 += rhs;
    }
}
