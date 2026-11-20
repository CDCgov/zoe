//! Other alignment output types besides [`Alignment`].

use crate::alignment::Alignment;
use std::ops::Range;

/// The output of an alignment algorithm, allowing unmapped and overflowed
/// alignments to be represented.
///
/// Often this is some [`Alignment`], but if no portion of the query mapped to
/// the reference, [`Unmapped`] is used. [`Overflowed`] is used when the score
/// exceeded numeric type range.
///
/// [`Some`]: MaybeAligned::Some
/// [`Overflowed`]: MaybeAligned::Overflowed
/// [`Unmapped`]: MaybeAligned::Unmapped
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum MaybeAligned<T> {
    /// Contains a successful alignment result
    Some(T),
    /// Indicates the alignment score overflowed the numeric type
    Overflowed,
    /// Indicates the sequence could not be mapped/aligned
    Unmapped,
}

impl<T> MaybeAligned<T> {
    /// Unwraps the alignment, consuming the [`MaybeAligned<T>`] and returning
    /// the contained [`Alignment<T>`].
    ///
    /// # Panics
    ///
    /// Panics if the value is [`MaybeAligned::Overflowed`] or
    /// [`MaybeAligned::Unmapped`].
    #[inline]
    #[must_use]
    pub fn unwrap(self) -> T {
        match self {
            MaybeAligned::Some(aln) => aln,
            MaybeAligned::Overflowed => panic!("Alignment score overflowed!"),
            MaybeAligned::Unmapped => panic!("Sequence could not be mapped!"),
        }
    }

    /// Gets the contained [`Alignment`] as an [`Option`].
    #[inline]
    #[must_use]
    pub fn get(self) -> Option<T> {
        match self {
            MaybeAligned::Some(aln) => Some(aln),
            _ => None,
        }
    }

    /// Converts from `&MaybeAligned<T>` to `MaybeAligned<&T>`.
    #[inline]
    #[must_use]
    pub const fn as_ref(&self) -> MaybeAligned<&T> {
        match *self {
            MaybeAligned::Some(ref x) => MaybeAligned::Some(x),
            MaybeAligned::Overflowed => MaybeAligned::Overflowed,
            MaybeAligned::Unmapped => MaybeAligned::Unmapped,
        }
    }

    /// Converts from `&mut MaybeAligned<T>` to `MaybeAligned<&mut T>`.
    #[inline]
    #[must_use]
    pub const fn as_mut(&mut self) -> MaybeAligned<&mut T> {
        match *self {
            MaybeAligned::Some(ref mut x) => MaybeAligned::Some(x),
            MaybeAligned::Overflowed => MaybeAligned::Overflowed,
            MaybeAligned::Unmapped => MaybeAligned::Unmapped,
        }
    }

    /// Returns the alignment if it did not overflow, otherwise calls `f` and
    /// returns the result.
    #[inline]
    #[must_use]
    pub fn or_else_overflowed(self, f: impl FnOnce() -> Self) -> Self {
        if let MaybeAligned::Overflowed = self { f() } else { self }
    }

    /// Maps a [`MaybeAligned<T>`] to [`MaybeAligned<U>`] by applying a function
    /// to the contained value. Leaves [`Overflowed`] and [`Unmapped`] variants
    /// unchanged.
    ///
    /// [`Overflowed`]: MaybeAligned::Overflowed
    /// [`Unmapped`]: MaybeAligned::Unmapped
    #[inline]
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> MaybeAligned<U> {
        match self {
            MaybeAligned::Some(value) => MaybeAligned::Some(f(value)),
            MaybeAligned::Overflowed => MaybeAligned::Overflowed,
            MaybeAligned::Unmapped => MaybeAligned::Unmapped,
        }
    }

    /// Returns [`Unmapped`] if the [`MaybeAligned`] is [`Unmapped`],
    /// [`Overflowed`] if it is [`Overflowed`], otherwise calls `f` with the
    /// wrapped value and returns the result.
    ///
    /// [`Overflowed`]: MaybeAligned::Overflowed
    /// [`Unmapped`]: MaybeAligned::Unmapped
    #[inline]
    #[must_use]
    pub fn and_then<U, F>(self, f: F) -> MaybeAligned<U>
    where
        F: FnOnce(T) -> MaybeAligned<U>, {
        match self {
            MaybeAligned::Some(x) => f(x),
            MaybeAligned::Overflowed => MaybeAligned::Overflowed,
            MaybeAligned::Unmapped => MaybeAligned::Unmapped,
        }
    }
}

impl<T> MaybeAligned<Alignment<T>>
where
    T: Copy,
{
    /// Gets the alignment for when the query and reference are both reversed.
    ///
    /// See [`Reversibility`](Alignment#reversibility) for more details.
    #[inline]
    #[must_use]
    pub fn to_reverse(&self) -> Self {
        self.as_ref().map(Alignment::to_reverse)
    }

    /// Converts an alignment in-place so that it represents the alignment
    /// between the reversed query and reference.
    ///
    /// See [`Reversibility`](Alignment#reversibility) for more details.
    #[inline]
    pub fn make_reverse(&mut self) {
        self.as_mut().map(Alignment::make_reverse);
    }
}

/// An alignment output for helper functions returning the score, an index in
/// the reference, and an index in the query.
///
/// For user-facing algorithms, [`ScoreStarts`] or [`ScoreEnds`] should be used
/// to make the meanings of the fields explicit. This is used only on internal
/// helper functions that can return either one.
pub(crate) struct ScoreIndices<T> {
    pub score:     T,
    pub ref_idx:   usize,
    pub query_idx: usize,
}

impl<T> ScoreIndices<T> {
    /// Converts the [`ScoreIndices`] to the more explicit [`ScoreStarts`].
    ///
    /// ## Validity
    ///
    /// This should only be called when the indices in `self` represent the
    /// start of the alignment in each sequence.
    #[inline]
    #[must_use]
    pub fn into_score_starts(self) -> ScoreStarts<T> {
        ScoreStarts {
            score:       self.score,
            ref_start:   self.ref_idx,
            query_start: self.query_idx,
        }
    }

    /// Converts the [`ScoreIndices`] to the more explicit [`ScoreEnds`].
    ///
    /// ## Validity
    ///
    /// This should only be called when the indices in `self` represent the end
    /// of the alignment in each sequence.
    #[inline]
    #[must_use]
    pub fn into_score_ends(self) -> ScoreEnds<T> {
        ScoreEnds {
            score:     self.score,
            ref_end:   self.ref_idx,
            query_end: self.query_idx,
        }
    }
}

/// An alignment output for algorithms returning the score, the starting index
/// in the reference, and the starting index in the query.
pub struct ScoreStarts<T> {
    /// The score of the alignment.
    pub score:       T,
    /// The starting index of the alignment within the reference.
    pub ref_start:   usize,
    /// The starting index of the alignment within the query.
    pub query_start: usize,
}

/// An alignment output for algorithms returning the score, the ending index in
/// the reference, and the ending index in the query.
pub struct ScoreEnds<T> {
    /// The score of the alignment.
    pub score:     T,
    /// The ending index of the alignment within the reference.
    pub ref_end:   usize,
    /// The ending index of the alignment within the query.
    pub query_end: usize,
}

/// An alignment output for algorithms returning the score, the range of indices
/// aligned in the reference, and the range of indices aligned in the query.
pub struct ScoreAndRanges<T> {
    /// The score of the alignment.
    pub score:       T,
    /// The range of indices in the alignment for the reference sequence.
    pub ref_range:   Range<usize>,
    /// The range of indices in the alignment for the query sequence.
    pub query_range: Range<usize>,
}
