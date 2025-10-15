//! ## Functions for aligning sequence data.
//!
//! *Zoe* supports alignment for DNA, protein, or any other sequence data.
//! Currently, two modes of alignment are available:
//!
//! - [Smith-Waterman]: Optimal local alignment in the [`sw`] module.
//! - [Needleman–Wunsch]: Optimal global alignment in the [`nw`] module.
//!
//! The latter currently only has a scalar implementation, so is suitable for
//! testing or applications which are not performance-critical.
//!
//! [Smith-Waterman]:
//!     https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
//! [Needleman–Wunsch]:
//!     https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

pub mod nw;
#[cfg(feature = "dev-phmm")]
pub mod phmm;
pub mod sw;

mod errors;
mod pairwise;
mod profile;
mod profile_set;
mod sneaky_snake;
mod std_traits;
mod types;

pub use errors::*;
pub use pairwise::*;
pub use profile::*;
pub use profile_set::*;
pub use sneaky_snake::*;
pub use types::*;

/// A trait implemented on structs containing data fields related to a query
/// sequence and a reference sequence. It then exposes a method for inverting
/// the struct (swapping the query and reference).
pub(crate) trait HasQueryAndRefData {
    /// Swaps the fields related to the query and reference sequences. This
    /// should include altering other fields as well, such as inverting a CIGAR
    /// string if present.
    fn invert(self) -> Self;
}

impl<T> HasQueryAndRefData for Alignment<T>
where
    T: Copy,
{
    #[inline]
    fn invert(self) -> Self {
        Alignment::invert(&self)
    }
}

impl<T> HasQueryAndRefData for ScoreIndices<T> {
    #[inline]
    fn invert(self) -> Self {
        Self {
            score:     self.score,
            ref_idx:   self.query_idx,
            query_idx: self.ref_idx,
        }
    }
}

impl<T> HasQueryAndRefData for ScoreStarts<T> {
    #[inline]
    fn invert(self) -> Self {
        Self {
            score:       self.score,
            ref_start:   self.query_start,
            query_start: self.ref_start,
        }
    }
}

impl<T> HasQueryAndRefData for ScoreEnds<T> {
    #[inline]
    fn invert(self) -> Self {
        Self {
            score:     self.score,
            ref_end:   self.query_end,
            query_end: self.ref_end,
        }
    }
}

impl<T> HasQueryAndRefData for ScoreAndRanges<T> {
    fn invert(self) -> Self {
        Self {
            score:       self.score,
            ref_range:   self.query_range,
            query_range: self.ref_range,
        }
    }
}

impl<T> HasQueryAndRefData for MaybeAligned<T>
where
    T: HasQueryAndRefData,
{
    #[inline]
    fn invert(self) -> Self {
        self.map(HasQueryAndRefData::invert)
    }
}

impl<T> HasQueryAndRefData for Option<T>
where
    T: HasQueryAndRefData,
{
    #[inline]
    fn invert(self) -> Self {
        self.map(HasQueryAndRefData::invert)
    }
}

impl<T, E> HasQueryAndRefData for Result<T, E>
where
    T: HasQueryAndRefData,
{
    #[inline]
    fn invert(self) -> Self {
        self.map(HasQueryAndRefData::invert)
    }
}

/// An enum wrapping sequence data to indicate whether it represents a query
/// sequence or a reference sequence.
///
/// Wrapping a sequence in [`SeqSrc`] is common for alignment algorithms in
/// *Zoe*, allowing the methods to know which sequence is the query and which is
/// the reference. There are multiple equivalent ways of creating a [`SeqSrc`],
/// including wrapping a sequence directly and using an extension trait. For
/// example:
///
/// ```
/// # use zoe::{
/// #     alignment::{AsSeqSrc, IntoSeqSrc, SeqSrc},
/// #     prelude::Nucleotides,
/// # };
/// let seq = Nucleotides::from(b"AGTCGCTTGCATTTAC");
///
/// // Wrap a reference to the sequence
/// let wrapped1 = SeqSrc::Query(&seq);
/// let wrapped2 = seq.as_query_src();
/// assert_eq!(wrapped1, wrapped2);
///
/// // Wrap the owned sequence
/// let wrapped3 = SeqSrc::Query(seq.clone());
/// let wrapped4 = seq.into_query_src();
/// assert_eq!(wrapped3, wrapped4)
/// ```
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum SeqSrc<T> {
    /// A query sequence of type `T`.
    Query(T),
    /// A reference sequence of type `T`.
    Reference(T),
}

impl<T> SeqSrc<T> {
    /// Given the non-profile sequence wrapped in [`SeqSrc`], constructs an
    /// alignment using closure `f`, automatically inverting the alignment if
    /// the sequence is [`SeqSrc::Query`] (since standalone alignment methods
    /// in *Zoe* assume the profile is build from the query).
    ///
    /// The closure `f` will take as argument whatever sequence was not used to
    /// construct the profile. However, since this function handles inverting
    /// the alignment at the end if needed, it is in most scenarios okay (for
    /// the sake of readability) to name the closure argument `reference`.
    #[inline]
    #[must_use]
    pub(crate) fn make_alignment<F, A>(self, f: F) -> A
    where
        F: FnOnce(T) -> A,
        A: HasQueryAndRefData, {
        let (seq, invert) = match self {
            // Non-profile sequence is the query, so make sure to invert!
            SeqSrc::Query(seq) => (seq, true),
            // Profile sequence is the query, so no further action needed.
            SeqSrc::Reference(seq) => (seq, false),
        };

        let alignment = f(seq);

        if invert { alignment.invert() } else { alignment }
    }

    /// Maps the stored sequence via a closure, preserving whether it is a
    /// [`Query`] or [`Reference`].
    ///
    /// [`Query`]: SeqSrc::Query
    /// [`Reference`]: SeqSrc::Reference
    #[inline]
    #[must_use]
    pub fn map<U, F>(self, f: F) -> SeqSrc<U>
    where
        F: FnOnce(T) -> U, {
        match self {
            SeqSrc::Query(seq) => SeqSrc::Query(f(seq)),
            SeqSrc::Reference(seq) => SeqSrc::Reference(f(seq)),
        }
    }

    /// Returns whether the sequence is a query sequence.
    #[inline]
    #[must_use]
    pub fn is_query(&self) -> bool {
        matches!(self, SeqSrc::Query(_))
    }

    /// Returns whether the sequence is a reference sequence.
    #[inline]
    #[must_use]
    pub fn is_reference(&self) -> bool {
        matches!(self, SeqSrc::Reference(_))
    }
}

/// A trait providing convenience methods to wrap a reference to a sequence in
/// [`SeqSrc`] in preparation for alignment.
///
/// See [`SeqSrc`] for an example and more details.
pub trait AsSeqSrc {
    /// Wraps a reference to the sequence in [`SeqSrc::Query`].
    ///
    /// To wrap the owned type, use [`into_query_src`]. See [`SeqSrc`] for an
    /// example and more details.
    ///
    /// [`into_query_src`]: IntoSeqSrc::into_query_src
    #[inline]
    #[must_use]
    fn as_query_src(&self) -> SeqSrc<&Self> {
        SeqSrc::Query(self)
    }

    /// Wraps a reference to the sequence in [`SeqSrc::Reference`].
    ///
    /// To wrap the owned type, use [`into_ref_src`]. See [`SeqSrc`] for an
    /// example and more details.
    ///
    /// [`into_ref_src`]: IntoSeqSrc::into_ref_src
    #[inline]
    #[must_use]
    fn as_ref_src(&self) -> SeqSrc<&Self> {
        SeqSrc::Reference(self)
    }
}

impl<T: AsRef<[u8]>> AsSeqSrc for T {}

/// A trait providing convenience methods to wrap a sequence in [`SeqSrc`] in
/// preparation for alignment.
///
/// See [`SeqSrc`] for an example and more details.
pub trait IntoSeqSrc: Sized {
    /// Wraps the sequence in [`SeqSrc::Query`].
    ///
    /// To wrap a reference to the sequence instead, use [`as_query_src`]. See
    /// [`SeqSrc`] for an example and more details.
    ///
    /// [`as_query_src`]: AsSeqSrc::as_query_src
    #[inline]
    #[must_use]
    fn into_query_src(self) -> SeqSrc<Self> {
        SeqSrc::Query(self)
    }

    /// Wraps the sequence in [`SeqSrc::Reference`].
    ///
    /// To wrap a reference to the sequence instead, use [`as_ref_src`]. See
    /// [`SeqSrc`] for an example and more details.
    ///
    /// [`as_ref_src`]: AsSeqSrc::as_ref_src
    #[inline]
    #[must_use]
    fn into_ref_src(self) -> SeqSrc<Self> {
        SeqSrc::Reference(self)
    }
}

impl<T: AsRef<[u8]>> IntoSeqSrc for T {}
