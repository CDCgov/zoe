use crate::{
    alignment::phmm::{CorePhmm, PrecomputedDomainModule, PrecomputedLocalModule, SemiLocalModule},
    data::views::IndexAdjustable,
};
use std::ops::{Bound, Range};

/// A trait for structures that can be indexed via a [`PhmmIndex`], such as
/// pHMMs and modules.
pub(crate) trait PhmmIndexable: Sized {
    /// Returns the number of pseudomatch states in the pHMM-related structure
    /// (either a match state for the reference, or BEGIN or END).
    fn num_pseudomatch(&self) -> usize;

    /// Returns the length of the reference which the pHMM-related structure
    /// represents.
    ///
    /// This is also the number of match states in the pHMM-related structure,
    /// excluding BEGIN and END.
    #[inline]
    fn seq_len(&self) -> usize {
        self.num_pseudomatch() - 2
    }

    /// Returns the [`PhmmIndex`] as a dynamic programming index.
    #[inline]
    fn get_dp_index(&self, j: impl PhmmIndex) -> usize {
        j.get_phmm_dp_index(self)
    }

    /// Returns the [`PhmmIndex`] as a sequence index (with respect to the
    /// reference represented by the pHMM).
    ///
    /// If the index corresponds to the BEGIN state, then `None` is returned
    /// since this does not correspond to a position in the reference.
    fn get_seq_index(&self, j: impl PhmmIndex) -> Option<usize> {
        self.get_dp_index(j).checked_sub(1)
    }

    /// Converts the [`PhmmIndex`] to a sequence index (with respect to the
    /// reference represented by the pHMM).
    ///
    /// If the index corresponds to the BEGIN state, then `None` is returned
    /// since this does not correspond to a position in the reference.
    fn to_seq_index(&self, j: impl PhmmIndex) -> Option<SeqIndex> {
        self.get_seq_index(j).map(SeqIndex)
    }

    /// Returns a range of dynamic programming indices from a start and end
    /// [`PhmmIndex`].
    fn get_dp_range(&self, start: Bound<impl PhmmIndex>, end: Bound<impl PhmmIndex>) -> Range<usize> {
        let start = match start {
            Bound::Included(start) => self.get_dp_index(start),
            // +1 due to converting Excluded to Included
            Bound::Excluded(start) => self.get_dp_index(start) + 1,
            Bound::Unbounded => 0,
        };
        let end = match end {
            // +1 due to converting Included to Excluded
            Bound::Included(end) => self.get_dp_index(end) + 1,
            Bound::Excluded(end) => self.get_dp_index(end),
            Bound::Unbounded => self.seq_len(),
        };

        start..end
    }

    /// Returns a range of sequence indices from a start and end [`PhmmIndex`].
    ///
    ///
    /// If either index corresponds to the BEGIN state, then this will be mapped
    /// to 0 (the same sequence index that [`FirstMatch`] corresponds to).
    fn get_seq_range(&self, start: Bound<impl PhmmIndex>, end: Bound<impl PhmmIndex>) -> Range<usize> {
        // -1 for converting dynamic programming index to sequence index
        self.get_dp_range(start, end).saturating_sub(1)
    }
}

/// A trait for structures that can be indexed via a [`QueryIndex`], such as
/// query sequence.
pub(crate) trait QueryIndexable: Sized {
    /// Returns the length of the query sequence.
    fn seq_len(&self) -> usize;

    /// Returns the number of possible states that could be in when traversing
    /// the query.
    ///
    /// Either none of the bases are consumed, or up to and including
    /// `query_len` bases could be consumed.
    fn dp_len(&self) -> usize {
        self.seq_len() + 1
    }

    /// Returns the [`QueryIndex`] as a dynamic programming index.
    fn get_dp_index(&self, i: impl QueryIndex) -> usize {
        i.get_query_dp_index(self)
    }

    /// Converts the [`QueryIndex`] to a dynamic programming index.
    fn to_dp_index(&self, i: impl QueryIndex) -> DpIndex {
        DpIndex(self.get_dp_index(i))
    }

    /// Returns a range of dynamic programming indices from a start and end
    /// [`QueryIndex`].
    fn get_dp_range(&self, start: Bound<impl QueryIndex>, end: Bound<impl QueryIndex>) -> Range<usize> {
        let start = match start {
            Bound::Included(start) => self.get_dp_index(start),
            Bound::Excluded(start) => self.get_dp_index(start) + 1,
            Bound::Unbounded => 0,
        };
        let end = match end {
            Bound::Included(end) => self.get_dp_index(end) + 1,
            Bound::Excluded(end) => self.get_dp_index(end),
            Bound::Unbounded => self.seq_len(),
        };

        start..end
    }

    /// Returns a range of sequence indices from a start and end [`QueryIndex`].
    ///
    /// If either index corresponds [`NoBases`], then this will be mapped to 0
    /// (the same sequence index that [`FirstBase`] corresponds to).
    fn get_seq_range(&self, start: Bound<impl QueryIndex>, end: Bound<impl QueryIndex>) -> Range<usize> {
        // -1 for converting dynamic programming index to sequence index
        self.get_dp_range(start, end).saturating_sub(1)
    }
}

/// A trait representing different ways to index into a pHMM-related data
/// structure.
///
/// By indexing data structures with a [`PhmmIndex`], it allows for better
/// readability and correctness by requiring specification of the type of index
/// (such as a dynamic programming index with [`DpIndex`] or a sequence index
/// with [`SeqIndex`]). It also allows special elements to be accessed with
/// [`Begin`], [`FirstMatch`], [`LastMatch`], and [`End`].
pub(crate) trait PhmmIndex: IndexOffset {
    /// Helper function for [`PhmmIndexable::get_dp_index`], allowing each index
    /// type to control how it gets coverted to a dynamic programming index
    fn get_phmm_dp_index(&self, v: &impl PhmmIndexable) -> usize;

    /// Hook to allow the `End` index literal to be detected separately than the
    /// rest (e.g., for [`get_layer`] in [`CorePhmm`])
    ///
    /// [`get_layer`]: CorePhmm::get_layer
    #[inline]
    fn is_end(&self) -> bool {
        false
    }
}

/// A trait representing different ways to index into a query-related data
/// structure.
///
/// By indexing data structures with a [`QueryIndex`], it allows for better
/// readability and correctness by requiring specification of the type of index
/// (such as a dynamic programming index with [`DpIndex`] or a sequence index
/// with [`SeqIndex`]). It also allows special elements to be accessed with
/// [`NoBases`], [`FirstBase`], and [`LastBase`].
pub(crate) trait QueryIndex: IndexOffset {
    /// Helper function for [`QueryIndexable::get_dp_index`], allowing each
    /// index type to control how it gets converted to a dynamic programming
    /// index.
    fn get_query_dp_index(&self, v: &impl QueryIndexable) -> usize;
}

/// A [`PhmmIndex`] or [`QueryIndex`] representing an index with respect to the
/// sequence, rather than with respect to the dynamic programming tables.
///
/// 0 represents the first position in the sequence (reference or query).
#[repr(transparent)]
#[derive(Clone, Copy)]
pub(crate) struct SeqIndex(pub usize);

/// A [`PhmmIndex`] or [`QueryIndex`] representing an index with respect to the
/// dynamic programming tables, rather than with respect to the sequence itself.
///
/// 1 represents the first position in the sequence (reference or query), while
/// 0 represents matching no bases or the BEGIN state of the pHMM.
#[repr(transparent)]
#[derive(Clone, Copy)]
pub(crate) struct DpIndex(pub(crate) usize);

/// A [`PhmmIndex`] representing the BEGIN state of the pHMM, as well as the
/// first layer (which holds the transition probabilities out of the BEGIN
/// state).
#[derive(Clone, Copy)]
pub(crate) struct Begin;

/// A [`PhmmIndex`] representing the first match state of the pHMM after BEGIN.
///
/// This corresponds to the first residue in the reference sequences.
#[derive(Clone, Copy)]
pub(crate) struct FirstMatch;

/// A [`PhmmIndex`] representing the last match state of the pHMM before END.
///
/// This corresponds to the last residue in the reference sequences.
#[derive(Clone, Copy)]
pub(crate) struct LastMatch;

/// A [`PhmmIndex`] representing the END state of the pHMM.
///
/// When used to get a layer from a pHMM with [`get_layer`], this is treated as
/// the same thing as [`LastMatch`] since the END state does not have its own
/// layer.
///
/// [`get_layer`]: CorePhmm::get_layer
#[derive(Clone, Copy)]
pub(crate) struct End;

/// A [`QueryIndex`] representing not matching any bases from the sequence
/// (dynamic programming index 0).
#[derive(Clone, Copy)]
#[allow(dead_code)]
pub(crate) struct NoBases;

/// A [`QueryIndex`] representing the first base in the sequence.
#[derive(Clone, Copy)]
#[allow(dead_code)]
pub(crate) struct FirstBase;

/// A [`QueryIndex`] representing the last base in the sequence.
#[derive(Clone, Copy)]
pub(crate) struct LastBase;

impl<T, const S: usize> PhmmIndexable for CorePhmm<T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        // The END state does not have an index in the CorePhmm, so we add 1
        self.0.len() + 1
    }
}

impl<T> PhmmIndexable for SemiLocalModule<T> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.0.len()
    }
}

impl QueryIndexable for &[u8] {
    fn seq_len(&self) -> usize {
        self.len()
    }
}

impl<T, const S: usize> PhmmIndexable for PrecomputedLocalModule<'_, T, S> {
    #[inline]
    fn num_pseudomatch(&self) -> usize {
        self.external_params.num_pseudomatch()
    }
}

impl<T, const S: usize> QueryIndexable for PrecomputedDomainModule<T, S> {
    #[inline]
    fn seq_len(&self) -> usize {
        self.0.len() - 1
    }
}

impl PhmmIndex for DpIndex {
    #[inline]
    fn get_phmm_dp_index(&self, _v: &impl PhmmIndexable) -> usize {
        self.0
    }
}

impl PhmmIndex for SeqIndex {
    #[inline]
    fn get_phmm_dp_index(&self, _v: &impl PhmmIndexable) -> usize {
        self.0 + 1
    }
}

impl PhmmIndex for Begin {
    #[inline]
    fn get_phmm_dp_index(&self, _v: &impl PhmmIndexable) -> usize {
        0
    }
}

impl PhmmIndex for FirstMatch {
    #[inline]
    fn get_phmm_dp_index(&self, _v: &impl PhmmIndexable) -> usize {
        1
    }
}

impl PhmmIndex for LastMatch {
    #[inline]
    fn get_phmm_dp_index(&self, v: &impl PhmmIndexable) -> usize {
        v.num_pseudomatch() - 2
    }
}

impl PhmmIndex for End {
    #[inline]
    fn get_phmm_dp_index(&self, v: &impl PhmmIndexable) -> usize {
        v.num_pseudomatch() - 1
    }

    #[inline]
    fn is_end(&self) -> bool {
        true
    }
}

impl QueryIndex for NoBases {
    #[inline]
    fn get_query_dp_index(&self, _v: &impl QueryIndexable) -> usize {
        0
    }
}

impl QueryIndex for FirstBase {
    #[inline]
    fn get_query_dp_index(&self, _v: &impl QueryIndexable) -> usize {
        1
    }
}

impl QueryIndex for SeqIndex {
    #[inline]
    fn get_query_dp_index(&self, _v: &impl QueryIndexable) -> usize {
        self.0 + 1
    }
}

impl QueryIndex for DpIndex {
    #[inline]
    fn get_query_dp_index(&self, _v: &impl QueryIndexable) -> usize {
        self.0
    }
}

impl QueryIndex for LastBase {
    #[inline]
    fn get_query_dp_index(&self, v: &impl QueryIndexable) -> usize {
        // The last value in a slice of length dp_len
        v.dp_len() - 1
    }
}

/// A trait enabling a [`PhmmIndex`] or [`QueryIndex`] to be offset in either
/// direction via the type system, without converting it to [`DpIndex`].
pub(crate) trait IndexOffset: Copy {
    /// Get the index before the current one.
    #[inline]
    fn prev_index(self) -> PrevOffset<Self> {
        PrevOffset(self, 1)
    }

    /// Get the index after the current one.
    #[inline]
    fn next_index(self) -> NextOffset<Self> {
        NextOffset(self, 1)
    }
}

impl IndexOffset for SeqIndex {}
impl IndexOffset for DpIndex {}
impl IndexOffset for Begin {}
impl IndexOffset for FirstMatch {}
impl IndexOffset for LastMatch {}
impl IndexOffset for End {}
impl IndexOffset for NoBases {}
impl IndexOffset for FirstBase {}
impl IndexOffset for LastBase {}
impl<I: IndexOffset> IndexOffset for PrevOffset<I> {}
impl<I: IndexOffset> IndexOffset for NextOffset<I> {}

/// A wrapper type for a [`PhmmIndex`] or [`QueryIndex`] for specifying the
/// layer before another one by a given offset.
#[derive(Clone, Copy)]
pub struct PrevOffset<I>(I, usize);

impl<I: PhmmIndex> PhmmIndex for PrevOffset<I> {
    #[inline]
    fn get_phmm_dp_index(&self, v: &impl PhmmIndexable) -> usize {
        self.0.get_phmm_dp_index(v) - self.1
    }
}

impl<I: QueryIndex> QueryIndex for PrevOffset<I> {
    #[inline]
    fn get_query_dp_index(&self, v: &impl QueryIndexable) -> usize {
        self.0.get_query_dp_index(v) - self.1
    }
}

/// A wrapper type for a [`PhmmIndex`] or [`QueryIndex`] for specifying the
/// layer after another one by a given offset.
#[derive(Clone, Copy)]
pub struct NextOffset<I>(I, usize);

impl<I: PhmmIndex> PhmmIndex for NextOffset<I> {
    #[inline]
    fn get_phmm_dp_index(&self, v: &impl PhmmIndexable) -> usize {
        self.0.get_phmm_dp_index(v) + self.1
    }
}

impl<I: QueryIndex> QueryIndex for NextOffset<I> {
    #[inline]
    fn get_query_dp_index(&self, v: &impl QueryIndexable) -> usize {
        self.0.get_query_dp_index(v) + self.1
    }
}
