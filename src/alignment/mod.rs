//! ## Functions for aligning sequence data.
//!
//! *Zoe* supports alignment for DNA, protein, or any other
//! sequence data.
//!
//! - [Smith-Waterman]: Optimal local alignment in the [`sw`] module.
//!
//! [Smith-Waterman]: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

#[cfg(feature = "dev-phmm")]
pub mod phmm;
pub mod sw;

mod errors;
mod profile;
mod profile_set;
mod state;
mod std_traits;

pub use errors::*;
pub use profile::*;
pub use profile_set::*;
pub use state::*;

use crate::data::cigar::Ciglet;
use std::ops::Range;

/// The output of an alignment algorithm.
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
    /// Unwraps the alignment, consuming the [`MaybeAligned<T>`] and returning the
    /// contained [`Alignment<T>`].
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

    /// Returns the alignment if it did not overflow, otherwise calls `f` and
    /// returns the result.
    #[inline]
    #[must_use]
    pub fn or_else_overflowed(self, f: impl FnOnce() -> Self) -> Self {
        if let MaybeAligned::Overflowed = self { f() } else { self }
    }

    /// Maps a `MaybeAligned<T>` to `MaybeAligned<U>` by applying a function to the contained value.
    /// Leaves `Overflowed` and `Unmapped` variants unchanged.
    #[inline]
    #[must_use]
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> MaybeAligned<U> {
        match self {
            MaybeAligned::Some(value) => MaybeAligned::Some(f(value)),
            MaybeAligned::Overflowed => MaybeAligned::Overflowed,
            MaybeAligned::Unmapped => MaybeAligned::Unmapped,
        }
    }

    /// Returns `Unmapped` if the `MaybeAligned` is `Unmapped`, `Overflowed` if it is
    /// `Overflowed`, otherwise calls `f` with the wrapped value and returns the
    /// result.
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

// For the `Alignment` struct below, both ranges are 0-based and end-exclusive.
// For a global alignment, the ranges will each be encompass the full length of
// the sequences. For local alignment, `query_range` does NOT include hard or
// soft clipped bases, even though `cigar` does contain this information.

/// A struct representing the information for an alignment, such as its score
/// and where in the sequences it occurs.
#[non_exhaustive]
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Alignment<T> {
    /// The score of the alignment
    pub score:       T,
    /// A range for the indices of the reference that are included in the
    /// alignment
    pub ref_range:   Range<usize>,
    /// A range for the indices of the query that are included in the alignment,
    /// already **excludes** clipped bases
    pub query_range: Range<usize>,
    /// States mapping reference and query.
    pub states:      AlignmentStates,
    /// The length of the reference sequence
    pub ref_len:     usize,
    /// The length of the query sequence
    pub query_len:   usize,
}

impl<T> Alignment<T> {
    /// Given the output of an alignment algorithm, generate the aligned
    /// sequences, using `-` as a gap character. The first output is the
    /// reference, and the second is the query.
    #[inline]
    #[must_use]
    pub fn get_aligned_seqs(&self, reference: &[u8], query: &[u8]) -> (Vec<u8>, Vec<u8>) {
        pairwise_align_with(reference, query, &self.states, self.ref_range.start)
    }

    /// Returns an iterator over the pairs of aligned bases in `reference` and
    /// `query`. `None` is used to represent gaps.
    ///
    /// This has equivalent behavior as [`get_aligned_seqs`] but as an iterator.
    /// [`get_aligned_seqs`] may offer faster performance, but
    /// [`get_aligned_iter`] avoids allocations and may offer a more convenient
    /// syntax for handling gaps.
    ///
    /// [`get_aligned_seqs`]: Alignment::get_aligned_seqs
    /// [`get_aligned_iter`]: Alignment::get_aligned_iter
    #[inline]
    #[must_use]
    pub fn get_aligned_iter<'a>(
        &'a self, reference: &'a [u8], query: &'a [u8],
    ) -> AlignmentIter<'a, std::iter::Copied<std::slice::Iter<'a, Ciglet>>> {
        AlignmentIter::new(reference, query, &self.states, self.ref_range.start)
    }
}

impl<T: Copy> Alignment<T> {
    /// Gets the alignment for when the query and reference are swapped.
    #[must_use]
    pub fn invert(&self) -> Self {
        let mut states = AlignmentStates::new();
        states.soft_clip(self.ref_range.start);
        let inverted_ciglets = self.states.0.iter().filter_map(|&ciglet| match ciglet.op {
            b'S' | b'H' => None,
            b'D' => Some(Ciglet {
                inc: ciglet.inc,
                op:  b'I',
            }),
            b'I' => Some(Ciglet {
                inc: ciglet.inc,
                op:  b'D',
            }),
            _ => Some(ciglet),
        });
        states.extend_from_ciglets(inverted_ciglets);
        states.soft_clip(self.ref_len - self.ref_range.end);

        Self {
            score: self.score,
            ref_range: self.query_range.clone(),
            query_range: self.ref_range.clone(),
            states,
            ref_len: self.query_len,
            query_len: self.ref_len,
        }
    }

    /// Gets the alignment for when the query and reference are both reversed
    /// (or reverse-complemented).
    #[must_use]
    pub fn to_reverse(&self) -> Self {
        let ref_range = (self.ref_len - self.ref_range.end)..(self.ref_len - self.ref_range.start - 1);
        let query_range = (self.query_len - self.query_range.end)..(self.query_len - self.query_range.start - 1);
        let states = self.states.to_reverse();

        Alignment {
            score: self.score,
            ref_range,
            query_range,
            states,
            ref_len: self.ref_len,
            query_len: self.query_len,
        }
    }

    /// Converts an alignment in-place so that it represents the alignment
    /// between the reversed (or reverse-complemented) query and reference.
    pub fn make_reverse(&mut self) {
        self.ref_range = (self.ref_len - self.ref_range.end)..(self.ref_len - self.ref_range.start - 1);
        self.query_range = (self.query_len - self.query_range.end)..(self.query_len - self.query_range.start - 1);
        self.states.make_reverse();
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{alignment::sw::sw_scalar_alignment, data::WeightMatrix, data::types::cigar::Cigar};

    #[test]
    fn alignment_invert() {
        const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
        const GAP_OPEN: i8 = -3;
        const GAP_EXTEND: i8 = -1;

        let reference: &[u8] = b"GGCCACAGGATTGAGC";
        let query: &[u8] = b"TCTCAGATTGCAGTTT";

        let profile = ScalarProfile::<5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        let alignment = sw_scalar_alignment(reference, &profile).unwrap();
        let invert_alignment = alignment.invert();

        assert_eq!(alignment.ref_range, 3..15);
        assert_eq!(alignment.query_range, 1..13);
        assert_eq!(alignment.states, Cigar::from_slice_unchecked("1S5M1D4M1I2M3S"));

        assert_eq!(invert_alignment.ref_range, 1..13);
        assert_eq!(invert_alignment.query_range, 3..15);
        assert_eq!(invert_alignment.states, Cigar::from_slice_unchecked("3S5M1I4M1D2M1S"));
    }
}
