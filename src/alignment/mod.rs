//! ## Functions for aligning sequence data.
//!
//! *Zoe* supports efficient local alignment for DNA, protein, or any other
//! sequence data via the Smith-Waterman algorithm. For generating the score
//! (useful for database searches and determining whether sequences are
//! related), this includes both a traditional implementation (see
//! [`sw_scalar_score`]) and a striped SIMD implementation (see
//! [`sw_simd_score`]). The function [`sw_scalar_alignment`] can be used to
//! generate the alignment itself.
//!
//! To calculate an alignment score or generate an alignment, follow these
//! steps. The first two are designed to be done at compile time, and the
//! subsequent two to be done at runtime.
//!
//! 1. Choose an alphabet. *Zoe* provides an easy interface to work with DNA,
//!    assuming the bases `ACGT` are case-insensitive and `N` is used as a
//!    catch-all. Otherwise, you can define your own alphabet using
//!    [`ByteIndexMap::new`].
//!
//! 2. Specify the [`WeightMatrix`] used for scoring matches and mismatches. For
//!    DNA, [`new_dna_matrix`] is a convenient constructor.
//!
//! 3. Build the query profile, with either [`ScalarProfile::new`] or
//!    [`StripedProfile::new`] (for SIMD). This step combines the query, matrix
//!    of weights, and gap open and gap extend penalties. This step also
//!    performs some basic checks, and if any of the inputs are invalid, a
//!    [`QueryProfileError`] is returned.
//!
//! 4. Use the query profile to align against any number of different
//!    references. Simply call the [`ScalarProfile::smith_waterman_score`] or
//!    [`StripedProfile::smith_waterman_score`] methods. The profile can be
//!    reused in multiple different calls.
//!
//! Below is an example using DNA. We use a match score of 4 and a mismatch
//! score of -2, as defined in `WEIGHTS`. We also choose to use a
//! `StripedProfile` (so that the SIMD algorithm is used) with integer type `i8`
//! and 32 lanes.
//! ```
//! # use zoe::{alignment::StripedProfile, data::WeightMatrix};
//! let reference: &[u8] = b"GGCCACAGGATTGAG";
//! let query: &[u8] = b"CTCAGATTG";
//!
//! const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
//! const GAP_OPEN: i8 = -3;
//! const GAP_EXTEND: i8 = -1;
//!
//! let profile = StripedProfile::<i8, 32, 5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let score = profile.smith_waterman_score(reference).unwrap();
//! assert_eq!(score, 27);
//! ```
//!
//! Below is an example using a different alphabet. Matches are given a score of
//! 1 and mismatches are given a score of -1.
//! ```
//! # use zoe::{
//! #     alignment::StripedProfile,
//! #     data::{WeightMatrix, ByteIndexMap},
//! # };
//! let reference: &[u8] = b"BDAACAABDDDB";
//! let query: &[u8] = b"AABDDAB";
//!
//! const MAPPING: ByteIndexMap<4> = ByteIndexMap::new(*b"ABCD", b'A');
//! const WEIGHTS: WeightMatrix<i8, 4> = WeightMatrix::new(&MAPPING, 1, -1, None);
//! const GAP_OPEN: i8 = -4;
//! const GAP_EXTEND: i8 = -2;
//!
//! let profile = StripedProfile::<i8, 32, 4>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let score = profile.smith_waterman_score(reference).unwrap();
//! assert_eq!(score, 5);
//! ```
//!
//! When using the SIMD algorithm, you must specify the number of lanes `N` and
//! integer type `T` (typically i8, i16, i32, or i64). If the integer type has
//! too small of a range, it is possible the alignment score will overflow (in
//! which case `None` is returned). A higher-level abstraction to avoid this is
//! [`LocalProfiles`] and [`SharedProfiles`].
//!
//! The former is designed for use within a single thread, while the latter
//! allows multiple threads to access it. Both store a set of lazily-evaluated
//! query profiles for `i8`, `i16`, `i32`, and `i64`. To create one of these,
//! you can call [`LocalProfiles::new_with_w256`],
//! [`SharedProfiles::new_with_w256`], or one of the other constructors. Then,
//! call [`LocalProfiles::smith_waterman_score_from_i8`],
//! [`SharedProfiles::smith_waterman_score_from_i8`], or one of the other
//! methods.
//!
//! When using DNA, you can also create a profile by using
//! [`Nucleotides::into_local_profile`] or [`Nucleotides::into_shared_profile`].
//!
//! Below is the previous example using DNA, but with this higher-level API:
//! ```
//! # use zoe::{data::WeightMatrix, prelude::Nucleotides};
//! let reference: &[u8] = b"GGCCACAGGATTGAG";
//! let query: Nucleotides = b"CTCAGATTG".into();
//!
//! const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
//! const GAP_OPEN: i8 = -3;
//! const GAP_EXTEND: i8 = -1;
//!
//! let profile = query.into_local_profile(&WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let score = profile.smith_waterman_score_from_i8(reference).unwrap();
//! assert_eq!(score, 27);
//! ```
//!
//! And similarly with the different alphabet:
//! ```
//! # use zoe::{
//! #     alignment::LocalProfiles,
//! #     data::{WeightMatrix, ByteIndexMap},
//! # };
//! let reference: &[u8] = b"BDAACAABDDDB";
//! let query: &[u8] = b"AABDDAB";
//!
//! const MAPPING: ByteIndexMap<4> = ByteIndexMap::new(*b"ABCD", b'A');
//! const WEIGHTS: WeightMatrix<i8, 4> = WeightMatrix::new(&MAPPING, 1, -1, None);
//! const GAP_OPEN: i8 = -4;
//! const GAP_EXTEND: i8 = -2;
//!
//! let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let score = profile.smith_waterman_score_from_i8(reference).unwrap();
//! assert_eq!(score, 5);
//! ```
//!
//! [`sw_scalar_score`]: sw::sw_scalar_score
//! [`sw_simd_score`]: sw::sw_simd_score
//! [`sw_scalar_alignment`]: sw::sw_scalar_alignment
//! [`ByteIndexMap::new`]: crate::data::ByteIndexMap::new
//! [`WeightMatrix`]: crate::data::WeightMatrix
//! [`WeightMatrix::new`]: crate::data::WeightMatrix::new
//! [`into_biased_matrix`]: crate::data::WeightMatrix::into_biased_matrix
//! [`new_biased_dna_matrix`]: crate::data::WeightMatrix::new_biased_dna_matrix
//! [`new_dna_matrix`]: crate::data::WeightMatrix::new_dna_matrix
//! [`QueryProfileError`]: crate::data::err::QueryProfileError
//! [`Nucleotides::into_local_profile`]:
//!     crate::data::types::nucleotides::Nucleotides::into_local_profile
//! [`Nucleotides::into_shared_profile`]:
//!     crate::data::types::nucleotides::Nucleotides::into_shared_profile
//! [`ScalarProfile::new`]: ScalarProfile::new
//! [`StripedProfile::new`]: StripedProfile::new
//! [`ScalarProfile::smith_waterman_score`]: ScalarProfile::smith_waterman_score
//! [`StripedProfile::smith_waterman_score`]:
//!    StripedProfile::smith_waterman_score
//! [`LocalProfiles`]: LocalProfiles
//! [`SharedProfiles`]: SharedProfiles
//! [`LocalProfiles::smith_waterman_score_from_i8`]:
//!     LocalProfiles::smith_waterman_score_from_i8
//! [`SharedProfiles::smith_waterman_score_from_i8`]:
//!     SharedProfiles::smith_waterman_score_from_i8

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

// For the `Alignment` struct below, both ranges are 0-based and end-exclusive.
// For a global alignment, the ranges will each be encompass the full length of
// the sequences. For local alignment, `query_range` does NOT include hard or
// soft clipped bases, even though `cigar` does contain this information.

/// The output of an alignment algorithm
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

    // Convenience function for coercing to scores to `u64` for fuzzing and
    // testing.
    #[inline]
    #[must_use]
    #[cfg(any(feature = "fuzzing", test))]
    pub fn as_u64(&self) -> Alignment<u64>
    where
        T: crate::math::AnyInt, {
        Alignment {
            score:       self.score.cast_as::<u64>(),
            ref_range:   self.ref_range.clone(),
            query_range: self.query_range.clone(),
            states:      self.states.clone(),
            ref_len:     self.ref_len,
            query_len:   self.query_len,
        }
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
        let alignment = sw_scalar_alignment(reference, &profile);
        let invert_alignment = alignment.invert();

        assert_eq!(alignment.ref_range, 3..15);
        assert_eq!(alignment.query_range, 1..13);
        assert_eq!(alignment.states, Cigar::from_slice_unchecked("1S5M1D4M1I2M3S"));

        assert_eq!(invert_alignment.ref_range, 1..13);
        assert_eq!(invert_alignment.query_range, 3..15);
        assert_eq!(invert_alignment.states, Cigar::from_slice_unchecked("3S5M1I4M1D2M1S"));
    }
}
