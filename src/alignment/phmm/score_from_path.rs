//! Functions for calculating the score of a given alignment in a pHMM.
//!
//! <div class="warning note">
//!
//! **Note**
//!
//! You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
//! use these functions.
//!
//! </div>

use crate::alignment::{
    AlignmentStates,
    phmm::{DomainPhmm, GlobalPhmm, LocalPhmm, PhmmError, PhmmNumber, SemiLocalPhmm},
};
use std::ops::Range;

impl<T: PhmmNumber, const S: usize> GlobalPhmm<T, S> {
    /// Get the score for a particular alignment.
    ///
    /// This is designed to give the exact same score as [`viterbi`] when the
    /// best `alignment` is passed, performing all arithmetic operations in the
    /// same order so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are M, =, X, I, and D.
    ///
    /// [`viterbi`]: GlobalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(&self, seq: Q, alignment: &AlignmentStates) -> Result<T, PhmmError> {
        self.visit_params(seq, alignment, |_| {})
    }
}

impl<T: PhmmNumber, const S: usize> LocalPhmm<T, S> {
    /// Get the best score for a particular alignment.
    ///
    /// The score may be infinite if no paths have a nonzero probability. This
    /// is designed to give the exact same score as [`viterbi`] when the best
    /// `path` is passed, performing all arithmetic operations in the same order
    /// so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are M, =, X, I, and D. If `path` does not
    /// correspond to a valid path through a [`LocalPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    ///
    /// [`viterbi`]: LocalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(
        &self, seq: Q, alignment: &AlignmentStates, ref_range: Range<usize>,
    ) -> Result<T, PhmmError> {
        self.visit_params(seq, alignment, ref_range, |_| {})
    }
}

impl<T: PhmmNumber, const S: usize> DomainPhmm<T, S> {
    /// Get the best score for a particular alignment.
    ///
    /// The score may be infinite if no paths have a nonzero probability. This
    /// is designed to give the exact same score as [`viterbi`] when the best
    /// `path` is passed, performing all arithmetic operations in the same order
    /// so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are M, =, X, I, and D. If `path` does not
    /// correspond to a valid path through a [`DomainPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    ///
    /// [`viterbi`]: DomainPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(&self, seq: Q, alignment: &AlignmentStates) -> Result<T, PhmmError> {
        self.visit_params(seq, alignment, |_| {})
    }
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<T, S> {
    /// Get the best score for a particular path. The score may be infinite if
    /// no paths have a nonzero probability.
    ///
    /// This is designed to give the exact same score as [`viterbi`] when the
    /// best `path` is passed, performing all arithmetic operations in the same
    /// order so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are M, =, X, I, and D. If `path` does not
    /// correspond to a valid path through a [`LocalPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    ///
    /// [`viterbi`]: LocalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(
        &self, seq: Q, alignment: &AlignmentStates, ref_range: Range<usize>,
    ) -> Result<T, PhmmError> {
        self.visit_params(seq, alignment, ref_range, |_| {})
    }
}
