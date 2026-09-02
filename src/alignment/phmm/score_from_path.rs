//! Functions for calculating the score of a given alignment in a pHMM.

use crate::alignment::{
    AlignmentStates,
    phmm::{
        DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
        traverse::alignment::{
            DomainAlignmentVisitor, DomainTraverseFromAlignError, GlobalAlignmentVisitor, GlobalTraverseFromAlignError,
            LocalAlignmentVisitor, LocalTraverseFromAlignError, SemiLocalAlignmentVisitor, SemiLocalTraverseFromAlignError,
        },
    },
};

impl<T: PhmmNumber, const S: usize> GlobalPhmm<T, S>
where
    T: 'static,
{
    /// Gets the score for a particular alignment.
    ///
    /// This is designed to give the exact same score as [`viterbi`] when the
    /// best `alignment` is passed, performing all arithmetic operations in the
    /// same order so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// If the alignment and/or starting coordinate are not valid for the pHMM,
    /// then an error describing the cause is returned.
    ///
    /// [`viterbi`]: GlobalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(
        &self, seq: Q, states: &AlignmentStates,
    ) -> Result<T, GlobalTraverseFromAlignError> {
        self.traverse(GlobalAlignmentVisitor::new(seq.as_ref(), states, self))
    }
}

impl<T: PhmmNumber, const S: usize> LocalPhmm<T, S>
where
    T: 'static,
{
    /// Gets the best score for a particular alignment.
    ///
    /// The starting 0-based index of the alignment within reference coordinates
    /// is passed. If the alignment is empty (i.e., all soft-clipping), then
    /// `ref_start` is ignored.
    ///
    /// The score may be infinite if no paths have a nonzero probability. This
    /// is designed to give the exact same score as [`viterbi`] when the best
    /// `path` is passed, performing all arithmetic operations in the same order
    /// so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// If the alignment and/or starting coordinate are not valid for the pHMM,
    /// then an error describing the cause is returned.
    ///
    /// [`viterbi`]: LocalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(
        &self, seq: Q, alignment: &AlignmentStates, ref_start: usize,
    ) -> Result<T, LocalTraverseFromAlignError> {
        self.traverse(LocalAlignmentVisitor::new(seq.as_ref(), alignment, ref_start, self)?)
    }
}

impl<T: PhmmNumber, const S: usize> DomainPhmm<T, S>
where
    T: 'static,
{
    /// Gets the best score for a particular alignment.
    ///
    /// The score may be infinite if no paths have a nonzero probability. This
    /// is designed to give the exact same score as [`viterbi`] when the best
    /// `path` is passed, performing all arithmetic operations in the same order
    /// so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// If the alignment and/or starting coordinate are not valid for the pHMM,
    /// then an error describing the cause is returned.
    ///
    /// [`viterbi`]: DomainPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(
        &self, seq: Q, states: &AlignmentStates,
    ) -> Result<T, DomainTraverseFromAlignError> {
        self.traverse(DomainAlignmentVisitor::new(seq.as_ref(), states, self)?)
    }
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<T, S>
where
    T: 'static,
{
    /// Gets the best score for a particular path.
    ///
    /// The starting 0-based index of the alignment within reference coordinates
    /// is passed. If the alignment is empty (i.e., the query sequence is empty
    /// and as such there are no states in the alignment), then `ref_start` is
    /// ignored.
    ///
    /// The score may be infinite if no paths have a nonzero probability. This
    /// is designed to give the exact same score as [`viterbi`] when the best
    /// `path` is passed, performing all arithmetic operations in the same order
    /// so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// If the alignment and/or starting coordinate are not valid for the pHMM,
    /// then an error describing the cause is returned.
    ///
    /// [`viterbi`]: LocalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(
        &self, seq: Q, states: &AlignmentStates, ref_start: usize,
    ) -> Result<T, SemiLocalTraverseFromAlignError> {
        self.traverse(SemiLocalAlignmentVisitor::new(seq.as_ref(), states, ref_start, self)?)
    }
}
