//! Errors and context for alignment pHMM visitors.
//!
//! Traversing pHMMs based on alignments involves lots of edge cases and
//! potentials for error. For example, the CIGAR string may include too many or
//! too few reference-consuming or query-consuming operations. The start of the
//! alignment within the reference may be out of bounds. An early exit from an
//! INSERT or DELETE state may be attempted.
//!
//! Furthermore, a distinct visitor for each pHMM type is required due to the
//! nuances and ambiguities posed by each separate type. Some of the behavior
//! can be reused/shared such as [`CoreAlignmentVisitor`] and
//! [`AlignmentVisitorWithExit`], but error handling makes this difficult since
//! we would like the context to include the full queries/alignments (not just
//! the parts corresponding to the core pHMM).
//!
//! To resolve this, the innermost reusable part ([`CoreAlignmentVisitor`]) also
//! contains a generic context. The outermost visitors create this context and
//! construct the inner visitors with it. Traits on this such as
//! [`CoreContextToErr`] and [`CoreContextWithExitToErr`] provide methods for
//! constructing the different errors from the context.
//!
//! In this way, errors originating deep within the visitors can be constructed
//! directly with full context, rather than needing a complicated error-handling
//! scheme where errors are transformed/overwritten with each layer of
//! abstraction.
//!
//! [`CoreAlignmentVisitor`]:
//!     crate::alignment::phmm::traverse::alignment::CoreAlignmentVisitor
//! [`AlignmentVisitorWithExit`]:
//!     crate::alignment::phmm::traverse::alignment::AlignmentVisitorWithExit

use crate::{
    alignment::{
        AlignmentStates,
        phmm::{
            InvalidModelError,
            indexing::DpIndex,
            state::{PhmmState, PhmmStateOrModule},
        },
    },
    data::{cigar::LenInAlignment, err::GetCode},
};
use std::{error::Error, fmt::Display};

/// The possible errors when using an alignment to traverse a [`GlobalPhmm`].
///
/// [`GlobalPhmm`]: crate::alignment::phmm::GlobalPhmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum GlobalTraverseFromAlignError {
    /// An error caused by an invalid model.
    InvalidModel(InvalidModelError),
    /// Unsupported CIGAR operation outside of `MDI=X`.
    InvalidCigarOp(u8),
    /// The length of the reference as implied by the model and the alignment
    /// disagree.
    ModelLenMismatch {
        ref_len:              usize,
        ref_len_in_alignment: usize,
    },
    /// The length of the query as implied by the sequence and the alignment
    /// disagree.
    QueryLenMismatch {
        query_len:              usize,
        query_len_in_alignment: usize,
    },
}

/// The possible errors when using an alignment to traverse a [`SemiLocalPhmm`].
///
/// [`SemiLocalPhmm`]: crate::alignment::phmm::SemiLocalPhmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum SemiLocalTraverseFromAlignError {
    /// An error caused by an invalid model.
    InvalidModel(InvalidModelError),
    /// Unsupported CIGAR operation outside of `MDI=X`.
    InvalidCigarOp(u8),
    /// The length of the reference as implied by the model is shorter than
    /// what is implied by the alignment.
    ModelLenMismatch {
        ref_len:              usize,
        ref_start:            usize,
        ref_len_in_alignment: usize,
    },
    /// The length of the query as implied by the sequence and the alignment
    /// disagree.
    QueryLenMismatch {
        query_len:              usize,
        query_len_in_alignment: usize,
    },
    /// Due to a semilocal transition into a layer other than the BEGIN or END
    /// layer of the core pHMM, a match operation in the alignment was expected,
    /// but it was not found.
    MissingMatchOp {
        ref_len:   usize,
        ref_coord: usize,
        op:        u8,
    },
    /// After directly entering the END state from the [`SemiLocalModule`] at
    /// the beginning of the pHMM, an operation in the CIGAR string was
    /// encountered.
    ///
    /// [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule
    OpAfterEnteringEnd { op: u8 },
    /// The starting index of the alignment within reference coordinates was out
    /// of bounds for the model.
    RefStartOutOfBounds { ref_start: usize, ref_len: usize },
    /// The alignment ended while within an insert or delete state, rather than
    /// a match state as expected.
    InvalidEarlyExit {
        ref_len: usize,
        layer:   DpIndex,
        state:   PhmmState,
    },
}

/// The possible errors when using an alignment to traverse a [`DomainPhmm`].
///
/// [`DomainPhmm`]: crate::alignment::phmm::DomainPhmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum DomainTraverseFromAlignError {
    /// An error caused by an invalid model.
    InvalidModel(InvalidModelError),
    /// Unsupported CIGAR operation outside of `MDIS=X`.
    InvalidCigarOp(u8),
    /// Found soft clipping in the middle of the alignment.
    InternalClipping,
    /// Two adjacent ciglets have the same operation.
    DuplicateOp,
    /// The length of the reference as implied by the model and the alignment
    /// disagree.
    ModelLenMismatch {
        ref_len:              usize,
        ref_len_in_alignment: usize,
    },
    /// The length of the query as implied by the sequence and the alignment
    /// disagree.
    QueryLenMismatch {
        query_len:              usize,
        query_len_in_alignment: usize,
        skipped_start:          usize,
        skipped_end:            usize,
    },
}

/// The possible errors when using an alignment to traverse a [`LocalPhmm`].
///
/// [`LocalPhmm`]: crate::alignment::phmm::LocalPhmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum LocalTraverseFromAlignError {
    /// An error caused by an invalid model.
    InvalidModel(InvalidModelError),
    /// Unsupported CIGAR operation outside of `MDIS=X`.
    InvalidCigarOp(u8),
    /// Found soft clipping in the middle of the alignment.
    InternalClipping,
    /// Two adjacent ciglets have the same operation.
    DuplicateOp,
    /// The length of the reference as implied by the model is shorter than
    /// what is implied by the alignment.
    ModelLenMismatch {
        ref_len:              usize,
        ref_start:            usize,
        ref_len_in_alignment: usize,
    },
    /// The length of the query as implied by the sequence and the alignment
    /// disagree.
    QueryLenMismatch {
        query_len:              usize,
        query_len_in_alignment: usize,
        skipped_start:          usize,
        skipped_end:            usize,
    },
    /// Due to a semilocal transition into a layer other than the BEGIN or END
    /// layer of the core pHMM, a match operation in the alignment was expected,
    /// but it was not found.
    MissingMatchOp {
        ref_len:   usize,
        ref_coord: usize,
        op:        u8,
    },
    /// After directly entering the END state from the [`SemiLocalModule`] at
    /// the beginning of the pHMM, an operation in the CIGAR string was
    /// encountered.
    ///
    /// [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule
    OpAfterEnteringEnd { op: u8 },
    /// The starting index of the alignment within reference coordinates was out
    /// of bounds for the model.
    RefStartOutOfBounds { ref_start: usize, ref_len: usize },
    /// The alignment ended within an insert or delete state or attempted to
    /// exit early from one of these states, rather than a match state as
    /// expected.
    InvalidEarlyExit {
        ref_len: usize,
        layer:   DpIndex,
        state:   PhmmState,
    },
}

impl Display for GlobalTraverseFromAlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            GlobalTraverseFromAlignError::InvalidModel(err) => write!(f, "{err}"),
            GlobalTraverseFromAlignError::InvalidCigarOp(op) => write!(
                f,
                "An invalid CIGAR operation was found ({}). Expected one of MDI=X.",
                op as char
            ),
            GlobalTraverseFromAlignError::ModelLenMismatch {
                ref_len,
                ref_len_in_alignment,
            } => write!(
                f,
                "The length of the reference as implied by the pHMM is {ref_len}, but the alignment implies a length of {ref_len_in_alignment}"
            ),
            GlobalTraverseFromAlignError::QueryLenMismatch {
                query_len,
                query_len_in_alignment,
            } => write!(
                f,
                "The length of the query is {query_len}, but the alignment implies the length should be {query_len_in_alignment}"
            ),
        }
    }
}

impl Display for SemiLocalTraverseFromAlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            SemiLocalTraverseFromAlignError::InvalidModel(err) => write!(f, "{err}"),
            SemiLocalTraverseFromAlignError::InvalidCigarOp(op) => write!(
                f,
                "An invalid CIGAR operation was found ({}). Expected one of MDI=X.",
                op as char
            ),
            SemiLocalTraverseFromAlignError::ModelLenMismatch {
                ref_len,
                ref_start,
                ref_len_in_alignment,
            } => {
                if let Some(total_ref) = ref_start.checked_add(ref_len_in_alignment) {
                    write!(
                        f,
                        "The length of the reference as implied by the pHMM is {ref_len}, but the alignment implies a length of at least {total_ref} (starting at {ref_start} and consuming {ref_len_in_alignment} residues)",
                    )
                } else {
                    write!(
                        f,
                        "The length of the reference as implied by the pHMM is {ref_len}, but the alignment implies a longer length (starting at {ref_start} and consuming {ref_len_in_alignment} residues)"
                    )
                }
            }
            SemiLocalTraverseFromAlignError::QueryLenMismatch {
                query_len,
                query_len_in_alignment,
            } => write!(
                f,
                "The length of the query is {query_len}, but the alignment implies the length should be {query_len_in_alignment}"
            ),
            SemiLocalTraverseFromAlignError::MissingMatchOp { ref_len, ref_coord, op } => write!(
                f,
                "The alignment enters the MATCH state corresponding to reference coordinate {ref_coord} (out of a total reference length of {ref_len}). Instead of finding a match, the operation {op} was found.",
                op = op as char
            ),
            SemiLocalTraverseFromAlignError::OpAfterEnteringEnd { op } => write!(
                f,
                "After directly entering the END state from the beginning of the pHMM, the operation {op} was found.",
                op = op as char
            ),
            SemiLocalTraverseFromAlignError::RefStartOutOfBounds { ref_start, ref_len } => write!(
                f,
                "The alignment starts at reference coordinate {ref_start}, but the total length of the reference is only {ref_len}"
            ),
            SemiLocalTraverseFromAlignError::InvalidEarlyExit { ref_len, layer, state } => {
                let ref_coord = layer.to_seq_index().map_or(0, |x| x.0);
                write!(
                    f,
                    "The alignment ends in the {state} state at reference coordinate {ref_coord} (out of a total reference length of {ref_len}), but an early exit is only permitted from a match state"
                )
            }
        }
    }
}

impl Display for DomainTraverseFromAlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            DomainTraverseFromAlignError::InvalidModel(err) => write!(f, "{err}"),
            DomainTraverseFromAlignError::InvalidCigarOp(op) => write!(
                f,
                "An invalid CIGAR operation was found ({}). Expected one of MDIS=X.",
                op as char
            ),
            DomainTraverseFromAlignError::InternalClipping => {
                write!(f, "Soft clipping was found in the middle of the alignment.")
            }
            DomainTraverseFromAlignError::DuplicateOp => {
                write!(f, "Two adjacent ciglets were found with duplicate operations.")
            }
            DomainTraverseFromAlignError::ModelLenMismatch {
                ref_len,
                ref_len_in_alignment,
            } => write!(
                f,
                "The length of the reference as implied by the pHMM is {ref_len}, but the alignment implies a length of {ref_len_in_alignment}"
            ),
            DomainTraverseFromAlignError::QueryLenMismatch {
                query_len,
                skipped_start,
                skipped_end,
                query_len_in_alignment,
            } => write!(
                f,
                "The length of the query is {query_len}, but the alignment implies the length should be {query_len_in_alignment} ({skipped_start} clipped at start, {skipped_end} clipped at end)"
            ),
        }
    }
}

impl Display for LocalTraverseFromAlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            LocalTraverseFromAlignError::InvalidModel(err) => write!(f, "{err}"),
            LocalTraverseFromAlignError::InvalidCigarOp(op) => write!(
                f,
                "An invalid CIGAR operation was found ({}). Expected one of MDIS=X.",
                op as char
            ),
            LocalTraverseFromAlignError::InternalClipping => {
                write!(f, "Soft clipping was found in the middle of the alignment.")
            }
            LocalTraverseFromAlignError::DuplicateOp => {
                write!(f, "Two adjacent ciglets were found with duplicate operations.")
            }
            LocalTraverseFromAlignError::ModelLenMismatch {
                ref_len,
                ref_start,
                ref_len_in_alignment,
            } => write!(
                f,
                "The length of the reference as implied by the pHMM is {ref_len}, but the alignment implies a length of at least {total_ref} (starting at {ref_start} and consuming {ref_len_in_alignment} residues)",
                total_ref = ref_start + ref_len_in_alignment
            ),
            LocalTraverseFromAlignError::QueryLenMismatch {
                query_len,
                query_len_in_alignment,
                skipped_start,
                skipped_end,
            } => write!(
                f,
                "The length of the query is {query_len}, but the alignment implies the length should be {query_len_in_alignment} ({skipped_start} clipped at start, {skipped_end} clipped at end)"
            ),
            LocalTraverseFromAlignError::MissingMatchOp { ref_len, ref_coord, op } => write!(
                f,
                "The alignment enters the MATCH state corresponding to reference coordinate {ref_coord} (out of a total reference length of {ref_len}). Instead of finding a match, the operation {op} was found.",
                op = op as char
            ),
            LocalTraverseFromAlignError::OpAfterEnteringEnd { op } => write!(
                f,
                "After directly entering the END state from the beginning of the pHMM, the operation {op} was found.",
                op = op as char
            ),
            LocalTraverseFromAlignError::RefStartOutOfBounds { ref_start, ref_len } => write!(
                f,
                "The alignment starts at reference coordinate {ref_start}, but the total length of the reference is only {ref_len}"
            ),
            LocalTraverseFromAlignError::InvalidEarlyExit { ref_len, layer, state } => {
                let ref_coord = layer.to_seq_index().map_or(0, |x| x.0);
                write!(
                    f,
                    "The alignment ends in the {state} state at reference coordinate {ref_coord} (out of a total reference length of {ref_len}), but an early exit is only permitted from a match state"
                )
            }
        }
    }
}

impl Error for GlobalTraverseFromAlignError {}
impl GetCode for GlobalTraverseFromAlignError {}
impl Error for SemiLocalTraverseFromAlignError {}
impl GetCode for SemiLocalTraverseFromAlignError {}
impl Error for DomainTraverseFromAlignError {}
impl GetCode for DomainTraverseFromAlignError {}
impl Error for LocalTraverseFromAlignError {}
impl GetCode for LocalTraverseFromAlignError {}

/// A trait unifying the context types for all pHMM types.
pub(super) trait CoreContextToErr {
    /// The type of error yielded by the pHMM's alignment visitor.
    type Error;

    /// Constructs an error for the length of the query as implied by the
    /// sequence and the alignment disagreeing.
    fn query_len_mismatch(&self) -> Self::Error;

    /// Constructs an error for the length of the reference as implied by the
    /// model and the alignment disagreeing.
    fn model_len_mismatch(&self) -> Self::Error;

    /// Constructs an error when a alignment operation is not legal for the
    /// given pHMM type (i.e., not corresponding to variants of [`PhmmState`]
    /// for global/semilocal or [`PhmmStateOrModule`] for domain/local).
    fn invalid_op(&self, op: u8) -> Self::Error;

    /// Constructs the appropriate error when `op` cannot be converted into
    /// [`PhmmState`] by [`choose_core_transition`].
    ///
    /// [`choose_core_transition`]:
    ///     crate::alignment::phmm::traverse::alignment::CoreAlignmentVisitor::choose_core_transition
    fn core_transition_op_error(&self, op: u8, layer: DpIndex, exiting: PhmmState) -> Self::Error;

    /// Constructs the appropriate error when the operation iterator runs out
    /// during [`choose_core_transition`].
    ///
    /// [`choose_core_transition`]:
    ///     crate::alignment::phmm::traverse::alignment::CoreAlignmentVisitor::choose_core_transition
    fn core_transition_no_op_error(&self, layer: DpIndex, exiting: PhmmState) -> Self::Error;

    /// Constructs the appropriate error when [`choose_end_or_insert`]
    /// encounters an operation that is not `I`.
    ///
    /// [`choose_end_or_insert`]:
    ///     crate::alignment::phmm::traverse::alignment::CoreAlignmentVisitor::choose_end_or_insert
    fn choose_end_or_insert_op_error(&self, op: u8, layer: DpIndex, exiting: PhmmState) -> Self::Error;
}

/// A trait unifying the context types for pHMMs with semilocal transitions
/// (semilocal and local pHMMs).
pub(super) trait CoreContextWithExitToErr: CoreContextToErr {
    /// Constructs an error when the operation following the last match state is
    /// not `I` or `S` (for `LocalAlignmentVisitor`).
    fn choose_end_insert_or_exit_op_error(&self, op: u8, layer: DpIndex) -> Self::Error;

    /// Constructs an error caused by a match operation not being present,
    /// despite entering into a layer other than the BEGIN or END layer from the
    /// starting module.
    fn no_match_after_enter(&self, op: u8) -> Self::Error;

    /// Constructs an error caused by the starting index of the alignment within
    /// reference coordinates being out of bounds for the model.
    fn ref_start_out_of_bounds(&self) -> Self::Error;
}

/// The context required to form a [`GlobalTraverseFromAlignError`].
pub(super) struct GlobalContext<'a> {
    /// The length of the query sequence.
    pub query_len: usize,
    /// The alignment used for traversal.
    pub states:    &'a AlignmentStates,
    /// The length of the reference as implied by the pHMM.
    pub ref_len:   usize,
}

/// The context required to form a [`SemiLocalTraverseFromAlignError`].
pub(super) struct SemiLocalContext<'a> {
    /// The length of the query sequence.
    pub query_len: usize,
    /// The alignment used for traversal.
    pub states:    &'a AlignmentStates,
    /// The length of the reference as implied by the pHMM.
    pub ref_len:   usize,
    /// The starting index of the alignment within reference coordinates.
    pub ref_start: usize,
}

/// The context required to form a [`DomainTraverseFromAlignError`].
pub(super) struct DomainContext<'a> {
    /// The length of the query sequence.
    pub query_len:     usize,
    /// The alignment used for traversal.
    pub states:        &'a AlignmentStates,
    /// The length of the reference as implied by the pHMM.
    pub ref_len:       usize,
    /// The number of bases in the query that were emitted by the module at the
    /// start of the pHMM.
    ///
    /// Most of the time this can be read from `states`, but in the case of an
    /// empty alignment, this is not possible.
    pub skipped_start: usize,
    /// The number of bases in the query that were emitted by the module at the
    /// end of the pHMM.
    ///
    /// Most of the time this can be read from `states`, but in the case of an
    /// empty alignment, this is not possible.
    pub skipped_end:   usize,
}

/// The context required to form a [`LocalTraverseFromAlignError`].
pub(super) struct LocalContext<'a> {
    /// The length of the query sequence.
    pub query_len:     usize,
    /// The alignment used for traversal.
    pub states:        &'a AlignmentStates,
    /// The length of the reference as implied by the pHMM.
    pub ref_len:       usize,
    /// The starting index of the alignment within reference coordinates.
    pub ref_start:     usize,
    /// The number of bases in the query that were emitted by the module at the
    /// start of the pHMM.
    ///
    /// Most of the time this can be read from `states`, but in the case of an
    /// empty alignment, this is not possible.
    pub skipped_start: usize,
    /// The number of bases in the query that were emitted by the module at the
    /// end of the pHMM.
    ///
    /// Most of the time this can be read from `states`, but in the case of an
    /// empty alignment, this is not possible.
    pub skipped_end:   usize,
}

impl CoreContextToErr for GlobalContext<'_> {
    type Error = GlobalTraverseFromAlignError;

    fn query_len_mismatch(&self) -> Self::Error {
        GlobalTraverseFromAlignError::QueryLenMismatch {
            query_len:              self.query_len,
            query_len_in_alignment: self.states.query_len_in_alignment(),
        }
    }

    fn model_len_mismatch(&self) -> Self::Error {
        GlobalTraverseFromAlignError::ModelLenMismatch {
            ref_len:              self.ref_len,
            ref_len_in_alignment: self.states.ref_len_in_alignment(),
        }
    }

    fn invalid_op(&self, op: u8) -> Self::Error {
        GlobalTraverseFromAlignError::InvalidCigarOp(op)
    }

    fn core_transition_op_error(&self, op: u8, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // A global pHMM only accepts MDI=X, which correspond exactly to
        // PhmmState, so core_transition_op_error is always an invalid_op error
        self.invalid_op(op)
    }

    fn core_transition_no_op_error(&self, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // GlobalAlignmentVisitor represents the full CIGAR string in its
        // operation iterator, and no early exit is permitted, so a missing
        // operation is always a model length mismatch
        self.model_len_mismatch()
    }

    fn choose_end_or_insert_op_error(&self, op: u8, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // Confirm that the operation is MDI=X, or return invalid operation
        if PhmmState::from_op(op).is_none() {
            return self.invalid_op(op);
        }

        // Found one of MD=X, which consumes another reference position.
        // However, since traversal is at the end of the model, this is a model
        // length mismatch
        self.model_len_mismatch()
    }
}

impl CoreContextToErr for SemiLocalContext<'_> {
    type Error = SemiLocalTraverseFromAlignError;

    fn query_len_mismatch(&self) -> Self::Error {
        SemiLocalTraverseFromAlignError::QueryLenMismatch {
            query_len:              self.query_len,
            query_len_in_alignment: self.states.query_len_in_alignment(),
        }
    }

    fn model_len_mismatch(&self) -> Self::Error {
        SemiLocalTraverseFromAlignError::ModelLenMismatch {
            ref_len:              self.ref_len,
            ref_start:            self.ref_start,
            ref_len_in_alignment: self.states.ref_len_in_alignment(),
        }
    }

    fn invalid_op(&self, op: u8) -> Self::Error {
        SemiLocalTraverseFromAlignError::InvalidCigarOp(op)
    }

    fn core_transition_op_error(&self, op: u8, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // A semilocal pHMM only accepts MDI=X, which correspond exactly to
        // PhmmState, so core_transition_op_error is always an invalid_op error
        self.invalid_op(op)
    }

    fn core_transition_no_op_error(&self, layer: DpIndex, exiting: PhmmState) -> Self::Error {
        // choose_core_transition is only called for a SemiLocalPhmm while in an
        // insert or delete state that is not in the last match layer, from
        // which early exit is not allowed
        SemiLocalTraverseFromAlignError::InvalidEarlyExit {
            ref_len: self.ref_len,
            layer,
            state: exiting,
        }
    }

    fn choose_end_or_insert_op_error(&self, op: u8, layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // For semilocal pHMMs, the only way early exiting is implied is via the
        // operation iterator ending, since soft clipping is not permitted.
        // However, the operation iterator did not end if this function is being
        // called, so it is equivalent to call
        // choose_end_insert_or_exit_op_error
        self.choose_end_insert_or_exit_op_error(op, layer)
    }
}

impl CoreContextWithExitToErr for SemiLocalContext<'_> {
    fn choose_end_insert_or_exit_op_error(&self, op: u8, _layer: DpIndex) -> Self::Error {
        // This function is only called when the operation is not I or S.
        // Confirm that the operation is in MD=X, or return invalid operation
        if PhmmState::from_op(op).is_none() {
            return self.invalid_op(op);
        }

        // Found one of MD=X, which consumes another reference position.
        // However, since traversal is at the end of the model, this is a model
        // length mismatch
        self.model_len_mismatch()
    }

    fn no_match_after_enter(&self, op: u8) -> Self::Error {
        SemiLocalTraverseFromAlignError::MissingMatchOp {
            ref_len: self.ref_len,
            ref_coord: self.ref_start,
            op,
        }
    }

    fn ref_start_out_of_bounds(&self) -> Self::Error {
        SemiLocalTraverseFromAlignError::RefStartOutOfBounds {
            ref_start: self.ref_start,
            ref_len:   self.ref_len,
        }
    }
}

impl CoreContextToErr for DomainContext<'_> {
    type Error = DomainTraverseFromAlignError;

    fn query_len_mismatch(&self) -> Self::Error {
        DomainTraverseFromAlignError::QueryLenMismatch {
            query_len:              self.query_len,
            skipped_start:          self.skipped_start,
            skipped_end:            self.skipped_end,
            query_len_in_alignment: self.states.query_len_in_alignment(),
        }
    }

    fn model_len_mismatch(&self) -> Self::Error {
        DomainTraverseFromAlignError::ModelLenMismatch {
            ref_len:              self.ref_len,
            ref_len_in_alignment: self.states.ref_len_in_alignment(),
        }
    }

    fn invalid_op(&self, op: u8) -> Self::Error {
        DomainTraverseFromAlignError::InvalidCigarOp(op)
    }

    fn core_transition_op_error(&self, op: u8, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // The operation is not MDI=X, so if it is also not S, then it is
        // invalid
        if op != b'S' {
            return self.invalid_op(op);
        }

        // DomainAlignmentVisitor strips the first and last clipping ciglet if
        // present during initialization. Hence, S is either internal clipping
        // or a duplicate adjacent clipping operation

        // Check for duplicate operation
        if self.states.as_slice().array_windows::<2>().any(|[c1, c2]| c1.op == c2.op) {
            return DomainTraverseFromAlignError::DuplicateOp;
        }

        // Only remaining possibility is internal clipping
        DomainTraverseFromAlignError::InternalClipping
    }

    fn core_transition_no_op_error(&self, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // Domain pHMMs do not allow early exit, so if the alignment ends
        // prematurely within the core pHMM, this is always a model length
        // mismatch
        self.model_len_mismatch()
    }

    fn choose_end_or_insert_op_error(&self, op: u8, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // Confirm that the operation is MDI=X, or return invalid operation
        if PhmmState::from_op(op).is_none() {
            return self.invalid_op(op);
        }

        // Found one of MD=X, which consumes another reference position.
        // However, since traversal is at the end of the model, this is a model
        // length mismatch
        self.model_len_mismatch()
    }
}

impl DomainContext<'_> {
    /// A helper function for generating the error when an extra operation is
    /// found in [`finalize`] for [`DomainAlignmentVisitor`].
    ///
    /// [`DomainAlignmentVisitor`]:
    ///     crate::alignment::phmm::traverse::alignment::DomainAlignmentVisitor
    /// [`finalize`]: crate::alignment::phmm::traverse::DomainVisitor::finalize
    pub(super) fn remaining_op_error(&self, op: u8) -> DomainTraverseFromAlignError {
        // Validate that the operation is in MDIS=X
        if PhmmStateOrModule::from_op(op).is_none() {
            return self.invalid_op(op);
        }

        // Confirm no duplicate operations, which could be the cause of an
        // unexpected soft clipping operation
        if self.states.as_slice().array_windows::<2>().any(|[c1, c2]| c1.op == c2.op) {
            return DomainTraverseFromAlignError::DuplicateOp;
        }

        // Confirm no internal soft clipping
        let mut inner_states = self.states.as_slice();
        inner_states.split_off_first();
        inner_states.split_off_last();
        if inner_states.iter().any(|ciglet| ciglet.op == b'S' && ciglet.inc > 0) {
            return DomainTraverseFromAlignError::InternalClipping;
        }

        match PhmmState::from_op(op) {
            Some(PhmmState::Insert) => self.query_len_mismatch(),
            Some(PhmmState::Match | PhmmState::Delete) => self.model_len_mismatch(),
            None => DomainTraverseFromAlignError::InternalClipping,
        }
    }
}

impl CoreContextToErr for LocalContext<'_> {
    type Error = LocalTraverseFromAlignError;

    fn query_len_mismatch(&self) -> Self::Error {
        LocalTraverseFromAlignError::QueryLenMismatch {
            query_len:              self.query_len,
            query_len_in_alignment: self.states.query_len_in_alignment(),
            skipped_start:          self.skipped_start,
            skipped_end:            self.skipped_end,
        }
    }

    fn model_len_mismatch(&self) -> Self::Error {
        LocalTraverseFromAlignError::ModelLenMismatch {
            ref_len:              self.ref_len,
            ref_start:            self.ref_start,
            ref_len_in_alignment: self.states.ref_len_in_alignment(),
        }
    }

    fn invalid_op(&self, op: u8) -> Self::Error {
        LocalTraverseFromAlignError::InvalidCigarOp(op)
    }

    fn core_transition_op_error(&self, op: u8, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // PhmmState::from_op rejected op, and for a local pHMM, S is the only
        // additional valid operation.
        if op != b'S' {
            return self.invalid_op(op);
        }

        // LocalAlignmentVisitor strips the first and last clipping ciglet if
        // present during initialization. Hence, S is either internal clipping
        // or a duplicate adjacent clipping operation

        // Check for duplicate operation
        if self.states.as_slice().array_windows::<2>().any(|[c1, c2]| c1.op == c2.op) {
            return LocalTraverseFromAlignError::DuplicateOp;
        }

        // Only remaining possibility is internal clipping
        LocalTraverseFromAlignError::InternalClipping
    }

    fn core_transition_no_op_error(&self, layer: DpIndex, exiting: PhmmState) -> Self::Error {
        // Either soft clipping occurred (which was stripped) or the alignment
        // ran out. Either case represents an invalid early exit from the core
        // pHMM
        LocalTraverseFromAlignError::InvalidEarlyExit {
            ref_len: self.ref_len,
            layer,
            state: exiting,
        }
    }

    fn choose_end_or_insert_op_error(&self, op: u8, _layer: DpIndex, _exiting: PhmmState) -> Self::Error {
        // Confirm that the operation is MDI=X, or return invalid operation
        if PhmmState::from_op(op).is_none() {
            return self.invalid_op(op);
        }

        // Found one of MD=X, which consumes another reference position.
        // However, since traversal is at the end of the model, this is a model
        // length mismatch
        self.model_len_mismatch()
    }
}

impl LocalContext<'_> {
    /// A helper function for generating the error when an extra operation is
    /// found in [`finalize`] for [`LocalAlignmentVisitor`].
    ///
    /// [`LocalAlignmentVisitor`]:
    ///     crate::alignment::phmm::traverse::alignment::LocalAlignmentVisitor
    /// [`finalize`]: crate::alignment::phmm::traverse::LocalVisitor::finalize
    pub(super) fn remaining_op_error(&self, op: u8) -> LocalTraverseFromAlignError {
        // Validate that the operation is in MDIS=X
        if PhmmStateOrModule::from_op(op).is_none() {
            return self.invalid_op(op);
        }

        // Confirm no duplicate operations, which could be the cause of an
        // unexpected soft clipping operation
        if self.states.as_slice().array_windows::<2>().any(|[c1, c2]| c1.op == c2.op) {
            return LocalTraverseFromAlignError::DuplicateOp;
        }

        // Confirm no internal soft clipping
        let mut inner_states = self.states.as_slice();
        inner_states.split_off_first();
        inner_states.split_off_last();
        if inner_states.iter().any(|ciglet| ciglet.op == b'S' && ciglet.inc > 0) {
            return LocalTraverseFromAlignError::InternalClipping;
        }

        match PhmmState::from_op(op) {
            Some(PhmmState::Insert) => self.query_len_mismatch(),
            Some(PhmmState::Match | PhmmState::Delete) => self.model_len_mismatch(),
            None => LocalTraverseFromAlignError::InternalClipping,
        }
    }
}

impl CoreContextWithExitToErr for LocalContext<'_> {
    fn choose_end_insert_or_exit_op_error(&self, op: u8, _layer: DpIndex) -> Self::Error {
        // This function is only called when the operation is not I or S.
        // Confirm that the operation is in MD=X, or return invalid operation
        if PhmmState::from_op(op).is_none() {
            return self.invalid_op(op);
        }

        // Found one of MD=X, which consumes another reference position.
        // However, since traversal is at the end of the model, this is a model
        // length mismatch
        self.model_len_mismatch()
    }

    fn no_match_after_enter(&self, op: u8) -> Self::Error {
        LocalTraverseFromAlignError::MissingMatchOp {
            ref_len: self.ref_len,
            ref_coord: self.ref_start,
            op,
        }
    }

    fn ref_start_out_of_bounds(&self) -> Self::Error {
        LocalTraverseFromAlignError::RefStartOutOfBounds {
            ref_start: self.ref_start,
            ref_len:   self.ref_len,
        }
    }
}
