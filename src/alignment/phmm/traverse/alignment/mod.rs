//! Visitors for traversing pHMMs based on an alignment.

use crate::{
    alignment::{
        AlignmentStates, CigletOpIter, NextCiglet, PeekOp,
        phmm::{
            DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
            components::{EmissionParams, TransitionParams},
            indexing::{Begin, DpIndex, End, FirstMatch, GetLayer, GetModule, PhmmIndex, PhmmIndexable, SeqIndex},
            modules::{DomainModule, SemiLocalModule, SemiLocalParams},
            state::{PhmmState, PhmmStateOrModule},
            traverse::{
                DomainVisitor, EndInsert, EndInsertExit, GetScoreDomain, GetScoreLocal, GlobalVisitor, LocalVisitor,
                ModuleLocation, SemiLocalVisitor,
            },
        },
    },
    data::{ByteIndexMap, cigar::Ciglet},
};
use std::ops::{Range, RangeInclusive};

mod error;

pub use error::*;

/// A helper struct for all alignment visitors providing the necessary state and
/// functions for handling transitions and emissions in the core pHMM.
///
/// The score of the alignment is tracked as well, which is needed to resolve
/// ambiguities in domain, semilocal, and local pHMMs (by choosing the best
/// path).
///
/// ## Parameters
///
/// - `'a`: The lifetime of the alignment and the query.
/// - `T`: The type of the parameter used by the pHMM.
struct CoreAlignmentVisitor<'a, T, C> {
    /// The current score observed so far.
    ///
    /// This is needed to resolve ambiguities in domain, semilocal, and local
    /// pHMMs (by choosing the best path). All methods on
    /// [`CoreAlignmentVisitor`] update the score immediately, which agrees with
    /// the order of the floating point operations used in the Viterbi
    /// algorithm. Wrappers around this struct will read and mutate this field
    /// as well, but they may need to group parameters together before adding
    /// them in order to use the correct order and avoid floating point error.
    score:     T,
    /// An iterator over the operations in the CIGAR string. This is advanced
    /// when transitioning between states.
    op_iter:   CigletOpIter<'a>,
    /// The index in the query that is to be emitted next.
    query_idx: usize,
    /// The aligned-against query sequence corresponding to the core pHMM.
    query:     &'a [u8],
    /// Any context needed to form error messages of the proper type.
    context:   C,
}

impl<'a, T, C> CoreAlignmentVisitor<'a, T, C>
where
    T: PhmmNumber,
    C: CoreContextToErr,
{
    /// Initializes a new [`CoreAlignmentVisitor`] with the provided
    /// information.
    #[must_use]
    fn new(query: &'a [u8], states: &'a [Ciglet], context: C) -> Self {
        Self {
            score: T::ZERO,
            op_iter: CigletOpIter::new(states),
            query_idx: 0,
            query,
            context,
        }
    }

    /// The shared implementation for [`GlobalVisitor::choose_emission`],
    /// [`SemiLocalVisitor::choose_emission`],
    /// [`DomainVisitor::choose_emission`], and
    /// [`LocalVisitor::choose_emission`].
    ///
    /// ## Errors
    ///
    /// [`query_len_mismatch`] is returned if there are no remaining residues in
    /// the query.
    ///
    /// [`query_len_mismatch`]: CoreContextToErr::query_len_mismatch
    fn choose_emission<const S: usize>(
        &mut self, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
    ) -> Result<usize, C::Error> {
        match self.query.get(self.query_idx) {
            Some(byte) => {
                self.query_idx += 1;
                let idx = map.to_index(*byte);
                self.score += params[idx];
                Ok(idx)
            }
            None => Err(self.context.query_len_mismatch()),
        }
    }

    /// The shared implementation for choosing a transition within the core
    /// pHMM.
    ///
    /// ## Errors
    ///
    /// - [`core_transition_no_op_error`] if there is not another operation to
    ///   consume
    /// - [`core_transition_op_error`] if the operation is not in `MDI=X`
    ///
    /// [`core_transition_no_op_error`]: CoreContextToErr::core_transition_no_op_error
    /// [`core_transition_op_error`]: CoreContextToErr::core_transition_op_error
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<PhmmState, C::Error> {
        let op = self
            .op_iter
            .next()
            .ok_or_else(|| self.context.core_transition_no_op_error(layer, exiting))?;

        let next_state = PhmmState::from_op(op).ok_or_else(|| self.context.core_transition_op_error(op, layer, exiting))?;
        self.score += params[(exiting, next_state)];
        Ok(next_state)
    }

    /// The shared implementation for choosing the end or insert state at the
    /// end of the core pHMM.
    ///
    /// ## Errors
    ///
    /// [`choose_end_or_insert_op_error`] is returned if another operation is
    /// present but is not `I`.
    ///
    /// [`choose_end_or_insert_op_error`]: CoreContextToErr::choose_end_or_insert_op_error
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<EndInsert, C::Error> {
        let Some(op) = self.op_iter.next() else {
            self.score += params[(exiting, PhmmState::Match)];
            return Ok(EndInsert::End);
        };

        if PhmmState::from_op(op) != Some(PhmmState::Insert) {
            return Err(self.context.choose_end_or_insert_op_error(op, layer, exiting));
        }

        self.score += params[(exiting, PhmmState::Insert)];
        Ok(EndInsert::Insert)
    }
}

/// A wrapper around [`CoreAlignmentVisitor`] providing the necessary state and
/// functions for handling models which allow arbitrary entry and exit into the
/// model's match states (local and semilocal).
///
/// [`CoreAlignmentVisitor::score`] is updated only for transitions/emissions
/// within the core pHMM, and not for semilocal transitions into the core pHMM
/// or out of it. This is to allow the wrapper visitor to properly group the
/// parameters before adding them.
///
/// ## Parameters
///
/// - '`a`: The lifetime of the alignment and the query.
/// - `T`: The type of the parameter used by the pHMM.
/// - `ALLOW_CLIPPING`: If `true`, clipping is allowed (in the case of a local
///   pHMM). If `false`, clipping is not allowed (in the case of a semilocal
///   pHMM).
struct AlignmentVisitorWithExit<'a, T, C, const ALLOW_CLIPPING: bool> {
    /// The state for handling the transition out of states besides the match
    /// state, as well as emissions.
    inner:     CoreAlignmentVisitor<'a, T, C>,
    /// The starting 0-based index of the alignment within reference
    /// coordinates, used to infer when the pHMM is entered.
    ///
    /// If the alignment is empty (i.e., the query sequence is empty and as such
    /// there are no states in the alignment), then `ref_start` is ignored.
    /// Otherwise, it is validated that this is strictly in `0..phmm.seq_len()`.
    ref_start: usize,
}

impl<'a, T, C, const ALLOW_CLIPPING: bool> AlignmentVisitorWithExit<'a, T, C, ALLOW_CLIPPING>
where
    T: PhmmNumber,
    C: CoreContextWithExitToErr,
{
    /// Initializes a new [`AlignmentVisitorWithExit`] with the provided
    /// information.
    ///
    /// The starting 0-based index of the alignment within reference coordinates
    /// is passed. If the alignment is empty (i.e., the query sequence is empty
    /// and as such there are no states in the alignment), then `ref_start` is
    /// ignored.
    ///
    /// ## Errors
    ///
    /// `ref_start` must be less than the reference coordinate length to which
    /// the pHMM corresponds if the alignment is non-empty.
    fn new<P>(query: &'a [u8], states: &'a [Ciglet], ref_start: usize, phmm: &P, context: C) -> Result<Self, C::Error>
    where
        P: PhmmIndexable, {
        if !states.is_empty() && ref_start >= phmm.seq_len() {
            return Err(context.ref_start_out_of_bounds());
        }

        Ok(Self {
            inner: CoreAlignmentVisitor::new(query, states, context),
            ref_start,
        })
    }

    /// The shared implementation for choosing an emission within the core pHMM.
    ///
    /// ## Errors
    ///
    /// [`query_len_mismatch`] is returned if there are no remaining residues in
    /// the query.
    ///
    /// [`query_len_mismatch`]: CoreContextToErr::query_len_mismatch
    fn choose_emission<const S: usize>(
        &mut self, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
    ) -> Result<usize, C::Error> {
        self.inner.choose_emission(params, map)
    }

    /// The shared implementation for choosing a transition within the core
    /// pHMM.
    ///
    /// ## Errors
    ///
    /// - [`core_transition_no_op_error`] if there is not another operation to
    ///   consume
    /// - [`core_transition_op_error`] if the operation is not in `MDI=X`
    ///
    /// [`core_transition_no_op_error`]: CoreContextToErr::core_transition_no_op_error
    /// [`core_transition_op_error`]: CoreContextToErr::core_transition_op_error
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<PhmmState, C::Error> {
        self.inner.choose_core_transition(layer, exiting, params)
    }

    /// The shared implementation for choosing the end or insert state at the
    /// end of the core pHMM.
    ///
    /// ## Errors
    ///
    /// [`choose_end_or_insert_op_error`] is returned if another operation is
    /// present but is not `I`.
    ///
    /// [`choose_end_or_insert_op_error`]: CoreContextToErr::choose_end_or_insert_op_error
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<EndInsert, C::Error> {
        self.inner.choose_end_or_insert(layer, exiting, params)
    }

    /// The shared implementation for choosing an transition within the core
    /// pHMM or exiting from it.
    ///
    /// ## Errors
    ///
    /// [`invalid_op`] is returned if the operation is not in `MDI=X` (or `S`,
    /// if `ALLOW_CLIPPING` is true).
    ///
    /// [`invalid_op`]: CoreContextToErr::invalid_op
    fn choose_core_transition_or_exit(
        &mut self, params: &TransitionParams<T>, _exit_param: T,
    ) -> Result<PhmmStateOrModule, C::Error> {
        let Some(op) = self.inner.op_iter.peek_op() else {
            return Ok(PhmmStateOrModule::Module);
        };

        let next_state = if ALLOW_CLIPPING {
            PhmmStateOrModule::from_op(op).ok_or_else(|| self.inner.context.invalid_op(op))?
        } else {
            PhmmState::from_op(op)
                .ok_or_else(|| self.inner.context.invalid_op(op))?
                .into()
        };

        if let Some(next_state) = PhmmState::get_from(next_state) {
            self.inner.op_iter.next();
            self.inner.score += params[(PhmmState::Match, next_state)];
        }

        Ok(next_state)
    }

    /// A helper function for `choose_end_insert_or_exit`, consisting of the
    /// shared implementation for deciding whether the enter the insert state at
    /// the end of the core pHMM, specifically from the last match state. If
    /// `true`, then the insert state is entered. If `false`, then traversal may
    /// enter the END state or may exit directly. This ambiguity is handled by
    /// the wrapping visitor.
    ///
    /// If `ALLOW_CLIPPING` is false, then the only way this can return
    /// `Ok(false)` is if all operations from the alignment have been consumed.
    ///
    /// ## Errors
    ///
    /// - [`invalid_op`] if the next operation is `S` and `ALLOW_CLIPPING` is
    ///   false.
    /// - [`choose_end_insert_or_exit_op_error`] if the operation is not `I` or
    ///   `S`
    ///
    /// [`invalid_op`]: CoreContextToErr::invalid_op
    /// [`choose_end_insert_or_exit_op_error`]:
    ///     CoreContextWithExitToErr::choose_end_insert_or_exit_op_error
    fn enter_insert_at_end(&mut self, layer: DpIndex, params: &TransitionParams<T>) -> Result<bool, C::Error> {
        let Some(op) = self.inner.op_iter.peek_op() else {
            return Ok(false);
        };

        match PhmmStateOrModule::from_op(op) {
            Some(PhmmStateOrModule::Insert) => {
                self.inner.op_iter.next();
                self.inner.score += params[(PhmmState::Match, PhmmState::Insert)];
                Ok(true)
            }
            Some(PhmmStateOrModule::Module) if ALLOW_CLIPPING => Ok(false),
            Some(PhmmStateOrModule::Module) => Err(self.inner.context.invalid_op(op)),
            Some(PhmmStateOrModule::Match | PhmmStateOrModule::Delete) | None => {
                Err(self.inner.context.choose_end_insert_or_exit_op_error(op, layer))
            }
        }
    }

    /// The shared implementation for entering the core pHMM.
    ///
    /// If the reference range starts at 0, then this poses ambiguity with
    /// whether the begin state was passed through.
    ///
    /// If no states in the pHMM are matched against, this also poses ambiguity.
    /// Either the alignment enters the begin state and then immediately exits,
    /// or it enters the end state and then immediately exits.
    ///
    /// Both ambiguities are resolved by picking the path with the least score.
    ///
    /// ## Errors
    ///
    /// - [`invalid_op`] if the operation is not in `MDIS=X`.
    /// - [`no_match_after_enter`] if the start of the alignment within the pHMM
    ///   is not the BEGIN or END layer (so a match state is entered), and the
    ///   operation is not in `M=X`.
    ///
    /// [`invalid_op`]: CoreContextToErr::invalid_op
    /// [`no_match_after_enter`]: CoreContextWithExitToErr::no_match_after_enter
    fn enter_core<P, const S: usize>(&mut self, module: &SemiLocalModule<T>, phmm: &P) -> Result<DpIndex, C::Error>
    where
        P: GetModule<End: SemiLocalParams<T>> + PhmmIndexable + GetLayer<T, S>, {
        let next_op = self.inner.op_iter.peek_op();

        let next_op = match next_op {
            // Handle an empty alignment, which will ignore ref_start
            None | Some(b'S') => {
                let score_through_begin =
                    self.inner.score + module.get_score(Begin) + phmm.end().semilocal_params().get_score(Begin);
                let score_through_end =
                    self.inner.score + module.get_score(End) + phmm.end().semilocal_params().get_score(End);

                let idx = if score_through_begin <= score_through_end {
                    Begin.to_dp_index()
                } else {
                    End.to_dp_index(phmm)
                };

                return Ok(idx);
            }
            Some(op) => op,
        };

        // Already handled S above, so PhmmState is fine
        let Some(state) = PhmmState::from_op(next_op) else {
            return Err(self.inner.context.invalid_op(next_op));
        };

        // Handle case with no ambiguity (entering after first match state).
        if self.ref_start > 0 {
            // If we are not entering End, we are entering match state with
            // emission, so we must consume an `M`
            if state != PhmmState::Match {
                return Err(self.inner.context.no_match_after_enter(next_op));
            }

            // Advance op_iter instead of just peeking
            self.inner.op_iter.next();

            return Ok(SeqIndex(self.ref_start).to_dp_index());
        }

        match state {
            PhmmState::Match => {
                let begin_to_first_match_param = phmm.layer(Begin).transition[(PhmmState::Match, PhmmState::Match)];

                let skip_begin_score = self.inner.score + module.get_score(FirstMatch);
                let through_begin_score = self.inner.score + module.get_score(Begin) + begin_to_first_match_param;

                let layer = if through_begin_score <= skip_begin_score {
                    Begin.to_dp_index()
                } else {
                    self.inner.op_iter.next();
                    FirstMatch.to_dp_index()
                };

                Ok(layer)
            }
            PhmmState::Insert | PhmmState::Delete => Ok(Begin.to_dp_index()),
        }
    }
}

/// A global pHMM visitor, following the path given by an alignment.
pub struct GlobalAlignmentVisitor<'a, T> {
    inner: CoreAlignmentVisitor<'a, T, GlobalContext<'a>>,
}

impl<'a, T> GlobalAlignmentVisitor<'a, T>
where
    T: PhmmNumber,
{
    /// Initializes a new [`GlobalAlignmentVisitor`] from the given query and
    /// [`AlignmentStates`].
    #[must_use]
    pub fn new<const S: usize>(query: &'a [u8], states: &'a AlignmentStates, phmm: &GlobalPhmm<T, S>) -> Self {
        let context = GlobalContext {
            query_len: query.len(),
            states,
            ref_len: phmm.seq_len(),
        };

        Self {
            inner: CoreAlignmentVisitor::new(query, states.as_slice(), context),
        }
    }
}

impl<T, const S: usize> GlobalVisitor<T, S> for GlobalAlignmentVisitor<'_, T>
where
    T: PhmmNumber,
{
    type Output = T;
    type Error = GlobalTraverseFromAlignError;

    /// Selects the index of the residue that is emitted at a match or insert
    /// state in the pHMM.
    ///
    /// ## Errors
    ///
    /// [`QueryLenMismatch`] is returned if there are no remaining residues in
    /// the query.
    ///
    /// [`QueryLenMismatch`]: GlobalTraverseFromAlignError::QueryLenMismatch
    fn choose_emission(
        &mut self, _layer: DpIndex, _state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        _phmm: &GlobalPhmm<T, S>,
    ) -> Result<usize, GlobalTraverseFromAlignError> {
        self.inner.choose_emission(params, map)
    }

    /// Selects the transition within the pHMM to take (what the next
    /// [`PhmmState`] to enter should be).
    ///
    /// When in the [`LastMatch`] layer, [`choose_end_or_insert`] is called
    /// instead.
    ///
    /// ## Errors
    ///
    /// If there is not another operation to consume, then [`ModelLenMismatch`]
    /// is returned. Otherwise, if the operation is not in `MDI=X`, then one of
    /// [`InvalidCigarOp`], [`QueryLenMismatch`], or [`ModelLenMismatch`] is
    /// returned.
    ///
    /// [`InvalidCigarOp`]: GlobalTraverseFromAlignError::InvalidCigarOp
    /// [`QueryLenMismatch`]: GlobalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: GlobalTraverseFromAlignError::ModelLenMismatch
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_end_or_insert`]: GlobalVisitor::choose_end_or_insert
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, _phmm: &GlobalPhmm<T, S>,
    ) -> Result<PhmmState, GlobalTraverseFromAlignError> {
        self.inner.choose_core_transition(layer, exiting, params)
    }

    /// Selects whether the traversal should enter the [`End`] state from any of
    /// the states in the [`LastMatch`] layer, or whether the final insert state
    /// should be entered.
    ///
    /// ## Errors
    ///
    /// - [`ModelLenMismatch`] if `M` or `D` is the next operation in the
    ///   alignment
    /// - [`InvalidCigarOp`] if the operation is not in `MDI=X`
    ///
    /// [`ModelLenMismatch`]: GlobalTraverseFromAlignError::ModelLenMismatch
    /// [`InvalidCigarOp`]: GlobalTraverseFromAlignError::InvalidCigarOp
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, _phmm: &GlobalPhmm<T, S>,
    ) -> Result<EndInsert, GlobalTraverseFromAlignError> {
        self.inner.choose_end_or_insert(layer, exiting, params)
    }

    /// Returns the final score, as well as checking that the traversal ended in
    /// a valid state.
    ///
    /// ## Errors
    ///
    /// Returns [`QueryLenMismatch`] if the full query was not consumed.
    ///
    /// [`QueryLenMismatch`]: GlobalTraverseFromAlignError::QueryLenMismatch
    fn finalize(mut self, _phmm: &GlobalPhmm<T, S>) -> Result<T, GlobalTraverseFromAlignError> {
        // It is known that op_iter will be empty since choose_end_or_insert
        // will continue to return Insert until either op_iter is empty or an
        // error is thrown
        debug_assert!(self.inner.op_iter.next().is_none());

        if self.inner.query_idx < self.inner.query.len() {
            Err(self.inner.context.query_len_mismatch())
        } else {
            Ok(self.inner.score)
        }
    }
}

/// A semilocal pHMM visitor, following the path given by an alignment.
///
/// Any ambiguities in the path are resolved by picking the one with the least
/// score.
pub struct SemiLocalAlignmentVisitor<'a, T> {
    inner: AlignmentVisitorWithExit<'a, T, SemiLocalContext<'a>, false>,
}

impl<'a, T> SemiLocalAlignmentVisitor<'a, T>
where
    T: PhmmNumber + 'static,
{
    /// Initializes a new [`SemiLocalAlignmentVisitor`] from the given query and
    /// alignment.
    ///
    /// The starting 0-based index of the alignment within reference coordinates
    /// is passed. If the alignment is empty (i.e., the query sequence is empty
    /// and as such there are no states in the alignment), then `ref_start` is
    /// ignored.
    ///
    /// ## Errors
    ///
    /// [`RefStartOutOfBounds`] is returned if `ref_start` is not less than the
    /// reference coordinate length to which the pHMM corresponds.
    ///
    /// [`RefStartOutOfBounds`]:
    ///     SemiLocalTraverseFromAlignError::RefStartOutOfBounds
    pub fn new<const S: usize>(
        query: &'a [u8], states: &'a AlignmentStates, ref_start: usize, phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<Self, SemiLocalTraverseFromAlignError> {
        let context = SemiLocalContext {
            query_len: query.len(),
            states,
            ref_len: phmm.seq_len(),
            ref_start,
        };

        Ok(Self {
            inner: AlignmentVisitorWithExit::new(query, states.as_slice(), ref_start, phmm, context)?,
        })
    }
}

impl<T, const S: usize> SemiLocalVisitor<T, S> for SemiLocalAlignmentVisitor<'_, T>
where
    T: PhmmNumber,
{
    type Output = T;
    type Error = SemiLocalTraverseFromAlignError;

    /// Selects the index of the residue that is emitted at a match or insert
    /// state in the pHMM.
    ///
    /// ## Errors
    ///
    /// [`QueryLenMismatch`] is returned if there are no remaining residues in
    /// the query.
    ///
    /// [`QueryLenMismatch`]: SemiLocalTraverseFromAlignError::QueryLenMismatch
    fn choose_emission(
        &mut self, _layer: DpIndex, _state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        _phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<usize, SemiLocalTraverseFromAlignError> {
        self.inner.choose_emission(params, map)
    }

    /// Selects the transition within the core pHMM to take (what the next
    /// [`PhmmState`] to enter should be) when in an insert or delete state.
    ///
    /// When in a match state, [`choose_core_transition_or_exit`] is called
    /// instead. When in the [`LastMatch`] layer, [`choose_end_or_insert`] is
    /// called instead.
    ///
    /// ## Errors
    ///
    /// - [`InvalidEarlyExit`] if there is not another operation to consume
    /// - [`InvalidCigarOp`] if the operation is not in `MDI=X`
    ///
    /// [`InvalidEarlyExit`]: SemiLocalTraverseFromAlignError::InvalidEarlyExit
    /// [`InvalidCigarOp`]: SemiLocalTraverseFromAlignError::InvalidCigarOp
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_core_transition_or_exit`]:
    ///     SemiLocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_or_insert`]: SemiLocalVisitor::choose_end_or_insert
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, _phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<PhmmState, SemiLocalTraverseFromAlignError> {
        self.inner.choose_core_transition(layer, exiting, params)
    }

    /// Selects whether the traversal should enter the [`End`] state from the
    /// insert or delete states in the [`LastMatch`] layer, or whether the final
    /// insert state should be entered.
    ///
    /// When in the match state of the [`LastMatch`] layer,
    /// [`choose_end_insert_or_exit`] is called instead.
    ///
    /// ## Errors
    ///
    /// - [`ModelLenMismatch`] if `M` or `D` is the next operation in the
    ///   alignment
    /// - [`InvalidCigarOp`] if the operation is not in `MDI=X`
    ///
    /// [`ModelLenMismatch`]: SemiLocalTraverseFromAlignError::ModelLenMismatch
    /// [`InvalidCigarOp`]: SemiLocalTraverseFromAlignError::InvalidCigarOp
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_end_insert_or_exit`]:
    ///     SemiLocalVisitor::choose_end_insert_or_exit
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, _phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<EndInsert, SemiLocalTraverseFromAlignError> {
        self.inner.choose_end_or_insert(layer, exiting, params)
    }

    /// Selects the transition within the core pHMM to take (what the next
    /// [`PhmmState`] to enter should be), or whether to exit from the pHMM
    /// early.
    ///
    /// This is only called when in a match state.
    ///
    /// The `score` field of the inner visitor is updated with the transition,
    /// including transitions out of the core pHMM.
    ///
    /// ## Errors
    ///
    /// [`InvalidCigarOp`] is returned if the operation is not in `MDI=X`
    ///
    /// [`InvalidCigarOp`]: SemiLocalTraverseFromAlignError::InvalidCigarOp
    fn choose_core_transition_or_exit(
        &mut self, _layer: DpIndex, params: &TransitionParams<T>, exit_param: T, _phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<PhmmStateOrModule, SemiLocalTraverseFromAlignError> {
        self.inner.choose_core_transition_or_exit(params, exit_param)
    }

    /// From the match state in the [`LastMatch`] layer, selects whether the
    /// traversal should enter the [`End`] state, enter the final insert state,
    /// or exit early from the pHMM.
    ///
    /// The `score` field of the inner visitor is updated with the transition,
    /// including transitions out of the core pHMM.
    ///
    /// ## Errors
    ///
    /// If the next operation is present but not equal to `I` or `S`, then one
    /// of [`InvalidCigarOp`], [`QueryLenMismatch`], or [`ModelLenMismatch`] is
    /// returned, depending on the cause.
    ///
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`InvalidCigarOp`]: SemiLocalTraverseFromAlignError::InvalidCigarOp
    /// [`QueryLenMismatch`]: SemiLocalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: SemiLocalTraverseFromAlignError::ModelLenMismatch
    fn choose_end_insert_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
        _phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<EndInsertExit, SemiLocalTraverseFromAlignError> {
        let enter_insert = self.inner.enter_insert_at_end(layer, params)?;

        if enter_insert {
            return Ok(EndInsertExit::Insert);
        }

        let score_with_exit = self.inner.inner.score + exit_param;
        let score_through_end = self.inner.inner.score + params[(PhmmState::Match, PhmmState::Match)] + exit_from_end_param;

        if score_through_end <= score_with_exit {
            self.inner.inner.score += params[(PhmmState::Match, PhmmState::Match)];
            Ok(EndInsertExit::End)
        } else {
            Ok(EndInsertExit::Exit)
        }
    }

    /// Selects the layer of the pHMM to enter from the [`SemiLocalModule`] at
    /// the start of the pHMM.
    ///
    /// ## Errors
    ///
    /// [`InvalidCigarOp`], [`QueryLenMismatch`], or [`ModelLenMismatch`] is
    /// returned if the next operation is not in `MDIS=X`. [`MissingMatchOp`] is
    /// returned if the start of the alignment within the pHMM is not the BEGIN
    /// or END layer (so a match state is entered), and the operation is not in
    /// `M=X`.
    ///
    /// [`InvalidCigarOp`]: SemiLocalTraverseFromAlignError::InvalidCigarOp
    /// [`QueryLenMismatch`]: SemiLocalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: SemiLocalTraverseFromAlignError::ModelLenMismatch
    /// [`MissingMatchOp`]: SemiLocalTraverseFromAlignError::MissingMatchOp
    fn enter_core(
        &mut self, module: &SemiLocalModule<T>, phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<DpIndex, SemiLocalTraverseFromAlignError> {
        let index = self.inner.enter_core(module, phmm)?;

        // Update the score for transitions into the core pHMM
        self.inner.inner.score += module.get_score(index);

        Ok(index)
    }

    /// Performs any behavior necessary given that the traversal is exiting the
    /// END state into the [`SemiLocalModule`] at the end of the pHMM.
    ///
    /// ## Errors
    ///
    /// This implementation is infallible.
    fn exit_core_from_end(
        &mut self, _layer: DpIndex, _exit_param: T, _phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<(), SemiLocalTraverseFromAlignError> {
        Ok(())
    }

    /// Performs any actions upon exiting the core pHMM.
    ///
    /// This is provided as a hook for actions that must be performed upon
    /// exiting the core pHMM. It will be called directly after
    /// [`choose_core_transition_or_exit`] returns
    /// [`PhmmStateOrModule::Module`], [`choose_end_insert_or_exit`] returns
    /// [`EndInsertExit::Exit`], or [`exit_core_from_end`] is called.
    ///
    /// ## Errors
    ///
    /// This implementation is infallible.
    ///
    /// [`choose_core_transition_or_exit`]:
    ///     SemiLocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_insert_or_exit`]:
    ///     SemiLocalVisitor::choose_end_insert_or_exit
    /// [`exit_core_from_end`]: SemiLocalVisitor::exit_core_from_end
    fn exit_core(
        &mut self, _layer_idx: DpIndex, exit_param: T, _phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<(), SemiLocalTraverseFromAlignError> {
        self.inner.inner.score += exit_param;
        Ok(())
    }

    /// Returns the final score, as well as checking that the traversal ended in
    /// a valid state.
    ///
    /// ## Errors
    ///
    /// If an additional operation is present, then one of [`InvalidCigarOp`],
    /// [`QueryLenMismatch`], or [`ModelLenMismatch`] is returned depending on
    /// the cause. If the full query wasn't consumed, then [`QueryLenMismatch`]
    /// is returned.
    ///
    /// [`InvalidCigarOp`]: SemiLocalTraverseFromAlignError::InvalidCigarOp
    /// [`QueryLenMismatch`]: SemiLocalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: SemiLocalTraverseFromAlignError::ModelLenMismatch
    fn finalize(
        mut self, _phmm: &SemiLocalPhmm<T, S>, _aligned_layers: RangeInclusive<DpIndex>,
    ) -> Result<T, SemiLocalTraverseFromAlignError> {
        // op_iter will normally be empty based on enter_insert_at_end's
        // guarantees, but if enter_core directly enters the END state, then we
        // could have problems
        if let Some(op) = self.inner.inner.op_iter.next() {
            // Validate the operation
            if PhmmState::from_op(op).is_none() {
                return Err(self.inner.inner.context.invalid_op(op));
            }

            return Err(SemiLocalTraverseFromAlignError::OpAfterEnteringEnd { op });
        }

        if self.inner.inner.query_idx != self.inner.inner.query.len() {
            return Err(self.inner.inner.context.query_len_mismatch());
        }

        Ok(self.inner.inner.score)
    }
}

/// A domain pHMM visitor, following the path given by an alignment.
pub struct DomainAlignmentVisitor<'a, T> {
    inner:          CoreAlignmentVisitor<'a, T, DomainContext<'a>>,
    begin_residues: &'a [u8],
    end_residues:   &'a [u8],
}

impl<'a, T> DomainAlignmentVisitor<'a, T>
where
    T: PhmmNumber + 'static,
{
    /// Initializes a new [`DomainAlignmentVisitor`] from the given query and
    /// [`AlignmentStates`].
    ///
    /// ## Errors
    ///
    /// Returns [`QueryLenMismatch`] if the query is not long enough based on
    /// the soft clipping present in `states`.
    ///
    /// [`QueryLenMismatch`]: DomainTraverseFromAlignError::QueryLenMismatch
    pub fn new<const S: usize>(
        mut query: &'a [u8], states: &'a AlignmentStates, phmm: &DomainPhmm<T, S>,
    ) -> Result<Self, DomainTraverseFromAlignError> {
        // skipped_start and skipped_end get initialized later
        let mut context = DomainContext {
            query_len: query.len(),
            states,
            ref_len: phmm.seq_len(),
            skipped_start: 0,
            skipped_end: 0,
        };

        let mut states = states.as_slice();

        let begin_residues_len = states.next_ciglet_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
        let end_residues_len = states.next_ciglet_back_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);

        // Finish initializing context
        context.skipped_start = begin_residues_len;
        context.skipped_end = end_residues_len;

        let (begin_residues, end_residues) = if let Some(end_start) = query.len().checked_sub(end_residues_len)
            && let Some(end_residues) = query.split_off(end_start..)
            && let Some(begin_residues) = query.split_off(..begin_residues_len)
        {
            (begin_residues, end_residues)
        } else {
            return Err(context.query_len_mismatch());
        };

        Ok(Self {
            inner: CoreAlignmentVisitor::new(query, states, context),
            begin_residues,
            end_residues,
        })
    }
}

impl<T, const S: usize> DomainVisitor<T, S> for DomainAlignmentVisitor<'_, T>
where
    T: PhmmNumber,
{
    type Output = T;
    type Error = DomainTraverseFromAlignError;

    /// Selects the index of the residue that is emitted at a match or insert
    /// state in the core pHMM.
    ///
    /// ## Errors
    ///
    /// [`QueryLenMismatch`] is returned if there are no remaining residues in
    /// the query.
    ///
    /// [`QueryLenMismatch`]: DomainTraverseFromAlignError::QueryLenMismatch
    fn choose_emission(
        &mut self, _layer: DpIndex, _state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        _phmm: &DomainPhmm<T, S>,
    ) -> Result<usize, DomainTraverseFromAlignError> {
        self.inner.choose_emission(params, map)
    }

    /// Selects the transition within the pHMM to take (what the next
    /// [`PhmmState`] to enter should be).
    ///
    /// When in the [`LastMatch`] layer, [`choose_end_or_insert`] is called
    /// instead.
    ///
    /// ## Errors
    ///
    /// If there is not another operation to consume, then [`ModelLenMismatch`]
    /// is returned. Otherwise, if the operation is not in `MDI=X`, then one of
    /// [`InvalidCigarOp`], [`DuplicateOp`], [`InternalClipping`],
    /// [`QueryLenMismatch`], or [`ModelLenMismatch`] is returned.
    ///
    /// [`ModelLenMismatch`]: DomainTraverseFromAlignError::ModelLenMismatch
    /// [`InvalidCigarOp`]: DomainTraverseFromAlignError::InvalidCigarOp
    /// [`DuplicateOp`]: DomainTraverseFromAlignError::DuplicateOp
    /// [`InternalClipping`]: DomainTraverseFromAlignError::InternalClipping
    /// [`QueryLenMismatch`]: DomainTraverseFromAlignError::QueryLenMismatch
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_end_or_insert`]: DomainVisitor::choose_end_or_insert
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, _phmm: &DomainPhmm<T, S>,
    ) -> Result<PhmmState, DomainTraverseFromAlignError> {
        self.inner.choose_core_transition(layer, exiting, params)
    }

    /// Selects whether the traversal should enter the [`End`] state from any of
    /// the states in the [`LastMatch`] layer, or whether the final insert state
    /// should be entered.
    ///
    /// ## Errors
    ///
    /// If another operation is present but is not `I`, one of
    /// [`InvalidCigarOp`], [`DuplicateOp`], [`InternalClipping`],
    /// [`QueryLenMismatch`], or [`ModelLenMismatch`] is returned.
    ///
    /// [`InvalidCigarOp`]: DomainTraverseFromAlignError::InvalidCigarOp
    /// [`DuplicateOp`]: DomainTraverseFromAlignError::DuplicateOp
    /// [`InternalClipping`]: DomainTraverseFromAlignError::InternalClipping
    /// [`QueryLenMismatch`]: DomainTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: DomainTraverseFromAlignError::ModelLenMismatch
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, _phmm: &DomainPhmm<T, S>,
    ) -> Result<EndInsert, DomainTraverseFromAlignError> {
        if self.inner.op_iter.peek_op() == Some(b'S') {
            self.inner.score += params[(exiting, PhmmState::Match)];
            return Ok(EndInsert::End);
        }

        self.inner.choose_end_or_insert(layer, exiting, params)
    }

    /// Selects the index of the residue that is emitted within either
    /// [`DomainModule`].
    ///
    /// ## Errors
    ///
    /// [`QueryLenMismatch`] is returned if there are no remaining residues in
    /// the query.
    ///
    /// [`QueryLenMismatch`]: DomainTraverseFromAlignError::QueryLenMismatch
    fn choose_domain_emission(
        &mut self, _params: &EmissionParams<T, S>, mapping: &ByteIndexMap<S>, loc: ModuleLocation, _phmm: &DomainPhmm<T, S>,
    ) -> Result<usize, DomainTraverseFromAlignError> {
        // This function does not update the score, since that was done
        // pre-emptively in enter_module_insert
        let byte = match loc {
            ModuleLocation::Begin => self.begin_residues.split_off_first(),
            ModuleLocation::End => self.end_residues.split_off_first(),
        }
        .ok_or_else(|| self.inner.context.query_len_mismatch())?;

        Ok(mapping.to_index(*byte))
    }

    /// Selects whether to enter the insert state within either
    /// [`DomainModule`].
    ///
    /// If `true` is returned, then the insert state is entered. If `false` is
    /// returned, traversal continues to the end of the [`DomainModule`].
    ///
    /// ## Errors
    ///
    /// This implementation is infallible.
    fn enter_module_insert(
        &mut self, _module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<bool, DomainTraverseFromAlignError> {
        // This function pre-emptively updates the score all at once for the
        // domain module. This is to ensure the order of floating point
        // operations is correct.
        match loc {
            ModuleLocation::Begin => {
                self.inner.score += phmm.get_begin_score(self.begin_residues);
                Ok(!self.begin_residues.is_empty())
            }
            ModuleLocation::End => {
                self.inner.score += phmm.get_end_score(self.end_residues);
                Ok(!self.end_residues.is_empty())
            }
        }
    }

    /// Selects whether to exit the insert state within either [`DomainModule`].
    ///
    /// If `true` is returned, then the insert state is exited and traversal
    /// continues to the end of the [`DomainModule`]. If `false` is returned,
    /// then traversal stays within the insert state.
    ///
    /// ## Errors
    ///
    /// This implementation is infallible.
    fn exit_module_insert(
        &mut self, _module: &DomainModule<T, S>, loc: ModuleLocation, _phmm: &DomainPhmm<T, S>,
    ) -> Result<bool, DomainTraverseFromAlignError> {
        // This function does not update the score, since that was done
        // pre-emptively in enter_module_insert
        match loc {
            ModuleLocation::Begin => Ok(self.begin_residues.is_empty()),
            ModuleLocation::End => Ok(self.end_residues.is_empty()),
        }
    }

    fn exiting_module(
        &mut self, _module: &DomainModule<T, S>, _loc: ModuleLocation, _phmm: &DomainPhmm<T, S>,
    ) -> Result<(), Self::Error> {
        // TODO: Likely want to do score updates here for clarity, rather than
        // in enter_module_insert
        Ok(())
    }

    /// Returns the final score, as well as checking that the traversal ended in
    /// a valid state.
    ///
    /// ## Errors
    ///
    /// If an additional operation is present, then one of [`InvalidCigarOp`],
    /// [`DuplicateOp`], [`InternalClipping`], [`QueryLenMismatch`], or
    /// [`ModelLenMismatch`] is returned depending on the cause. If the full
    /// query wasn't consumed, then [`QueryLenMismatch`] is returned.
    ///
    /// [`InvalidCigarOp`]: DomainTraverseFromAlignError::InvalidCigarOp
    /// [`DuplicateOp`]: DomainTraverseFromAlignError::DuplicateOp
    /// [`InternalClipping`]: DomainTraverseFromAlignError::InternalClipping
    /// [`QueryLenMismatch`]: DomainTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: DomainTraverseFromAlignError::ModelLenMismatch
    fn finalize(
        mut self, _phmm: &DomainPhmm<T, S>, _aligned_seq: Range<SeqIndex>,
    ) -> Result<T, DomainTraverseFromAlignError> {
        if let Some(op) = self.inner.op_iter.next() {
            return Err(self.inner.context.remaining_op_error(op));
        }

        if self.inner.query_idx != self.inner.query.len() {
            return Err(self.inner.context.query_len_mismatch());
        }

        Ok(self.inner.score)
    }
}

/// A local pHMM visitor, following the path given by an alignment.
///
/// Any ambiguities in the path are resolved by picking the one with the least
/// score.
pub struct LocalAlignmentVisitor<'a, T> {
    /// The inner visitor combining functionality for semilocal and local pHMMs.
    inner:          AlignmentVisitorWithExit<'a, T, LocalContext<'a>, true>,
    /// The residues remaining which must be emit from the module at the
    /// beginning of the pHMM.
    begin_residues: &'a [u8],
    /// The residues remaining which must be emit from the module at the end of
    /// the pHMM.
    end_residues:   &'a [u8],
    /// Any pre-computed information if the alignment is empty (all soft
    /// clipping).
    empty_info:     Option<EmptyInfo>,
}

struct EmptyInfo {
    through_begin: bool,
}

impl<'a, T> LocalAlignmentVisitor<'a, T>
where
    T: PhmmNumber,
{
    /// Initializes a new [`LocalAlignmentVisitor`] from the given query and
    /// alignment.
    ///
    /// The starting 0-based index of the alignment within reference coordinates
    /// is passed. If the alignment is empty (i.e., all soft-clipping), then
    /// `ref_start` is ignored.
    ///
    /// ## Errors
    ///
    /// [`RefStartOutOfBounds`] is returned if `ref_start` is not less than the
    /// reference coordinate length to which the pHMM corresponds.
    ///
    /// [`RefStartOutOfBounds`]:
    ///     LocalTraverseFromAlignError::RefStartOutOfBounds
    pub fn new<const S: usize>(
        mut query: &'a [u8], states: &'a AlignmentStates, ref_start: usize, phmm: &LocalPhmm<T, S>,
    ) -> Result<Self, LocalTraverseFromAlignError> {
        // skipped_start and skipped_end get initialized later
        let mut context = LocalContext {
            query_len: query.len(),
            states,
            ref_len: phmm.seq_len(),
            ref_start,
            skipped_start: 0,
            skipped_end: 0,
        };

        let mut states = states.as_slice();

        // Determine whether the alignment is entirely clipping
        let all_clipping_states = &[Ciglet {
            op:  b'S',
            inc: query.len(),
        }];

        let all_clipping = states == all_clipping_states || (query.is_empty() && states.is_empty());

        let (begin_residues, end_residues, empty_info) = if all_clipping {
            // In case all the possible scorings yield INFINITY, by default we
            // pick to put all emissions in the module at the beginning.
            let mut best_begin_residues: &[u8] = query;
            let mut best_end_residues: &[u8] = b"";
            let mut best_through_begin = true;
            let mut best_score = T::INFINITY;

            for begin_residues_len in 0..=query.len() {
                let (begin_residues, end_residues) = query.split_at(begin_residues_len);
                for (through_begin, through_state) in [(true, Begin.to_dp_index()), (false, End.to_dp_index(phmm))] {
                    let begin_score = phmm.get_begin_score(begin_residues, through_state);
                    let end_score = phmm.get_end_score(end_residues, through_state);
                    let score = begin_score + end_score;

                    if score < best_score {
                        best_begin_residues = begin_residues;
                        best_end_residues = end_residues;
                        best_through_begin = through_begin;
                        best_score = score;
                    }
                }
            }

            let empty_info = EmptyInfo {
                through_begin: best_through_begin,
            };

            // Clear states and query
            states = &[];
            query = &[];

            context.skipped_start = best_begin_residues.len();
            context.skipped_end = best_end_residues.len();

            (best_begin_residues, best_end_residues, Some(empty_info))
        } else {
            let begin_residues_len = states.next_ciglet_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
            let end_residues_len = states.next_ciglet_back_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);

            context.skipped_start = begin_residues_len;
            context.skipped_end = end_residues_len;

            let (begin_residues, end_residues) = if let Some(end_start) = query.len().checked_sub(end_residues_len)
                && let Some(end_residues) = query.split_off(end_start..)
                && let Some(begin_residues) = query.split_off(..begin_residues_len)
            {
                (begin_residues, end_residues)
            } else {
                return Err(context.query_len_mismatch());
            };

            (begin_residues, end_residues, None)
        };

        let inner = AlignmentVisitorWithExit::new(query, states, ref_start, phmm, context)?;

        Ok(Self {
            inner,
            begin_residues,
            end_residues,
            empty_info,
        })
    }

    fn compute_end_module_score<const S: usize>(&self, exit_param: T, phmm: &LocalPhmm<T, S>) -> T {
        let domain_score = phmm.end().domain_params.get_end_score(self.end_residues, phmm.mapping());
        domain_score + exit_param
    }
}

impl<T, const S: usize> LocalVisitor<T, S> for LocalAlignmentVisitor<'_, T>
where
    T: PhmmNumber,
{
    type Output = T;
    type Error = LocalTraverseFromAlignError;

    /// Selects the index of the residue that is emitted at a match or insert
    /// state in the pHMM.
    ///
    /// ## Errors
    ///
    /// [`QueryLenMismatch`] is returned if there are no remaining residues in
    /// the query.
    ///
    /// [`QueryLenMismatch`]: LocalTraverseFromAlignError::QueryLenMismatch
    fn choose_emission(
        &mut self, _layer: DpIndex, _state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        _phmm: &LocalPhmm<T, S>,
    ) -> Result<usize, LocalTraverseFromAlignError> {
        self.inner.choose_emission(params, map)
    }

    /// Selects the transition within the core pHMM to take (what the next
    /// [`PhmmState`] to enter should be) when in an insert or delete state.
    ///
    /// When in a match state, [`choose_core_transition_or_exit`] is called
    /// instead. When in the [`LastMatch`] layer, [`choose_end_or_insert`] is
    /// called instead.
    ///
    /// ## Errors
    ///
    /// If there is not another operation to consume, then [`InvalidEarlyExit`]
    /// is returned. Otherwise, if the operation is not in `MDI=X`, then one of
    /// [`InvalidCigarOp`], [`DuplicateOp`], [`InternalClipping`],
    /// [`QueryLenMismatch`], or [`ModelLenMismatch`] is returned.
    ///
    /// [`InvalidEarlyExit`]: LocalTraverseFromAlignError::InvalidEarlyExit
    /// [`InvalidCigarOp`]: LocalTraverseFromAlignError::InvalidCigarOp
    /// [`DuplicateOp`]: LocalTraverseFromAlignError::DuplicateOp
    /// [`InternalClipping`]: LocalTraverseFromAlignError::InternalClipping
    /// [`QueryLenMismatch`]: LocalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: LocalTraverseFromAlignError::ModelLenMismatch
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_core_transition_or_exit`]:
    ///     LocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_or_insert`]: LocalVisitor::choose_end_or_insert
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, _phmm: &LocalPhmm<T, S>,
    ) -> Result<PhmmState, LocalTraverseFromAlignError> {
        self.inner.choose_core_transition(layer, exiting, params)
    }

    /// Selects whether the traversal should enter the [`End`] state from the
    /// insert or delete states in the [`LastMatch`] layer, or whether the final
    /// insert state should be entered.
    ///
    /// When in the match state of the [`LastMatch`] layer,
    /// [`choose_end_insert_or_exit`] is called instead.
    ///
    /// ## Errors
    ///
    /// If another operation is present but is not `I`, then [`InvalidCigarOp`],
    /// [`DuplicateOp`], [`InternalClipping`], [`QueryLenMismatch`], or
    /// [`ModelLenMismatch`] is returned depending on the cause.
    ///
    /// [`InvalidCigarOp`]: LocalTraverseFromAlignError::InvalidCigarOp
    /// [`DuplicateOp`]: LocalTraverseFromAlignError::DuplicateOp
    /// [`InternalClipping`]: LocalTraverseFromAlignError::InternalClipping
    /// [`QueryLenMismatch`]: LocalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: LocalTraverseFromAlignError::ModelLenMismatch
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_end_insert_or_exit`]: LocalVisitor::choose_end_insert_or_exit
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, _phmm: &LocalPhmm<T, S>,
    ) -> Result<EndInsert, LocalTraverseFromAlignError> {
        if self.inner.inner.op_iter.peek_op() == Some(b'S') {
            self.inner.inner.score += params[(exiting, PhmmState::Match)];
            return Ok(EndInsert::End);
        }

        self.inner.choose_end_or_insert(layer, exiting, params)
    }

    /// Selects the transition within the core pHMM to take (what the next
    /// [`PhmmState`] to enter should be), or whether to exit from the pHMM
    /// early.
    ///
    /// This is only called when in a match state.
    ///
    /// ## Errors
    ///
    /// If the operation is not in `MDIS=X`, then one of [`InvalidCigarOp`],
    /// [`DuplicateOp`], [`InternalClipping`], [`QueryLenMismatch`], or
    /// [`ModelLenMismatch`] is returned, depending on the cause.
    ///
    /// [`InvalidCigarOp`]: LocalTraverseFromAlignError::InvalidCigarOp
    /// [`DuplicateOp`]: LocalTraverseFromAlignError::DuplicateOp
    /// [`InternalClipping`]: LocalTraverseFromAlignError::InternalClipping
    /// [`QueryLenMismatch`]: LocalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: LocalTraverseFromAlignError::ModelLenMismatch
    fn choose_core_transition_or_exit(
        &mut self, _layer: DpIndex, params: &TransitionParams<T>, exit_param: T, _phmm: &LocalPhmm<T, S>,
    ) -> Result<PhmmStateOrModule, LocalTraverseFromAlignError> {
        let out = if self.empty_info.is_some() {
            PhmmStateOrModule::Module
        } else {
            self.inner.choose_core_transition_or_exit(params, exit_param)?
        };

        Ok(out)
    }

    /// From the match state in the [`LastMatch`] layer, selects whether the
    /// traversal should enter the [`End`] state, enter the final insert state,
    /// or exit early from the pHMM.
    ///
    /// ## Errors
    ///
    /// If the next operation is present but not equal to `I` or `S`, then one
    /// of [`InvalidCigarOp`], [`DuplicateOp`], [`InternalClipping`],
    /// [`QueryLenMismatch`], or [`ModelLenMismatch`] is returned, depending on
    /// the cause.
    ///
    /// [`InvalidCigarOp`]: LocalTraverseFromAlignError::InvalidCigarOp
    /// [`DuplicateOp`]: LocalTraverseFromAlignError::DuplicateOp
    /// [`InternalClipping`]: LocalTraverseFromAlignError::InternalClipping
    /// [`QueryLenMismatch`]: LocalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: LocalTraverseFromAlignError::ModelLenMismatch
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`End`]: crate::alignment::phmm::indexing::End
    fn choose_end_insert_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
        phmm: &LocalPhmm<T, S>,
    ) -> Result<EndInsertExit, LocalTraverseFromAlignError> {
        if self.empty_info.is_some() {
            return Ok(EndInsertExit::Exit);
        }

        if self.inner.enter_insert_at_end(layer, params)? {
            return Ok(EndInsertExit::Insert);
        }

        let score_with_exit = self.inner.inner.score + self.compute_end_module_score(exit_param, phmm);
        let score_through_end = self.inner.inner.score
            + params[(PhmmState::Match, PhmmState::Match)]
            + self.compute_end_module_score(exit_from_end_param, phmm);

        if score_through_end <= score_with_exit {
            self.inner.inner.score += params[(PhmmState::Match, PhmmState::Match)];
            Ok(EndInsertExit::End)
        } else {
            Ok(EndInsertExit::Exit)
        }
    }

    /// Selects the layer of the pHMM to enter from the [`SemiLocalModule`] at
    /// the start of the pHMM.
    ///
    /// The `score` for the entire [`LocalModule`] at the start of the pHMM is
    /// added in this method.
    ///
    /// ## Errors
    ///
    /// [`InvalidCigarOp`] is returned if the operation is not in `MDI=XS`.
    ///
    /// [`InvalidCigarOp`]: LocalTraverseFromAlignError::InvalidCigarOp
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    fn enter_core(
        &mut self, module: &SemiLocalModule<T>, phmm: &LocalPhmm<T, S>,
    ) -> Result<DpIndex, LocalTraverseFromAlignError> {
        let layer = if let Some(empty_info) = &self.empty_info {
            if empty_info.through_begin {
                Begin.to_dp_index()
            } else {
                End.to_dp_index(phmm)
            }
        } else {
            self.inner.enter_core(module, phmm)?
        };

        self.inner.inner.score += phmm.get_begin_semilocal_score(layer);
        Ok(layer)
    }

    /// Performs any behavior necessary given that the traversal is exiting the
    /// END state into the [`LocalModule`] at the end of the pHMM.
    ///
    /// ## Errors
    ///
    /// This implementation is infallible.
    ///
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    fn exit_core_from_end(
        &mut self, _layer: DpIndex, _exit_param: T, _phmm: &LocalPhmm<T, S>,
    ) -> Result<(), LocalTraverseFromAlignError> {
        Ok(())
    }

    /// Selects the index of the residue that is emitted within either
    /// [`LocalModule`].
    ///
    /// ## Errors
    ///
    /// [`QueryLenMismatch`] is returned if there are no remaining residues in
    /// the query.
    ///
    /// [`QueryLenMismatch`]: DomainTraverseFromAlignError::QueryLenMismatch
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    fn choose_local_emission(
        &mut self, _params: &EmissionParams<T, S>, mapping: &ByteIndexMap<S>, loc: ModuleLocation, _phmm: &LocalPhmm<T, S>,
    ) -> Result<usize, LocalTraverseFromAlignError> {
        let byte = match loc {
            ModuleLocation::Begin => self.begin_residues.split_off_first(),
            ModuleLocation::End => self.end_residues.split_off_first(),
        }
        .ok_or_else(|| self.inner.inner.context.query_len_mismatch())?;

        Ok(mapping.to_index(*byte))
    }

    /// Selects whether to enter the insert state within either [`LocalModule`].
    ///
    /// If `true` is returned, then the insert state is entered. If `false` is
    /// returned, traversal continues to the end of the [`LocalModule`].
    ///
    /// If this is called for the module at the start of the pHMM, then the
    /// parameter for the domain module is stored to `domain_start_param`. If
    /// this is called for the module at the end of the pHMM, then both
    /// `exit_param` and the parameter for the domain module are accumulated
    /// into the `score` field of the inner visitor.
    ///
    /// ## Errors
    ///
    /// This implementation is infallible.
    ///
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    fn enter_module_insert(
        &mut self, _module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<bool, LocalTraverseFromAlignError> {
        match loc {
            ModuleLocation::Begin => {
                self.inner.inner.score += phmm.get_begin_domain_score(self.begin_residues);
                Ok(!self.begin_residues.is_empty())
            }
            ModuleLocation::End => Ok(!self.end_residues.is_empty()),
        }
    }

    /// Selects whether to exit the insert state within either [`LocalModule`].
    ///
    /// If `true` is returned, then the insert state is exited and traversal
    /// continues to the end of the [`LocalModule`]. If `false` is returned,
    /// then traversal stays within the insert state.
    ///
    /// ## Errors
    ///
    /// This implementation is infallible.
    ///
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    fn exit_module_insert(
        &mut self, _module: &DomainModule<T, S>, loc: ModuleLocation, _phmm: &LocalPhmm<T, S>,
    ) -> Result<bool, LocalTraverseFromAlignError> {
        match loc {
            ModuleLocation::Begin => Ok(self.begin_residues.is_empty()),
            ModuleLocation::End => Ok(self.end_residues.is_empty()),
        }
    }

    fn exiting_domain_module(
        &mut self, _module: &DomainModule<T, S>, _loc: ModuleLocation, _phmm: &LocalPhmm<T, S>,
    ) -> Result<(), Self::Error> {
        // TODO: Likely want to do score updates here for clarity, rather than
        // in enter_module_insert
        Ok(())
    }

    fn exit_core(
        &mut self, _layer_idx: DpIndex, exit_param: T, phmm: &LocalPhmm<T, S>,
    ) -> Result<(), LocalTraverseFromAlignError> {
        self.inner.inner.score += self.compute_end_module_score(exit_param, phmm);
        Ok(())
    }

    /// Returns the final score, as well as checking that the traversal ended in
    /// a valid state.
    ///
    /// ## Errors
    ///
    /// If an additional operation is present, then one of [`InvalidCigarOp`],
    /// [`DuplicateOp`], [`InternalClipping`], [`QueryLenMismatch`], or
    /// [`ModelLenMismatch`] is returned depending on the cause. If the full
    /// query wasn't consumed, then [`QueryLenMismatch`] is returned.
    ///
    /// [`InvalidCigarOp`]: LocalTraverseFromAlignError::InvalidCigarOp
    /// [`DuplicateOp`]: LocalTraverseFromAlignError::DuplicateOp
    /// [`InternalClipping`]: LocalTraverseFromAlignError::InternalClipping
    /// [`QueryLenMismatch`]: LocalTraverseFromAlignError::QueryLenMismatch
    /// [`ModelLenMismatch`]: LocalTraverseFromAlignError::ModelLenMismatch
    fn finalize(
        mut self, _phmm: &LocalPhmm<T, S>, _aligned_layers: RangeInclusive<DpIndex>, _aligned_seq: Range<SeqIndex>,
    ) -> Result<T, LocalTraverseFromAlignError> {
        if let Some(op) = self.inner.inner.op_iter.next() {
            return Err(self.inner.inner.context.remaining_op_error(op));
        }

        if self.inner.inner.query_idx != self.inner.inner.query.len() {
            return Err(self.inner.inner.context.query_len_mismatch());
        }

        Ok(self.inner.inner.score)
    }
}

impl PhmmState {
    /// Converts a CIGAR-style operation to a [`PhmmState`], return `None` if it
    /// is not in `MDI=X`.
    #[inline]
    pub(crate) fn from_op(op: u8) -> Option<Self> {
        match op {
            b'M' | b'=' | b'X' => Some(PhmmState::Match),
            b'D' => Some(PhmmState::Delete),
            b'I' => Some(PhmmState::Insert),
            _ => None,
        }
    }
}

impl PhmmStateOrModule {
    /// Converts a CIGAR-style operation to a [`PhmmState`], return `None` if it
    /// is not in `MDIS=X`.
    #[inline]
    pub(crate) fn from_op(op: u8) -> Option<Self> {
        match op {
            b'M' | b'=' | b'X' => Some(PhmmStateOrModule::Match),
            b'D' => Some(PhmmStateOrModule::Delete),
            b'I' => Some(PhmmStateOrModule::Insert),
            b'S' => Some(PhmmStateOrModule::Module),
            _ => None,
        }
    }
}
