//! The visitor traits for each pHMM type.

use crate::{
    alignment::phmm::{
        DomainPhmm, GlobalPhmm, LocalPhmm, SemiLocalPhmm,
        components::{EmissionParams, TransitionParams},
        indexing::{DpIndex, SeqIndex},
        modules::{DomainModule, SemiLocalModule},
        state::{PhmmState, PhmmStateOrModule},
        traverse::{EndInsert, EndInsertExit, ModuleLocation},
    },
    data::ByteIndexMap,
};
use std::ops::{Range, RangeInclusive};

/// A visitor over a [`GlobalPhmm`], for use with [`traverse_global_phmm`].
///
/// [`GlobalPhmm`]: crate::alignment::phmm::models::GlobalPhmm
/// [`traverse_global_phmm`]:
///     crate::alignment::phmm::traverse::traverse_global_phmm
pub trait GlobalVisitor<T, const S: usize> {
    /// The type output upon finalization of traversal.
    type Output;

    /// The type of error used by the visitor, if fallible traversal is
    /// necessary.
    type Error;

    /// Selects the index of the residue that is emitted at a match or insert
    /// state in the pHMM.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn choose_emission(
        &mut self, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        phmm: &GlobalPhmm<T, S>,
    ) -> Result<usize, Self::Error>;

    /// Selects the transition within the pHMM to take (what the next
    /// [`PhmmState`] to enter should be).
    ///
    /// When in the [`LastMatch`] layer, [`choose_end_or_insert`] is called
    /// instead.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_end_or_insert`]: GlobalVisitor::choose_end_or_insert
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &GlobalPhmm<T, S>,
    ) -> Result<PhmmState, Self::Error>;

    /// Selects whether the traversal should enter the [`End`] state from any of
    /// the states in the [`LastMatch`] layer, or whether the final insert state
    /// should be entered.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &GlobalPhmm<T, S>,
    ) -> Result<EndInsert, Self::Error>;

    /// Finalizes traversal and performs any checking for the visitor.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn finalize(self, phmm: &GlobalPhmm<T, S>) -> Result<Self::Output, Self::Error>;
}

/// A visitor over a [`SemiLocalPhmm`], for use with
/// [`traverse_semilocal_phmm`].
///
/// [`SemiLocalPhmm`]: crate::alignment::phmm::SemiLocalPhmm
/// [`traverse_semilocal_phmm`]:
///     crate::alignment::phmm::traverse::traverse_semilocal_phmm
pub trait SemiLocalVisitor<T, const S: usize> {
    /// The type output upon finalization of traversal.
    type Output;

    /// The type of error used by the visitor, if fallible traversal is
    /// necessary.
    type Error;

    /// Selects the index of the residue that is emitted at a match or insert
    /// state in the pHMM.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn choose_emission(
        &mut self, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<usize, Self::Error>;

    /// Selects the transition within the core pHMM to take (what the next
    /// [`PhmmState`] to enter should be) when in an insert or delete state.
    ///
    /// When in a match state, [`choose_core_transition_or_exit`] is called
    /// instead. When in the [`LastMatch`] layer, [`choose_end_or_insert`] is
    /// called instead.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_core_transition_or_exit`]:
    ///     SemiLocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_or_insert`]: SemiLocalVisitor::choose_end_or_insert
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<PhmmState, Self::Error>;

    /// Selects whether the traversal should enter the [`End`] state from the
    /// insert or delete states in the [`LastMatch`] layer, or whether the final
    /// insert state should be entered.
    ///
    /// When in the match state of the [`LastMatch`] layer,
    /// [`choose_end_insert_or_exit`] is called instead.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_end_insert_or_exit`]:
    ///     SemiLocalVisitor::choose_end_insert_or_exit
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<EndInsert, Self::Error>;

    /// Selects the transition within the core pHMM to take (what the next
    /// [`PhmmState`] to enter should be), or whether to exit from the pHMM
    /// early.
    ///
    /// This is only called when in a match state.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn choose_core_transition_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<PhmmStateOrModule, Self::Error>;

    /// From the match state in the [`LastMatch`] layer, selects whether the
    /// traversal should enter the [`End`] state, enter the final insert state,
    /// or exit early from the pHMM.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    fn choose_end_insert_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
        phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<EndInsertExit, Self::Error>;

    /// Selects the layer of the pHMM to enter from the [`SemiLocalModule`] at
    /// the start of the pHMM.
    ///
    /// ## Validity
    ///
    /// This must return an index that is in-range for the given pHMM, otherwise
    /// [`traverse_semilocal_phmm`] may panic.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`traverse_semilocal_phmm`]:
    ///     crate::alignment::phmm::traverse::traverse_semilocal_phmm
    fn enter_core(&mut self, module: &SemiLocalModule<T>, phmm: &SemiLocalPhmm<T, S>) -> Result<DpIndex, Self::Error>;

    /// Performs any behavior necessary given that the traversal is exiting the
    /// END state into the [`SemiLocalModule`] at the end of the pHMM.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn exit_core_from_end(&mut self, layer: DpIndex, exit_param: T, phmm: &SemiLocalPhmm<T, S>) -> Result<(), Self::Error>;

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
    /// See the implementor for documentation of possible errors.
    ///
    /// [`choose_core_transition_or_exit`]:
    ///     SemiLocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_insert_or_exit`]:
    ///     SemiLocalVisitor::choose_end_insert_or_exit
    /// [`exit_core_from_end`]: SemiLocalVisitor::exit_core_from_end
    fn exit_core(&mut self, layer_idx: DpIndex, exit_param: T, phmm: &SemiLocalPhmm<T, S>) -> Result<(), Self::Error>;

    /// Finalizes traversal and performs any checking for the visitor.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn finalize(
        self, phmm: &SemiLocalPhmm<T, S>, aligned_layers: RangeInclusive<DpIndex>,
    ) -> Result<Self::Output, Self::Error>;
}

/// A visitor over a [`DomainPhmm`], for use with [`traverse_domain_phmm`].
///
/// [`DomainPhmm`]: crate::alignment::phmm::models::DomainPhmm
/// [`traverse_domain_phmm`]:
///     crate::alignment::phmm::traverse::traverse_domain_phmm
pub trait DomainVisitor<T, const S: usize> {
    /// The type output upon finalization of traversal.
    type Output;

    /// The type of error used by the visitor, if fallible traversal is
    /// necessary.
    type Error;

    /// Selects the index of the residue that is emitted at a match or insert
    /// state in the core pHMM.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn choose_emission(
        &mut self, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        phmm: &DomainPhmm<T, S>,
    ) -> Result<usize, Self::Error>;

    /// Selects the transition within the pHMM to take (what the next
    /// [`PhmmState`] to enter should be).
    ///
    /// When in the [`LastMatch`] layer, [`choose_end_or_insert`] is called
    /// instead.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_end_or_insert`]: DomainVisitor::choose_end_or_insert
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &DomainPhmm<T, S>,
    ) -> Result<PhmmState, Self::Error>;

    /// Selects whether the traversal should enter the [`End`] state from any of
    /// the states in the [`LastMatch`] layer, or whether the final insert state
    /// should be entered.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &DomainPhmm<T, S>,
    ) -> Result<EndInsert, Self::Error>;

    /// Selects the index of the residue that is emitted within either
    /// [`DomainModule`].
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn choose_domain_emission(
        &mut self, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<usize, Self::Error>;

    /// Selects whether to enter the insert state within either
    /// [`DomainModule`].
    ///
    /// If `true` is returned, then the insert state is entered. If `false` is
    /// returned, traversal continues to the end of the [`DomainModule`].
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn enter_module_insert(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<bool, Self::Error>;

    /// Selects whether to exit the insert state within either [`DomainModule`].
    ///
    /// If `true` is returned, then the insert state is exited and traversal
    /// continues to the end of the [`DomainModule`]. If `false` is returned,
    /// then traversal stays within the insert state.
    ///
    /// ## Validity
    ///
    /// This function should eventually return `true` to ensure
    /// [`traverse_domain_phmm`] is not stuck in an infinite loop.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`traverse_domain_phmm`]:
    ///     crate::alignment::phmm::traverse::traverse_domain_phmm
    fn exit_module_insert(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<bool, Self::Error>;

    /// A hook called when the end of a [`DomainModule`] is reached.
    ///
    /// This is called directly after [`exit_module_insert`] returns `true` or
    /// [`enter_module_insert`] returns `false`.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`exit_module_insert`]: DomainVisitor::exit_module_insert
    /// [`enter_module_insert`]: DomainVisitor::enter_module_insert
    fn exiting_module(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<(), Self::Error>;

    /// Finalizes traversal and performs any checking for the visitor.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn finalize(self, phmm: &DomainPhmm<T, S>, query_range: Range<SeqIndex>) -> Result<Self::Output, Self::Error>;
}

/// A visitor over a [`LocalPhmm`], for use with [`traverse_local_phmm`].
///
/// [`LocalPhmm`]: crate::alignment::phmm::models::LocalPhmm
/// [`traverse_local_phmm`]:
///     crate::alignment::phmm::traverse::traverse_local_phmm
pub trait LocalVisitor<T, const S: usize> {
    /// The type output upon finalization of traversal.
    type Output;

    /// The type of error used by the visitor, if fallible traversal is
    /// necessary.
    type Error;

    /// Selects the index of the residue that is emitted at a match or insert
    /// state in the pHMM.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.    
    fn choose_emission(
        &mut self, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        phmm: &LocalPhmm<T, S>,
    ) -> Result<usize, Self::Error>;

    /// Selects the transition within the core pHMM to take (what the next
    /// [`PhmmState`] to enter should be) when in an insert or delete state.
    ///
    /// When in a match state, [`choose_core_transition_or_exit`] is called
    /// instead. When in the [`LastMatch`] layer, [`choose_end_or_insert`] is
    /// called instead.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_core_transition_or_exit`]:
    ///     LocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_or_insert`]: LocalVisitor::choose_end_or_insert
    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &LocalPhmm<T, S>,
    ) -> Result<PhmmState, Self::Error>;

    /// Selects whether the traversal should enter the [`End`] state from the
    /// insert or delete states in the [`LastMatch`] layer, or whether the final
    /// insert state should be entered.
    ///
    /// When in the match state of the [`LastMatch`] layer,
    /// [`choose_end_insert_or_exit`] is called instead.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`choose_end_insert_or_exit`]: LocalVisitor::choose_end_insert_or_exit
    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &LocalPhmm<T, S>,
    ) -> Result<EndInsert, Self::Error>;

    /// Selects the transition within the core pHMM to take (what the next
    /// [`PhmmState`] to enter should be), or whether to exit from the pHMM
    /// early.
    ///
    /// This is only called when in a match state.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn choose_core_transition_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, phmm: &LocalPhmm<T, S>,
    ) -> Result<PhmmStateOrModule, Self::Error>;

    /// From the match state in the [`LastMatch`] layer, selects whether the
    /// traversal should enter the [`End`] state, enter the final insert state,
    /// or exit early from the pHMM.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`LastMatch`]: crate::alignment::phmm::indexing::LastMatch
    /// [`End`]: crate::alignment::phmm::indexing::End
    fn choose_end_insert_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
        phmm: &LocalPhmm<T, S>,
    ) -> Result<EndInsertExit, Self::Error>;

    /// Selects the layer of the pHMM to enter from the [`SemiLocalModule`] at
    /// the start of the pHMM.
    ///
    /// ## Validity
    ///
    /// This must return an index that is in-range for the given pHMM, otherwise
    /// [`traverse_local_phmm`] may panic.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`traverse_local_phmm`]:
    ///     crate::alignment::phmm::traverse::traverse_local_phmm
    fn enter_core(&mut self, module: &SemiLocalModule<T>, phmm: &LocalPhmm<T, S>) -> Result<DpIndex, Self::Error>;

    /// Performs any behavior necessary given that the traversal is exiting the
    /// END state into the [`SemiLocalModule`] at the end of the pHMM.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn exit_core_from_end(&mut self, layer: DpIndex, exit_param: T, phmm: &LocalPhmm<T, S>) -> Result<(), Self::Error>;

    /// Selects the index of the residue that is emitted within either
    /// [`LocalModule`].
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    fn choose_local_emission(
        &mut self, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<usize, Self::Error>;

    /// Selects whether to enter the insert state within either [`LocalModule`].
    ///
    /// If `true` is returned, then the insert state is entered. If `false` is
    /// returned, traversal continues to the end of the [`LocalModule`].
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    fn enter_module_insert(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<bool, Self::Error>;

    /// Selects whether to exit the insert state within either [`LocalModule`].
    ///
    /// If `true` is returned, then the insert state is exited and traversal
    /// continues to the end of the [`LocalModule`]. If `false` is returned,
    /// then traversal stays within the insert state.
    ///
    /// ## Validity
    ///
    /// This function should eventually return `true` to ensure
    /// [`traverse_local_phmm`] is not stuck in an infinite loop.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    /// [`traverse_local_phmm`]:
    ///     crate::alignment::phmm::traverse::traverse_local_phmm
    fn exit_module_insert(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<bool, Self::Error>;

    /// A hook called when the end of a [`DomainModule`] is reached.
    ///
    /// This is called directly after [`exit_module_insert`] returns `true` or
    /// [`enter_module_insert`] returns `false`.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    ///
    /// [`exit_module_insert`]: DomainVisitor::exit_module_insert
    /// [`enter_module_insert`]: DomainVisitor::enter_module_insert
    fn exiting_domain_module(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<(), Self::Error>;

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
    /// See the implementor for documentation of possible errors.
    ///
    /// [`choose_core_transition_or_exit`]:
    ///     LocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_insert_or_exit`]: LocalVisitor::choose_end_insert_or_exit
    /// [`exit_core_from_end`]: LocalVisitor::exit_core_from_end
    fn exit_core(&mut self, layer_idx: DpIndex, exit_param: T, phmm: &LocalPhmm<T, S>) -> Result<(), Self::Error>;

    /// Finalizes traversal and performs any checking for the visitor.
    ///
    /// ## Errors
    ///
    /// See the implementor for documentation of possible errors.
    fn finalize(
        self, phmm: &LocalPhmm<T, S>, aligned_layers: RangeInclusive<DpIndex>, query_range: Range<SeqIndex>,
    ) -> Result<Self::Output, Self::Error>;
}
