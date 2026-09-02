//! Visitors for recording the Viterbi score achieved while traversing a pHMM.

use crate::{
    alignment::phmm::{
        DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
        components::{EmissionParams, TransitionParams},
        indexing::{DpIndex, SeqIndex},
        modules::{DomainModule, SemiLocalModule},
        state::{PhmmState, PhmmStateOrModule},
        traverse::{DomainVisitor, EndInsert, EndInsertExit, GlobalVisitor, LocalVisitor, ModuleLocation, SemiLocalVisitor},
    },
    data::ByteIndexMap,
};
use std::ops::{Range, RangeInclusive};

/// A wrapper around a visitor that tracks the score of the path taken.
///
/// This wrapper does not produce any errors of its own. The only errors are
/// those yielded by the inner visitor.
///
/// ## Parameters
///
/// - `V` - The type of the wrapped visitor.
/// - `T` - The parameter type used by the pHMM.
pub struct ScoreVisitor<V, T> {
    /// The visitor used to decide the traversal.
    inner_visitor:   V,
    /// The score of the traversal.
    ///
    /// This is updated immediately for core transitions/emissions. Domain
    /// emissions are accumulated into `module_inserted`, and the score is added
    /// here only when the [`DomainModule`] exits. The transition into the core
    /// pHMM is added immediately for a [`SemiLocalModule`] or [`LocalModule`].
    /// The transition out of the core pHMM is added immediately for a
    /// [`SemiLocalPhmm`] in [`exit_core`], or for [`LocalPhmm`], it is stored
    /// in `exit_param` and added when the [`DomainModule`] exits.
    ///
    /// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
    /// [`SemiLocalPhmm`]: crate::alignment::phmm::SemiLocalPhmm
    /// [`exit_core`]:
    ///     crate::alignment::phmm::traverse::SemiLocalVisitor::exit_core
    /// [`LocalPhmm`]: crate::alignment::phmm::LocalPhmm
    score:           T,
    /// In the case of a [`DomainPhmm`] or [`LocalPhmm`], this holds the
    /// residues that are emitted by the current [`DomainModule`]. Otherwise,
    /// this is empty and does not cause an allocation.
    ///
    /// [`DomainPhmm`]: crate::alignment::phmm::DomainPhmm
    /// [`LocalPhmm`]: crate::alignment::phmm::LocalPhmm
    module_inserted: Vec<u8>,
    /// In the case of a [`LocalPhmm`], this holds the parameter used to exit
    /// the core pHMM. Otherwise, this is zero.
    ///
    /// [`LocalPhmm`]: crate::alignment::phmm::LocalPhmm
    exit_param:      T,
}

impl<V, T> ScoreVisitor<V, T>
where
    T: PhmmNumber,
{
    /// Constructs a new [`ScoreVisitor`] that tracks the score as traversal
    /// using `inner_visitor` takes place.
    pub fn new(inner_visitor: V) -> Self {
        Self {
            inner_visitor,
            score: T::ZERO,
            module_inserted: Vec::new(),
            exit_param: T::ZERO,
        }
    }
}

/// The output of a [`ScoreVisitor`] wrapper, containing the original output and
/// the calculated score.
pub struct WithScore<A, T> {
    pub output: A,
    pub score:  T,
}

impl<V, T, const S: usize> GlobalVisitor<T, S> for ScoreVisitor<V, T>
where
    T: PhmmNumber,
    V: GlobalVisitor<T, S>,
{
    type Output = WithScore<V::Output, T>;
    type Error = V::Error;

    fn choose_emission(
        &mut self, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        phmm: &GlobalPhmm<T, S>,
    ) -> Result<usize, Self::Error> {
        let idx = self.inner_visitor.choose_emission(layer, state, params, map, phmm)?;
        self.score += params[idx];
        Ok(idx)
    }

    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &GlobalPhmm<T, S>,
    ) -> Result<PhmmState, Self::Error> {
        let next_state = self.inner_visitor.choose_core_transition(layer, exiting, params, phmm)?;
        self.score += params[(exiting, next_state)];
        Ok(next_state)
    }

    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &GlobalPhmm<T, S>,
    ) -> Result<EndInsert, Self::Error> {
        let next_state = self.inner_visitor.choose_end_or_insert(layer, exiting, params, phmm)?;
        let param = match next_state {
            EndInsert::End => params[(exiting, PhmmState::Match)],
            EndInsert::Insert => params[(exiting, PhmmState::Insert)],
        };
        self.score += param;
        Ok(next_state)
    }

    fn finalize(self, phmm: &GlobalPhmm<T, S>) -> Result<Self::Output, Self::Error> {
        Ok(WithScore {
            output: self.inner_visitor.finalize(phmm)?,
            score:  self.score,
        })
    }
}

impl<V, T, const S: usize> DomainVisitor<T, S> for ScoreVisitor<V, T>
where
    T: PhmmNumber,
    V: DomainVisitor<T, S>,
{
    type Output = WithScore<V::Output, T>;
    type Error = V::Error;

    fn choose_emission(
        &mut self, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        phmm: &DomainPhmm<T, S>,
    ) -> Result<usize, Self::Error> {
        let idx = self.inner_visitor.choose_emission(layer, state, params, map, phmm)?;
        self.score += params[idx];
        Ok(idx)
    }

    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &DomainPhmm<T, S>,
    ) -> Result<PhmmState, Self::Error> {
        let next_state = self.inner_visitor.choose_core_transition(layer, exiting, params, phmm)?;
        self.score += params[(exiting, next_state)];
        Ok(next_state)
    }

    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &DomainPhmm<T, S>,
    ) -> Result<EndInsert, Self::Error> {
        let next_state = self.inner_visitor.choose_end_or_insert(layer, exiting, params, phmm)?;
        let param = match next_state {
            EndInsert::End => params[(exiting, PhmmState::Match)],
            EndInsert::Insert => params[(exiting, PhmmState::Insert)],
        };
        self.score += param;
        Ok(next_state)
    }

    fn choose_domain_emission(
        &mut self, params: &EmissionParams<T, S>, mapping: &ByteIndexMap<S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<usize, Self::Error> {
        // This function does not update the score, since that is done later in
        // exit_module_insert. However, it does update module_inserted.
        let idx = self.inner_visitor.choose_domain_emission(params, mapping, loc, phmm)?;
        self.module_inserted.push(mapping.byte_keys()[idx]);
        Ok(idx)
    }

    fn enter_module_insert(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<bool, Self::Error> {
        // This function does not update the score, since that is done later in
        // exiting_module
        self.inner_visitor.enter_module_insert(module, loc, phmm)
    }

    fn exit_module_insert(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<bool, Self::Error> {
        // This function does not update the score, since that is done later in
        // exiting_module
        let exit_module_insert = self.inner_visitor.exit_module_insert(module, loc, phmm)?;

        Ok(exit_module_insert)
    }

    fn exiting_module(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &DomainPhmm<T, S>,
    ) -> Result<(), Self::Error> {
        // Update score all at once to have correct order of floating point
        // operations
        match loc {
            ModuleLocation::Begin => self.score += module.get_begin_score(&self.module_inserted, phmm.mapping()),
            ModuleLocation::End => self.score += module.get_end_score(&self.module_inserted, phmm.mapping()),
        }
        self.module_inserted.clear();

        self.inner_visitor.exiting_module(module, loc, phmm)
    }

    fn finalize(self, phmm: &DomainPhmm<T, S>, query_range: Range<SeqIndex>) -> Result<Self::Output, Self::Error> {
        Ok(WithScore {
            output: self.inner_visitor.finalize(phmm, query_range)?,
            score:  self.score,
        })
    }
}

impl<V, T, const S: usize> SemiLocalVisitor<T, S> for ScoreVisitor<V, T>
where
    T: PhmmNumber,
    V: SemiLocalVisitor<T, S>,
{
    type Output = WithScore<V::Output, T>;
    type Error = V::Error;

    fn choose_emission(
        &mut self, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<usize, Self::Error> {
        let idx = self.inner_visitor.choose_emission(layer, state, params, map, phmm)?;
        self.score += params[idx];
        Ok(idx)
    }

    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<PhmmState, Self::Error> {
        let next_state = self.inner_visitor.choose_core_transition(layer, exiting, params, phmm)?;
        self.score += params[(exiting, next_state)];
        Ok(next_state)
    }

    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<EndInsert, Self::Error> {
        let next_state = self.inner_visitor.choose_end_or_insert(layer, exiting, params, phmm)?;
        let param = match next_state {
            EndInsert::End => params[(exiting, PhmmState::Match)],
            EndInsert::Insert => params[(exiting, PhmmState::Insert)],
        };
        self.score += param;
        Ok(next_state)
    }

    fn choose_core_transition_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<PhmmStateOrModule, Self::Error> {
        let next_state = self
            .inner_visitor
            .choose_core_transition_or_exit(layer, params, exit_param, phmm)?;

        if let Some(state) = PhmmState::get_from(next_state) {
            self.score += params[(PhmmState::Match, state)];
        }

        Ok(next_state)
    }

    fn choose_end_insert_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
        phmm: &SemiLocalPhmm<T, S>,
    ) -> Result<EndInsertExit, Self::Error> {
        let next_state =
            self.inner_visitor
                .choose_end_insert_or_exit(layer, params, exit_param, exit_from_end_param, phmm)?;

        match next_state {
            EndInsertExit::End => self.score += params[(PhmmState::Match, PhmmState::Match)],
            EndInsertExit::Insert => self.score += params[(PhmmState::Match, PhmmState::Insert)],
            EndInsertExit::Exit => {}
        }

        Ok(next_state)
    }

    fn enter_core(&mut self, module: &SemiLocalModule<T>, phmm: &SemiLocalPhmm<T, S>) -> Result<DpIndex, Self::Error> {
        let layer = self.inner_visitor.enter_core(module, phmm)?;
        self.score += module.get_score(layer);
        Ok(layer)
    }

    fn exit_core_from_end(&mut self, layer: DpIndex, exit_param: T, phmm: &SemiLocalPhmm<T, S>) -> Result<(), Self::Error> {
        self.inner_visitor.exit_core_from_end(layer, exit_param, phmm)?;
        Ok(())
    }

    fn exit_core(&mut self, layer_idx: DpIndex, exit_param: T, phmm: &SemiLocalPhmm<T, S>) -> Result<(), Self::Error> {
        self.inner_visitor.exit_core(layer_idx, exit_param, phmm)?;
        self.score += exit_param;
        Ok(())
    }

    fn finalize(
        self, phmm: &SemiLocalPhmm<T, S>, aligned_layers: RangeInclusive<DpIndex>,
    ) -> Result<Self::Output, Self::Error> {
        Ok(WithScore {
            output: self.inner_visitor.finalize(phmm, aligned_layers)?,
            score:  self.score,
        })
    }
}

impl<V, T, const S: usize> LocalVisitor<T, S> for ScoreVisitor<V, T>
where
    T: PhmmNumber,
    V: LocalVisitor<T, S>,
{
    type Output = WithScore<V::Output, T>;
    type Error = V::Error;

    fn choose_emission(
        &mut self, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
        phmm: &LocalPhmm<T, S>,
    ) -> Result<usize, Self::Error> {
        let idx = self.inner_visitor.choose_emission(layer, state, params, map, phmm)?;
        self.score += params[idx];
        Ok(idx)
    }

    fn choose_core_transition(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &LocalPhmm<T, S>,
    ) -> Result<PhmmState, Self::Error> {
        let next_state = self.inner_visitor.choose_core_transition(layer, exiting, params, phmm)?;
        self.score += params[(exiting, next_state)];
        Ok(next_state)
    }

    fn choose_end_or_insert(
        &mut self, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>, phmm: &LocalPhmm<T, S>,
    ) -> Result<EndInsert, Self::Error> {
        let next_state = self.inner_visitor.choose_end_or_insert(layer, exiting, params, phmm)?;
        let param = match next_state {
            EndInsert::End => params[(exiting, PhmmState::Match)],
            EndInsert::Insert => params[(exiting, PhmmState::Insert)],
        };
        self.score += param;
        Ok(next_state)
    }

    fn choose_core_transition_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, phmm: &LocalPhmm<T, S>,
    ) -> Result<PhmmStateOrModule, Self::Error> {
        let next_state = self
            .inner_visitor
            .choose_core_transition_or_exit(layer, params, exit_param, phmm)?;

        if let Some(state) = PhmmState::get_from(next_state) {
            self.score += params[(PhmmState::Match, state)];
        }

        Ok(next_state)
    }

    fn choose_end_insert_or_exit(
        &mut self, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
        phmm: &LocalPhmm<T, S>,
    ) -> Result<EndInsertExit, Self::Error> {
        let next_state =
            self.inner_visitor
                .choose_end_insert_or_exit(layer, params, exit_param, exit_from_end_param, phmm)?;

        match next_state {
            EndInsertExit::End => self.score += params[(PhmmState::Match, PhmmState::Match)],
            EndInsertExit::Insert => self.score += params[(PhmmState::Match, PhmmState::Insert)],
            EndInsertExit::Exit => {}
        }

        Ok(next_state)
    }

    fn enter_core(&mut self, module: &SemiLocalModule<T>, phmm: &LocalPhmm<T, S>) -> Result<DpIndex, Self::Error> {
        let layer = self.inner_visitor.enter_core(module, phmm)?;
        self.score += module.get_score(layer);
        Ok(layer)
    }

    fn exit_core_from_end(&mut self, layer: DpIndex, exit_param: T, phmm: &LocalPhmm<T, S>) -> Result<(), Self::Error> {
        self.inner_visitor.exit_core_from_end(layer, exit_param, phmm)?;
        Ok(())
    }

    fn exit_core(&mut self, layer_idx: DpIndex, exit_param: T, phmm: &LocalPhmm<T, S>) -> Result<(), Self::Error> {
        self.inner_visitor.exit_core(layer_idx, exit_param, phmm)?;
        self.exit_param = exit_param;
        Ok(())
    }

    fn choose_local_emission(
        &mut self, params: &EmissionParams<T, S>, mapping: &ByteIndexMap<S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<usize, Self::Error> {
        // This function does not update the score, since that is done later in
        // exit_module_insert. However, it does update module_inserted.
        let idx = self.inner_visitor.choose_local_emission(params, mapping, loc, phmm)?;
        self.module_inserted.push(mapping.byte_keys()[idx]);
        Ok(idx)
    }

    fn enter_module_insert(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<bool, Self::Error> {
        // This function does not update the score, since that is done later in
        // exit_module_insert
        self.inner_visitor.enter_module_insert(module, loc, phmm)
    }

    fn exit_module_insert(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<bool, Self::Error> {
        // This function does not update the score, since that is done later in
        // exiting_module
        let exit_module_insert = self.inner_visitor.exit_module_insert(module, loc, phmm)?;

        Ok(exit_module_insert)
    }

    fn exiting_domain_module(
        &mut self, module: &DomainModule<T, S>, loc: ModuleLocation, phmm: &LocalPhmm<T, S>,
    ) -> Result<(), Self::Error> {
        // Update score all at once to have correct order of floating point
        // operations
        let domain_score = match loc {
            ModuleLocation::Begin => module.get_begin_score(&self.module_inserted, phmm.mapping()),
            ModuleLocation::End => module.get_end_score(&self.module_inserted, phmm.mapping()),
        };
        self.score += match loc {
            ModuleLocation::Begin => domain_score,
            ModuleLocation::End => domain_score + self.exit_param,
        };
        self.module_inserted.clear();

        self.inner_visitor.exiting_domain_module(module, loc, phmm)
    }

    fn finalize(
        self, phmm: &LocalPhmm<T, S>, aligned_layers: RangeInclusive<DpIndex>, query_range: Range<SeqIndex>,
    ) -> Result<Self::Output, Self::Error> {
        Ok(WithScore {
            output: self.inner_visitor.finalize(phmm, aligned_layers, query_range)?,
            score:  self.score,
        })
    }
}
