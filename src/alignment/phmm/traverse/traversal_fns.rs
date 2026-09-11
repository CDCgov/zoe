//! Functions for traversing pHMMs with a given visitor.

use crate::alignment::phmm::{
    DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
    indexing::{AlnIndex, AlnIndexable, Begin, DpIndex, End, GetLayer, GetMapping, GetModule, SeqIndex},
    modules::{DomainParams, SemiLocalParams},
    state::PhmmState,
    traverse::{
        DomainVisitor, EndInsert, EndInsertExit, GlobalVisitor, LocalVisitor, ModuleLocation, SemiLocalVisitor, VisitCore,
        VisitCoreOrExit, VisitDomainModule,
    },
};
use std::ops::{Range, RangeInclusive};

/// Traverses the core pHMM using a given visitor, for use with [`GlobalPhmm`]
/// or [`DomainPhmm`] (no early exit is permitted).
///
/// ## Parameters
///
/// - `V`: The type of the visitor
/// - `P`: The type of pHMM being traversed
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
///
/// [`CorePhmm`]: crate::alignment::phmm::components::CorePhmm
fn traverse_core_phmm<V, P, T, const S: usize>(
    phmm: &P, mut num_emitted: usize, visitor: &mut V,
) -> Result<usize, P::Error>
where
    P: GetMapping<S> + GetLayer<T, S> + VisitCore<V, T, S>, {
    let (mut layer, mut remaining_layers) = phmm.layers().split_first();
    let mut layer_idx = Begin.to_dp_index();
    let mut state = PhmmState::Match;

    while let Some((next_layer, rest)) = remaining_layers.split_first() {
        let params = &layer.transition;
        let next_state = phmm.choose_core_transition(visitor, layer_idx, state, params)?;

        match next_state {
            PhmmState::Match => {
                let emit_params = &layer.emission_match;
                layer_idx = layer_idx.next_index(&phmm);
                layer = next_layer;
                remaining_layers = rest;

                phmm.choose_emission(visitor, layer_idx, next_state, emit_params, phmm.mapping())?;
                num_emitted += 1;
            }
            PhmmState::Delete => {
                layer_idx = layer_idx.next_index(&phmm);
                layer = next_layer;
                remaining_layers = rest;
            }
            PhmmState::Insert => {
                let emit_params = &layer.emission_insert;

                phmm.choose_emission(visitor, layer_idx, next_state, emit_params, phmm.mapping())?;
                num_emitted += 1;
            }
        }

        state = next_state;
    }

    // Handle the last match layer
    loop {
        let params = &layer.transition;
        let next_state = phmm.choose_end_or_insert(visitor, layer_idx, state, params)?;
        match next_state {
            EndInsert::End => return Ok(num_emitted),
            EndInsert::Insert => {
                let emit_params = &layer.emission_insert;
                state = PhmmState::Insert;

                phmm.choose_emission(visitor, layer_idx, state, emit_params, phmm.mapping())?;
                num_emitted += 1;
            }
        }
    }
}

struct TraverseCorePhmmOrExitOutput<T> {
    exit_layer:     DpIndex,
    exit_param:     T,
    aligned_layers: RangeInclusive<DpIndex>,
    query_range:    Range<SeqIndex>,
}

/// Traverses the core pHMM using a given visitor, for use with
/// [`SemiLocalPhmm`] and [`LocalPhmm`] (early exit is permitted).
///
/// ## Errors
///
/// See the visitor's implementations for [`choose_emission`],
/// [`choose_core_transition`], [`choose_core_transition_or_exit`],
/// [`choose_end_or_insert`], [`choose_end_insert_or_exit`], and
/// [`exit_core_from_end`].
///
/// [`choose_emission`]: VisitCoreOrExit::choose_emission
/// [`choose_core_transition`]: VisitCoreOrExit::choose_core_transition
/// [`choose_core_transition_or_exit`]:
///     VisitCoreOrExit::choose_core_transition_or_exit
/// [`choose_end_or_insert`]: VisitCoreOrExit::choose_end_or_insert
/// [`choose_end_insert_or_exit`]: VisitCoreOrExit::choose_end_insert_or_exit
/// [`exit_core_from_end`]: VisitCoreOrExit::exit_core_from_end
fn traverse_core_phmm_or_exit<P, T, V, const S: usize>(
    phmm: &P, enter_layer: DpIndex, num_emitted_begin_module: usize, visitor: &mut V,
) -> Result<TraverseCorePhmmOrExitOutput<T>, P::Error>
where
    P: AlnIndexable + GetLayer<T, S> + GetMapping<S> + GetModule<End: SemiLocalParams<T>> + VisitCoreOrExit<V, T, S>,
    T: PhmmNumber + 'static, {
    let mut layer_idx = enter_layer;
    let mut num_emitted = num_emitted_begin_module;

    if layer_idx.eq_index(End, &phmm) {
        phmm.exit_core_from_end(visitor, layer_idx, phmm.end().semilocal_params().get_score(layer_idx))?;
        return Ok(TraverseCorePhmmOrExitOutput {
            exit_layer:     layer_idx,
            exit_param:     phmm.end().semilocal_params().get_score(layer_idx),
            aligned_layers: enter_layer..=layer_idx,
            query_range:    SeqIndex(num_emitted_begin_module)..SeqIndex(num_emitted),
        });
    }

    // Get the current layer (since we are not at End), any remaining
    // layers, and optionally the previous layer if we are not at Begin.
    let (prev_layer, mut layer, mut remaining_layers) = if let Some((before, layer_and_after)) =
        phmm.split_layers_at(layer_idx)
        && let Some((layer, remaining_layers)) = layer_and_after.split_first()
    {
        (before.last(), layer, remaining_layers)
    } else {
        // TODO: Convert this into an actual error variant, so now we need a
        // trait for visitor errors...
        panic!("The requested layer for entering the pHMM is out of bounds")
    };

    // The emission parameters are stored in the previous layer. If there is
    // no previous layer, then we are in Begin, which has no emission.
    if let Some(emit_params) = prev_layer.map(|layer| &layer.emission_match) {
        phmm.choose_emission(visitor, layer_idx, PhmmState::Match, emit_params, phmm.mapping())?;
        num_emitted += 1;
    }

    let mut state = PhmmState::Match;

    loop {
        // If the first branch doesn't run, then the current layer is the last
        // match layer
        if let Some((next_layer, rest)) = remaining_layers.split_first() {
            let params = &layer.transition;

            let next_state = if state == PhmmState::Match {
                let exit_param = phmm.end().semilocal_params().get_score(layer_idx);
                let next_state = phmm.choose_core_transition_or_exit(visitor, layer_idx, params, exit_param)?;

                match PhmmState::get_from(next_state) {
                    Some(next_state) => next_state,
                    None => break,
                }
            } else {
                phmm.choose_core_transition(visitor, layer_idx, state, params)?
            };

            match next_state {
                PhmmState::Match => {
                    let emit_params = &layer.emission_match;
                    layer_idx = layer_idx.next_index(&phmm);
                    layer = next_layer;
                    remaining_layers = rest;

                    phmm.choose_emission(visitor, layer_idx, next_state, emit_params, phmm.mapping())?;
                    num_emitted += 1;
                }
                PhmmState::Delete => {
                    layer_idx = layer_idx.next_index(&phmm);
                    layer = next_layer;
                    remaining_layers = rest;
                }
                PhmmState::Insert => {
                    phmm.choose_emission(visitor, layer_idx, next_state, &layer.emission_insert, phmm.mapping())?;
                    num_emitted += 1;
                }
            }

            state = next_state;
        } else if state == PhmmState::Match {
            let exit_param = phmm.end().semilocal_params().get_score(layer_idx);
            let exit_from_end_param = phmm.end().semilocal_params().get_score(End);
            let next_state =
                phmm.choose_end_insert_or_exit(visitor, layer_idx, &layer.transition, exit_param, exit_from_end_param)?;
            match next_state {
                EndInsertExit::End => {
                    layer_idx = End.to_dp_index(phmm);
                    phmm.exit_core_from_end(visitor, layer_idx, phmm.end().semilocal_params().get_score(layer_idx))?;
                    break;
                }
                EndInsertExit::Insert => {
                    let emit_params = &layer.emission_insert;
                    state = PhmmState::Insert;

                    phmm.choose_emission(visitor, layer_idx, state, emit_params, phmm.mapping())?;
                    num_emitted += 1;
                }
                EndInsertExit::Exit => break,
            }
        } else {
            let next_state = phmm.choose_end_or_insert(visitor, layer_idx, state, &layer.transition)?;
            match next_state {
                EndInsert::End => {
                    layer_idx = End.to_dp_index(phmm);
                    phmm.exit_core_from_end(visitor, layer_idx, phmm.end().semilocal_params().get_score(layer_idx))?;
                    break;
                }
                EndInsert::Insert => {
                    state = PhmmState::Insert;

                    phmm.choose_emission(visitor, layer_idx, state, &layer.emission_insert, phmm.mapping())?;
                    num_emitted += 1;
                }
            }
        }
    }

    Ok(TraverseCorePhmmOrExitOutput {
        exit_layer:     layer_idx,
        exit_param:     phmm.end().semilocal_params().get_score(layer_idx),
        aligned_layers: enter_layer..=layer_idx,
        query_range:    SeqIndex(num_emitted_begin_module)..SeqIndex(num_emitted),
    })
}

/// Traverses the [`DomainModule`] using a given visitor, for use with
/// [`DomainPhmm`] and [`LocalPhmm`].
///
/// ## Errors
///
/// See the visitor's implementations for [`choose_domain_emission`],
/// [`enter_module_insert`], and [`exit_module_insert`].
///
/// [`choose_domain_emission`]: VisitDomainModule::choose_domain_emission
/// [`enter_module_insert`]: VisitDomainModule::enter_module_insert
/// [`exit_module_insert`]: VisitDomainModule::exit_module_insert
/// [`DomainModule`]: crate::alignment::phmm::modules::DomainModule
fn traverse_domain_module<V, P, T, const S: usize>(
    phmm: &P, loc: ModuleLocation, visitor: &mut V,
) -> Result<TraverseDomainModuleOutput, P::Error>
where
    P: GetMapping<S> + GetModule<Begin: DomainParams<T, S>, End: DomainParams<T, S>> + VisitDomainModule<V, T, S>, {
    let module = match loc {
        ModuleLocation::Begin => phmm.begin().domain_params(),
        ModuleLocation::End => phmm.end().domain_params(),
    };

    let mut num_emitted = 0;

    if !phmm.enter_module_insert(visitor, module, loc)? {
        phmm.exiting_domain_module(visitor, module, loc)?;
        return Ok(TraverseDomainModuleOutput { num_emitted });
    }

    loop {
        phmm.choose_domain_emission(visitor, &module.background_emission, phmm.mapping(), loc)?;
        num_emitted += 1;

        if phmm.exit_module_insert(visitor, module, loc)? {
            break;
        }
    }

    phmm.exiting_domain_module(visitor, module, loc)?;

    Ok(TraverseDomainModuleOutput { num_emitted })
}

struct TraverseDomainModuleOutput {
    num_emitted: usize,
}

impl<T, const S: usize> GlobalPhmm<T, S> {
    /// Performs full traversal of the [`GlobalPhmm`], making decisions using
    /// the provided `visitor`.
    ///
    /// ## Errors
    ///
    /// See the visitor's implementations for [`choose_emission`],
    /// [`choose_core_transition`], [`choose_end_or_insert`], and [`finalize`].
    ///
    /// [`choose_emission`]: GlobalVisitor::choose_emission
    /// [`choose_core_transition`]: GlobalVisitor::choose_core_transition
    /// [`choose_end_or_insert`]: GlobalVisitor::choose_end_or_insert
    /// [`finalize`]: GlobalVisitor::finalize
    pub fn traverse<V>(&self, mut visitor: V) -> Result<V::Output, V::Error>
    where
        V: GlobalVisitor<T, S>, {
        traverse_core_phmm(self, 0, &mut visitor)?;
        visitor.finalize(self)
    }
}

impl<T, const S: usize> DomainPhmm<T, S> {
    /// Performs full traversal of the [`DomainPhmm`], making decisions using
    /// the provided `visitor`.
    ///
    /// ## Errors
    ///
    /// See the visitor's implementations for [`choose_domain_emission`],
    /// [`enter_module_insert`], [`exit_module_insert`], [`choose_emission`],
    /// [`choose_core_transition`], [`choose_end_or_insert`], and [`finalize`].
    ///
    /// [`choose_domain_emission`]: DomainVisitor::choose_domain_emission
    /// [`enter_module_insert`]: DomainVisitor::enter_module_insert
    /// [`exit_module_insert`]: DomainVisitor::exit_module_insert
    /// [`choose_emission`]: DomainVisitor::choose_emission
    /// [`choose_core_transition`]: DomainVisitor::choose_core_transition
    /// [`choose_end_or_insert`]: DomainVisitor::choose_end_or_insert
    /// [`finalize`]: DomainVisitor::finalize
    pub fn traverse<V>(&self, mut visitor: V) -> Result<V::Output, V::Error>
    where
        V: DomainVisitor<T, S>, {
        let TraverseDomainModuleOutput { num_emitted } = traverse_domain_module(self, ModuleLocation::Begin, &mut visitor)?;

        let query_start = SeqIndex(num_emitted);

        let num_emitted = traverse_core_phmm(self, num_emitted, &mut visitor)?;

        let query_range = query_start..SeqIndex(num_emitted);
        traverse_domain_module(self, ModuleLocation::End, &mut visitor)?;

        visitor.finalize(self, query_range.clone())
    }
}

impl<T, const S: usize> SemiLocalPhmm<T, S> {
    /// Performs full traversal of the [`SemiLocalPhmm`], making decisions using
    /// the provided `visitor`.
    ///
    /// ## Errors
    ///
    /// See the visitor's implementations for [`enter_core`],
    /// [`choose_emission`], [`choose_core_transition`],
    /// [`choose_core_transition_or_exit`], [`choose_end_or_insert`],
    /// [`choose_end_insert_or_exit`], [`exit_core_from_end`], and [`finalize`].
    ///
    /// [`enter_core`]: SemiLocalVisitor::enter_core
    /// [`choose_emission`]: SemiLocalVisitor::choose_emission
    /// [`choose_core_transition`]: SemiLocalVisitor::choose_core_transition
    /// [`choose_core_transition_or_exit`]:
    ///     SemiLocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_or_insert`]: SemiLocalVisitor::choose_end_or_insert
    /// [`choose_end_insert_or_exit`]:
    ///     SemiLocalVisitor::choose_end_insert_or_exit
    /// [`exit_core_from_end`]: SemiLocalVisitor::exit_core_from_end
    /// [`finalize`]: SemiLocalVisitor::finalize
    pub fn traverse<V>(&self, mut visitor: V) -> Result<V::Output, V::Error>
    where
        V: SemiLocalVisitor<T, S>,
        T: PhmmNumber + 'static, {
        let layer = visitor.enter_core(self.begin().semilocal_params(), self)?;

        let TraverseCorePhmmOrExitOutput {
            exit_layer,
            exit_param,
            aligned_layers,
            ..
        } = traverse_core_phmm_or_exit(self, layer, 0, &mut visitor)?;

        visitor.exit_core(exit_layer, exit_param, self)?;
        visitor.finalize(self, aligned_layers)
    }
}

impl<T, const S: usize> LocalPhmm<T, S> {
    /// Performs full traversal of the [`LocalPhmm`], making decisions using the
    /// provided `visitor`.
    ///
    /// ## Errors
    ///
    /// See the visitor's implementations for [`choose_local_emission`],
    /// [`enter_module_insert`], [`exit_module_insert`], [`enter_core`],
    /// [`choose_emission`], [`choose_core_transition`],
    /// [`choose_core_transition_or_exit`], [`choose_end_or_insert`],
    /// [`choose_end_insert_or_exit`], [`exit_core_from_end`], and [`finalize`].
    ///
    /// [`choose_local_emission`]: LocalVisitor::choose_local_emission
    /// [`enter_module_insert`]: LocalVisitor::enter_module_insert
    /// [`exit_module_insert`]: LocalVisitor::exit_module_insert
    /// [`enter_core`]: LocalVisitor::enter_core
    /// [`choose_emission`]: LocalVisitor::choose_emission
    /// [`choose_core_transition`]: LocalVisitor::choose_core_transition
    /// [`choose_core_transition_or_exit`]:
    ///     LocalVisitor::choose_core_transition_or_exit
    /// [`choose_end_or_insert`]: LocalVisitor::choose_end_or_insert
    /// [`choose_end_insert_or_exit`]: LocalVisitor::choose_end_insert_or_exit
    /// [`exit_core_from_end`]: LocalVisitor::exit_core_from_end
    /// [`finalize`]: LocalVisitor::finalize
    pub fn traverse<V>(&self, mut visitor: V) -> Result<V::Output, V::Error>
    where
        V: LocalVisitor<T, S>,
        T: PhmmNumber + 'static, {
        let TraverseDomainModuleOutput {
            num_emitted: num_emitted_begin,
        } = traverse_domain_module(self, ModuleLocation::Begin, &mut visitor)?;

        let enter_layer = visitor.enter_core(&self.begin().semilocal_params, self)?;

        let TraverseCorePhmmOrExitOutput {
            exit_layer,
            exit_param,
            aligned_layers,
            query_range,
        } = traverse_core_phmm_or_exit(self, enter_layer, num_emitted_begin, &mut visitor)?;

        visitor.exit_core(exit_layer, exit_param, self)?;

        traverse_domain_module(self, ModuleLocation::End, &mut visitor)?;

        visitor.finalize(self, aligned_layers, query_range)
    }
}
