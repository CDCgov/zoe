//! The drivers for use in traversal code.

use crate::{
    alignment::phmm::{
        DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
        components::{EmissionParams, LayerParams},
        indexing::{Begin, DpIndex, End, GetLayer, GetMapping, GetModule, PhmmIndex, PhmmIndexable, SeqIndex},
        modules::{DomainModule, DomainParams, SemiLocalParams},
        state::{DomainModuleState, PhmmState},
        traverse::{
            DomainVisitor, EndInsert, EndInsertExit, GlobalVisitor, LocalVisitor, ModuleLocation, SemiLocalVisitor,
            VisitCore, VisitCoreOrExit, VisitDomainModule,
        },
        views::{
            AsGlobalView, DomainPhmmView, GetLayerView, GetModuleView, GlobalPhmmView, LocalPhmmView, SemiLocalPhmmView,
        },
    },
    data::views::AsView,
};
use std::ops::{Range, RangeInclusive};

/// A driver for managing traversal of the [`CorePhmm`] within a [`GlobalPhmm`]
/// or [`DomainPhmm`] (no early exit is permitted).
///
/// ## Parameters
///
/// - `'a`: The lifetime of the pHMM view
/// - `P` - The view type of the pHMM being traversed (either [`GlobalPhmmView`]
///   or [`DomainPhmmView`])
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
///
/// [`CorePhmm`]: crate::alignment::phmm::models::CorePhmm
struct CorePhmmDriver<'a, P, T, const S: usize> {
    /// A view of the pHMM whose core will be traversed.
    phmm:             P,
    /// The current layer index of the traversal within `phmm`.
    layer_idx:        DpIndex,
    /// The current layer at `layer_idx`.
    layer:            &'a LayerParams<T, S>,
    /// The layers after `layer_idx`.
    remaining_layers: &'a [LayerParams<T, S>],
    /// The current state of the traversal within `phmm`.
    state:            PhmmState,
    /// If `Some`, then the next step is to emit a residue from the specified
    /// parameters.
    emit:             Option<&'a EmissionParams<T, S>>,
    /// The number of emitted residues from the pHMM.
    num_emitted:      usize,
}

impl<'a, P, T, const S: usize> CorePhmmDriver<'a, P, T, S>
where
    P: PhmmIndexable + GetMapping<S> + GetLayerView<'a, T, S> + AsGlobalView<T, S>,
{
    /// Constructs a new [`CorePhmmDriver`] from the specified pHMM view and the
    /// current number of emitted residues (in the case of a [`DomainModule`]).
    ///
    /// The traversal begins at the BEGIN state.
    fn new(phmm: P, num_emitted: usize) -> Self {
        let (layer, remaining_layers) = phmm.split_first_layer_view();

        Self {
            layer_idx: phmm.to_dp_index(Begin),
            layer,
            remaining_layers,
            phmm,
            state: PhmmState::Match,
            emit: None,
            num_emitted,
        }
    }

    /// Advances the [`CorePhmmDriver`] by one step using the provided `visitor`
    /// to select between options. Returns `true` if there is still more of the
    /// core pHMM to traverse (e.g., for use in a while loop).
    ///
    /// ## Errors
    ///
    /// See the visitor's implementations for [`choose_emission`],
    /// [`choose_core_transition`], and [`choose_end_or_insert`].
    ///
    /// [`choose_emission`]: VisitCore::choose_emission
    /// [`choose_core_transition`]: VisitCore::choose_core_transition
    /// [`choose_end_or_insert`]: VisitCore::choose_end_or_insert
    fn advance<V>(&mut self, visitor: &mut V) -> Result<bool, P::Error>
    where
        P: VisitCore<V, T, S> + GetLayer<T, S> + Copy, {
        if let Some(params) = self.emit {
            self.phmm
                .choose_emission(visitor, self.layer_idx, self.state, params, self.phmm.mapping())?;
            self.emit = None;
            self.num_emitted += 1;
        } else if let Some((next_layer, rest)) = self.remaining_layers.split_first() {
            let params = &self.layer.transition;
            let next_state = self
                .phmm
                .choose_core_transition(visitor, self.layer_idx, self.state, params)?;
            self.move_to_state(next_state, next_layer, rest);
        } else {
            // The current layer is LastMatch
            let params = &self.layer.transition;
            let next_state = self.phmm.choose_end_or_insert(visitor, self.layer_idx, self.state, params)?;
            match next_state {
                EndInsert::End => return Ok(false),
                EndInsert::Insert => self.move_to_insert(),
            }
        }

        Ok(true)
    }

    /// Transitions the driver to the specified state.
    ///
    /// This sets `emit` if the state emits a residue, and it also increments
    /// `layer_idx` (and updates `layer` and `remaining_layers`) if the next
    /// state is not [`PhmmState::Insert`].
    fn move_to_state(&mut self, state: PhmmState, next_layer: &'a LayerParams<T, S>, rest: &'a [LayerParams<T, S>]) {
        match state {
            PhmmState::Match => {
                self.emit = Some(&self.layer.emission_match);
                self.layer_idx = self.layer_idx.next_index(&self.phmm);
                self.layer = next_layer;
                self.remaining_layers = rest;
            }
            PhmmState::Delete => {
                self.layer_idx = self.layer_idx.next_index(&self.phmm);
                self.layer = next_layer;
                self.remaining_layers = rest;
            }
            PhmmState::Insert => return self.move_to_insert(),
        }
        self.state = state;
    }

    /// Transitions the driver to the insert state.
    ///
    /// This does not increment `layer_idx`.
    fn move_to_insert(&mut self) {
        self.emit = Some(&self.layer.emission_insert);
        self.state = PhmmState::Insert;
    }
}

/// A driver for managing traversal of the [`CorePhmm`] within a
/// [`SemiLocalPhmm`] or [`LocalPhmm`] (early exit is permitted).
///
/// ## Parameters
///
/// - `'a`: The lifetime of the pHMM view
/// - `P` - The view type of the pHMM being traversed (either
///   [`SemiLocalPhmmView`] or [`LocalPhmmView`])
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
///
/// [`CorePhmm`]: crate::alignment::phmm::models::CorePhmm
struct CorePhmmOrExitDriver<'a, P, T, const S: usize> {
    /// A view of the pHMM whose core will be traversed.
    phmm:             P,
    /// The current layer of the traversal within `phmm`.
    ///
    /// Once [`advance`] returns false, this means the core pHMM has been
    /// exited, and this remains equal to the layer index that was exited from.
    ///
    /// [`advance`]: CorePhmmOrExitDriver::advance
    layer_idx:        DpIndex,
    /// The current layer at `layer_idx`. This is `None` if `layer_idx` is
    /// [`End`].
    layer:            Option<&'a LayerParams<T, S>>,
    /// The layers after `layer_idx`.
    remaining_layers: &'a [LayerParams<T, S>],
    /// The current state of the traversal within `phmm`.
    state:            PhmmState,
    /// If `Some`, then the next step is to emit a residue from the specified
    /// parameters.
    emit:             Option<&'a EmissionParams<T, S>>,
    /// The number of emitted residues from the pHMM.
    num_emitted:      usize,
}

impl<'a, P, T, const S: usize> CorePhmmOrExitDriver<'a, P, T, S>
where
    P: PhmmIndexable + GetLayerView<'a, T, S> + GetMapping<S> + GetModule<End: SemiLocalParams<T>>,
    T: PhmmNumber,
{
    /// Constructs a new [`CorePhmmOrExitDriver`] from the specified pHMM view
    /// and the current number of emitted residues (in the case of a
    /// [`DomainModule`]). The traversal begins at specified `layer`.
    ///
    /// ## Panics
    ///
    /// If the `layer_idx` is out of bounds, this method panics.
    fn new(phmm: P, layer_idx: DpIndex, num_emitted: usize) -> Self {
        if layer_idx.eq_index(End, &phmm) {
            return Self {
                phmm,
                layer_idx,
                layer: None,
                remaining_layers: &[],
                state: PhmmState::Match,
                emit: None,
                num_emitted,
            };
        }

        // Get the current layer (since we are not at End), any remaining
        // layers, and optionally the previous layer if we are not at Begin.
        let (prev_layer, layer, remaining_layers) = if let Some((before, layer_and_after)) =
            phmm.split_layers_at_view(layer_idx)
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
        let emit = prev_layer.map(|layer| &layer.emission_match);

        Self {
            phmm,
            layer_idx,
            layer: Some(layer),
            remaining_layers,
            state: PhmmState::Match,
            emit,
            num_emitted,
        }
    }

    /// Advances the [`CorePhmmOrExitDriver`] by one step using the provided
    /// `visitor` to select between options. Returns `true` if there is still
    /// more of the core pHMM to traverse (e.g., for use in a while loop).
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
    /// [`choose_end_insert_or_exit`]:
    ///     VisitCoreOrExit::choose_end_insert_or_exit
    /// [`exit_core_from_end`]: VisitCoreOrExit::exit_core_from_end
    fn advance<V>(&mut self, visitor: &mut V) -> Result<bool, P::Error>
    where
        P: VisitCoreOrExit<V, T, S> + Copy, {
        if let Some(params) = self.emit {
            self.phmm
                .choose_emission(visitor, self.layer_idx, self.state, params, self.phmm.mapping())?;
            self.emit = None;
            self.num_emitted += 1;
            return Ok(true);
        }

        let Some(layer) = self.layer else {
            // The current layer is End
            self.phmm.exit_core_from_end(
                visitor,
                self.layer_idx,
                self.phmm.end().semilocal_params().get_score(self.layer_idx),
            )?;
            return Ok(false);
        };

        if let Some((next_layer, rest)) = self.remaining_layers.split_first() {
            let params = &layer.transition;

            if self.state == PhmmState::Match {
                let exit_param = self.phmm.end().semilocal_params().get_score(self.layer_idx);
                let next_state = self
                    .phmm
                    .choose_core_transition_or_exit(visitor, self.layer_idx, params, exit_param)?;
                if let Some(state) = PhmmState::get_from(next_state) {
                    self.move_to_state(state, layer, next_layer, rest);
                } else {
                    return Ok(false);
                }
            } else {
                let next_state = self
                    .phmm
                    .choose_core_transition(visitor, self.layer_idx, self.state, params)?;
                self.move_to_state(next_state, layer, next_layer, rest);
            }
        } else {
            // The current layer is LastMatch

            let params = &layer.transition;

            if self.state == PhmmState::Match {
                let exit_param = self.phmm.end().semilocal_params().get_score(self.layer_idx);
                let exit_from_end_param = self.phmm.end().semilocal_params().get_score(End);
                let next_state =
                    self.phmm
                        .choose_end_insert_or_exit(visitor, self.layer_idx, params, exit_param, exit_from_end_param)?;
                match next_state {
                    EndInsertExit::End => {
                        self.layer_idx = self.phmm.to_dp_index(End);
                        self.layer = None;
                    }
                    EndInsertExit::Insert => self.move_to_insert(layer),
                    EndInsertExit::Exit => return Ok(false),
                }
            } else {
                let next_state = self.phmm.choose_end_or_insert(visitor, self.layer_idx, self.state, params)?;
                match next_state {
                    EndInsert::End => {
                        self.layer_idx = self.phmm.to_dp_index(End);
                        self.layer = None;
                    }
                    EndInsert::Insert => self.move_to_insert(layer),
                }
            }
        }

        Ok(true)
    }

    /// Transitions the driver to the specified state. This sets `emit` if the
    /// state emits a residue, and it also increments `layer` if the next state
    /// is not [`PhmmState::Insert`].
    fn move_to_state(
        &mut self, state: PhmmState, current_layer: &'a LayerParams<T, S>, next_layer: &'a LayerParams<T, S>,
        rest: &'a [LayerParams<T, S>],
    ) {
        match state {
            PhmmState::Match => {
                self.emit = Some(&current_layer.emission_match);
                self.layer_idx = self.layer_idx.next_index(&self.phmm);
                self.layer = Some(next_layer);
                self.remaining_layers = rest;
            }
            PhmmState::Delete => {
                self.layer_idx = self.layer_idx.next_index(&self.phmm);
                self.layer = Some(next_layer);
                self.remaining_layers = rest;
            }
            PhmmState::Insert => return self.move_to_insert(current_layer),
        }
        self.state = state;
    }

    /// Transitions the driver to the insert state.
    ///
    /// This does not increment `layer_idx`.
    fn move_to_insert(&mut self, current_layer: &'a LayerParams<T, S>) {
        self.emit = Some(&current_layer.emission_insert);
        self.state = PhmmState::Insert;
    }
}

/// A driver for managing traversal of a [`DomainModule`] on either end of a
/// [`DomainPhmm`] or [`LocalPhmm`] (since [`DomainModule`] is used internally
/// by [`LocalModule`]).
///
/// ## Parameters
///
/// - `'a`: The lifetime of the pHMM view
/// - `P` - The view type of the pHMM being traversed (either [`DomainPhmmView`]
///   or [`LocalPhmmView`])
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
///
/// [`LocalModule`]: crate::alignment::phmm::modules::LocalModule
struct DomainModuleDriver<'a, P, T, const S: usize> {
    /// The module being traversed.
    module:      &'a DomainModule<T, S>,
    /// A view of the full pHMM.
    phmm:        P,
    /// The current state of the driver within the module.
    state:       DomainModuleState,
    /// The location of the module within the pHMM.
    loc:         ModuleLocation,
    /// The number of residues that have been emitted from the module so far.
    num_emitted: usize,
    /// Whether the next step of the driver is to emit a residue from the
    /// module.
    emit:        bool,
}

impl<'a, P, T, const S: usize> DomainModuleDriver<'a, P, T, S>
where
    P: GetMapping<S> + GetModuleView<'a, Begin: DomainParams<T, S>, End: DomainParams<T, S>> + 'a,
{
    /// Constructs a new [`DomainModuleDriver`] from the specified pHMM view and
    /// module location. The traversal begins at the beginning of the selected
    /// module.
    fn new(phmm: P, loc: ModuleLocation) -> Self {
        let module = match loc {
            ModuleLocation::Begin => phmm.begin_view().domain_params(),
            ModuleLocation::End => phmm.end_view().domain_params(),
        };

        Self {
            module,
            phmm,
            state: DomainModuleState::Begin,
            loc,
            num_emitted: 0,
            emit: false,
        }
    }

    // TODO: We actually may be doing two steps. Can this be reworked more
    // generally?
    /// Advances the [`DomainModuleDriver`] by one step using the provided
    /// `visitor` to select between options. Returns `true` if there is still
    /// more of the module to traverse (e.g., for use in a while loop).
    ///
    /// This returns `false` when [`enter_module_insert`] returns `false` or
    /// [`exit_module_insert`] returns `true`.
    ///
    /// ## Errors
    ///
    /// See the visitor's implementations for [`choose_domain_emission`],
    /// [`enter_module_insert`], and [`exit_module_insert`].
    ///
    /// [`choose_domain_emission`]: VisitDomainModule::choose_domain_emission
    /// [`enter_module_insert`]: VisitDomainModule::enter_module_insert
    /// [`exit_module_insert`]: VisitDomainModule::exit_module_insert
    fn advance<V>(&mut self, visitor: &mut V) -> Result<bool, P::Error>
    where
        P: VisitDomainModule<V, T, S> + Copy, {
        if self.emit {
            self.phmm
                .choose_domain_emission(visitor, &self.module.background_emission, self.phmm.mapping(), self.loc)?;
            self.num_emitted += 1;
            self.emit = false;
            Ok(true)
        } else {
            match self.state {
                DomainModuleState::Begin => {
                    if self.phmm.enter_module_insert(visitor, self.module, self.loc)? {
                        self.state = DomainModuleState::Insert;
                        self.emit = true;
                        Ok(true)
                    } else {
                        self.state = DomainModuleState::End;
                        self.phmm.exiting_domain_module(visitor, self.module, self.loc)?;
                        Ok(false)
                    }
                }
                DomainModuleState::Insert => {
                    if self.phmm.exit_module_insert(visitor, self.module, self.loc)? {
                        self.state = DomainModuleState::End;
                        self.phmm.exiting_domain_module(visitor, self.module, self.loc)?;
                        Ok(false)
                    } else {
                        self.emit = true;
                        Ok(true)
                    }
                }
                DomainModuleState::End => Ok(false),
            }
        }
    }
}

/// The internal state of the [`GlobalPhmmDriver`] (which internally is a state
/// machine).
enum GlobalPhmmDriverState<'a, T, const S: usize> {
    /// Traversal is currently within the core pHMM.
    Core(CorePhmmDriver<'a, GlobalPhmmView<'a, T, S>, T, S>),
    /// Traversal has finished, but a final call to [`finalize`] still needs to
    /// occur.
    ///
    /// [`finalize`]: GlobalVisitor::finalize
    Finalize(GlobalPhmmView<'a, T, S>),
}

/// A driver for managing traversal of a [`GlobalPhmm`].
///
/// ## Parameters
///
/// - `'a`: The lifetime of the pHMM reference
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
pub struct GlobalPhmmDriver<'a, T, const S: usize> {
    /// The internal state of the driver (viewed as a state machine).
    state: GlobalPhmmDriverState<'a, T, S>,
}

impl<'a, T, const S: usize> GlobalPhmmDriver<'a, T, S>
where
    T: Clone + 'static,
{
    /// Initializes a new [`GlobalPhmmDriver`] for traversing the specified
    /// pHMM, starting from the BEGIN state.
    #[must_use]
    pub fn new(phmm: &'a GlobalPhmm<T, S>) -> Self {
        Self {
            state: GlobalPhmmDriverState::Core(CorePhmmDriver::new(phmm.as_view(), 0)),
        }
    }

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
    pub fn run<V>(mut self, mut visitor: V) -> Result<V::Output, V::Error>
    where
        V: GlobalVisitor<T, S>, {
        loop {
            match &mut self.state {
                GlobalPhmmDriverState::Core(driver) => {
                    if !driver.advance(&mut visitor)? {
                        self.state = GlobalPhmmDriverState::Finalize(driver.phmm);
                    }
                }
                GlobalPhmmDriverState::Finalize(phmm) => {
                    return visitor.finalize(*phmm);
                }
            }
        }
    }
}

/// The internal state of the [`DomainPhmmDriver`] (which internally is a state
/// machine).
enum DomainPhmmDriverState<'a, T, const S: usize> {
    /// Traversal is current within the [`DomainModule`] at the start of the
    /// pHMM.
    StartModule(DomainModuleDriver<'a, DomainPhmmView<'a, T, S>, T, S>),
    /// Traversal is currently within the core pHMM.
    Core {
        driver:      CorePhmmDriver<'a, DomainPhmmView<'a, T, S>, T, S>,
        query_start: SeqIndex,
    },
    /// Traversal is current within the [`DomainModule`] at the end of the pHMM.
    EndModule {
        driver:      DomainModuleDriver<'a, DomainPhmmView<'a, T, S>, T, S>,
        query_range: Range<SeqIndex>,
    },
    /// Traversal has finished, but a final call to [`finalize`] still needs to
    /// occur.
    ///
    /// [`finalize`]: DomainVisitor::finalize
    Finalize {
        phmm:        DomainPhmmView<'a, T, S>,
        query_range: Range<SeqIndex>,
    },
}

/// A driver for managing traversal of a [`DomainPhmm`].
///
/// ## Parameters
///
/// - `'a`: The lifetime of the pHMM reference
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
pub struct DomainPhmmDriver<'a, T, const S: usize> {
    /// The internal state of the driver (viewed as a state machine).
    state: DomainPhmmDriverState<'a, T, S>,
}

impl<'a, T, const S: usize> DomainPhmmDriver<'a, T, S>
where
    T: Copy + 'static,
{
    /// Initializes a new [`DomainPhmmDriver`] for traversing the specified
    /// pHMM, starting from the domain module at the beginning of the pHMM.
    #[must_use]
    pub fn new(phmm: &'a DomainPhmm<T, S>) -> Self {
        Self {
            state: DomainPhmmDriverState::StartModule(DomainModuleDriver::new(phmm.as_view(), ModuleLocation::Begin)),
        }
    }

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
    pub fn run<V>(mut self, mut visitor: V) -> Result<V::Output, V::Error>
    where
        V: DomainVisitor<T, S>, {
        loop {
            match &mut self.state {
                DomainPhmmDriverState::StartModule(driver) => {
                    if !driver.advance(&mut visitor)? {
                        self.state = DomainPhmmDriverState::Core {
                            driver:      CorePhmmDriver::new(driver.phmm, driver.num_emitted),
                            query_start: SeqIndex(driver.num_emitted),
                        };
                    }
                }
                DomainPhmmDriverState::Core { driver, query_start } => {
                    if !driver.advance(&mut visitor)? {
                        self.state = DomainPhmmDriverState::EndModule {
                            driver:      DomainModuleDriver::new(driver.phmm, ModuleLocation::End),
                            query_range: *query_start..SeqIndex(driver.num_emitted),
                        };
                    }
                }
                DomainPhmmDriverState::EndModule { driver, query_range } => {
                    if !driver.advance(&mut visitor)? {
                        self.state = DomainPhmmDriverState::Finalize {
                            phmm:        driver.phmm,
                            query_range: query_range.clone(),
                        };
                    }
                }
                DomainPhmmDriverState::Finalize { phmm, query_range } => {
                    return visitor.finalize(*phmm, query_range.clone());
                }
            }
        }
    }
}

/// The internal state of the [`SemiLocalPhmmDriver`] (which internally is a
/// state machine).
enum SemiLocalPhmmDriverState<'a, T, const S: usize> {
    /// Traversal is current within the [`SemiLocalModule`] at the start of the
    /// pHMM (the next step will be picking a layer of the core pHMM to enter).
    ///
    /// [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule
    StartModule,
    /// Traversal is currently within the core pHMM.
    Core {
        driver:        CorePhmmOrExitDriver<'a, SemiLocalPhmmView<'a, T, S>, T, S>,
        entered_layer: DpIndex,
    },
    /// Traversal is exiting the pHMM, so call the exit hook.
    Exiting {
        phmm:             SemiLocalPhmmView<'a, T, S>,
        exit_from:        DpIndex,
        exit_param:       T,
        traversed_layers: RangeInclusive<DpIndex>,
    },
    /// Traversal has finished, but a final call to [`finalize`] still needs to
    /// occur.
    ///
    /// [`finalize`]: SemiLocalVisitor::finalize
    Finalize {
        phmm:             SemiLocalPhmmView<'a, T, S>,
        traversed_layers: RangeInclusive<DpIndex>,
    },
}

/// A driver for managing traversal of a [`SemiLocalPhmm`].
///
/// ## Parameters
///
/// - `'a`: The lifetime of the pHMM reference
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
pub struct SemiLocalPhmmDriver<'a, T, const S: usize> {
    /// A reference to the pHMM being traversed (since the first state of
    /// [`SemiLocalPhmmDriverState`] does not hold it, but subsequent states
    /// require it).
    phmm:  &'a SemiLocalPhmm<T, S>,
    /// The internal state of the driver (viewed as a state machine).
    state: SemiLocalPhmmDriverState<'a, T, S>,
}

impl<'a, T, const S: usize> SemiLocalPhmmDriver<'a, T, S>
where
    T: PhmmNumber + 'static,
{
    /// Initializes a new [`SemiLocalPhmmDriver`] for traversing the specified
    /// pHMM, starting from the semilocal module at the beginning of the pHMM.
    #[must_use]
    pub fn new(phmm: &'a SemiLocalPhmm<T, S>) -> Self {
        Self {
            phmm,
            state: SemiLocalPhmmDriverState::StartModule,
        }
    }

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
    pub fn run<V>(mut self, mut visitor: V) -> Result<V::Output, V::Error>
    where
        V: SemiLocalVisitor<T, S>, {
        loop {
            match &mut self.state {
                SemiLocalPhmmDriverState::StartModule => {
                    let layer = visitor.enter_core(self.phmm.begin().semilocal_params(), self.phmm.as_view())?;
                    self.state = SemiLocalPhmmDriverState::Core {
                        driver:        CorePhmmOrExitDriver::new(self.phmm.as_view(), layer, 0),
                        entered_layer: layer,
                    };
                }
                SemiLocalPhmmDriverState::Core { driver, entered_layer } => {
                    if !driver.advance(&mut visitor)? {
                        self.state = SemiLocalPhmmDriverState::Exiting {
                            phmm:             driver.phmm,
                            exit_from:        driver.layer_idx,
                            exit_param:       driver.phmm.end().semilocal_params().get_score(driver.layer_idx),
                            traversed_layers: *entered_layer..=driver.layer_idx,
                        };
                    }
                }
                SemiLocalPhmmDriverState::Exiting {
                    phmm,
                    exit_from,
                    exit_param,
                    traversed_layers,
                } => {
                    visitor.exit_core(*exit_from, *exit_param, *phmm)?;
                    self.state = SemiLocalPhmmDriverState::Finalize {
                        phmm:             *phmm,
                        traversed_layers: traversed_layers.clone(),
                    };
                }
                SemiLocalPhmmDriverState::Finalize { phmm, traversed_layers } => {
                    // TODO: layers_traversed or traversed_layers
                    return visitor.finalize(*phmm, traversed_layers.clone());
                }
            }
        }
    }
}

// TODO: Naming of states doesn't match semilocal

/// The internal state of the [`LocalPhmmDriver`] (which internally is a state
/// machine).
enum LocalPhmmDriverState<'a, T, const S: usize> {
    /// Traversal is current within the [`DomainModule`] at the start of the
    /// pHMM.
    StartModule(DomainModuleDriver<'a, LocalPhmmView<'a, T, S>, T, S>),
    /// Traversal is current within the [`SemiLocalModule`] at the start of the
    /// pHMM (the next step will be picking a layer of the core pHMM to enter).
    ///
    /// [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule
    ToCore { num_emitted_begin: usize },
    /// Traversal is currently within the core pHMM.
    Core {
        driver:      CorePhmmOrExitDriver<'a, LocalPhmmView<'a, T, S>, T, S>,
        ref_start:   DpIndex,
        query_start: SeqIndex,
    },
    /// Traversal is exiting the pHMM, so call the exit hook.
    Exiting {
        phmm:        LocalPhmmView<'a, T, S>,
        exit_from:   DpIndex,
        exit_param:  T,
        ref_range:   Range<DpIndex>,
        query_range: Range<SeqIndex>,
    },
    /// Traversal is current within the [`DomainModule`] at the end of the pHMM.
    EndModule {
        driver:      DomainModuleDriver<'a, LocalPhmmView<'a, T, S>, T, S>,
        ref_range:   Range<DpIndex>,
        query_range: Range<SeqIndex>,
    },
    /// Traversal has finished, but a final call to [`finalize`] still needs to
    /// occur.
    ///
    /// [`finalize`]: LocalVisitor::finalize
    Finalize {
        ref_range:   Range<DpIndex>,
        query_range: Range<SeqIndex>,
    },
}

/// A driver for managing traversal of a [`LocalPhmm`].
///
/// ## Parameters
///
/// - `'a`: The lifetime of the pHMM reference
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
pub struct LocalPhmmDriver<'a, T, const S: usize> {
    /// A reference to the pHMM being traversed (since
    /// [`LocalPhmmDriverState::ToCore`] does not hold it, but subsequent states
    /// require it).
    phmm:  &'a LocalPhmm<T, S>,
    /// The internal state of the driver (viewed as a state machine).
    state: LocalPhmmDriverState<'a, T, S>,
}

impl<'a, T, const S: usize> LocalPhmmDriver<'a, T, S>
where
    T: PhmmNumber + 'static,
{
    /// Initializes a new [`LocalPhmmDriver`] for traversing the specified pHMM,
    /// starting from the local module at the beginning of the pHMM.
    #[must_use]
    pub fn new(phmm: &'a LocalPhmm<T, S>) -> Self {
        Self {
            phmm,
            state: LocalPhmmDriverState::StartModule(DomainModuleDriver::new(phmm.as_view(), ModuleLocation::Begin)),
        }
    }

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
    pub fn run<V>(mut self, mut visitor: V) -> Result<V::Output, V::Error>
    where
        V: LocalVisitor<T, S>, {
        loop {
            match &mut self.state {
                LocalPhmmDriverState::StartModule(driver) => {
                    if !driver.advance(&mut visitor)? {
                        self.state = LocalPhmmDriverState::ToCore {
                            num_emitted_begin: driver.num_emitted,
                        };
                    }
                }
                LocalPhmmDriverState::ToCore { num_emitted_begin } => {
                    let layer = visitor.enter_core(&self.phmm.begin().semilocal_params, self.phmm.as_view())?;
                    self.state = LocalPhmmDriverState::Core {
                        driver:      CorePhmmOrExitDriver::new(self.phmm.as_view(), layer, *num_emitted_begin),
                        ref_start:   layer,
                        query_start: SeqIndex(*num_emitted_begin),
                    };
                }
                LocalPhmmDriverState::Core {
                    driver,
                    ref_start,
                    query_start,
                } => {
                    if !driver.advance(&mut visitor)? {
                        self.state = LocalPhmmDriverState::Exiting {
                            phmm:        driver.phmm,
                            exit_from:   driver.layer_idx,
                            exit_param:  driver.phmm.end().semilocal_params().get_score(driver.layer_idx),
                            ref_range:   *ref_start..driver.layer_idx.next_index(&driver.phmm),
                            query_range: *query_start..SeqIndex(driver.num_emitted),
                        };
                    }
                }
                LocalPhmmDriverState::Exiting {
                    phmm,
                    exit_from,
                    exit_param,
                    ref_range,
                    query_range,
                } => {
                    visitor.exit_core(*exit_from, *exit_param, *phmm)?;
                    self.state = LocalPhmmDriverState::EndModule {
                        driver:      DomainModuleDriver::new(self.phmm.as_view(), ModuleLocation::End),
                        ref_range:   ref_range.clone(),
                        query_range: query_range.clone(),
                    };
                }
                LocalPhmmDriverState::EndModule {
                    driver,
                    ref_range,
                    query_range,
                } => {
                    if !driver.advance(&mut visitor)? {
                        self.state = LocalPhmmDriverState::Finalize {
                            ref_range:   ref_range.clone(),
                            query_range: query_range.clone(),
                        };
                    }
                }
                LocalPhmmDriverState::Finalize { ref_range, query_range } => {
                    return visitor.finalize(self.phmm.as_view(), ref_range.clone(), query_range.clone());
                }
            }
        }
    }
}
