//! Implementations for examining the parameters encountered when following a
//! given alignment through a pHMM.
//!
//! See [`GlobalPhmm::visit_params`], [`LocalPhmm::visit_params`],
//! [`SemiLocalPhmm::visit_params`], and [`DomainPhmm::visit_params`].
//!
//! <div class="warning note">
//!
//! **Note**
//!
//! You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
//! use these methods.
//!
//! </div>

use crate::{
    alignment::{
        StatesSequence,
        phmm::{
            CorePhmm, DomainPhmm, GetLayer, GetModule, GlobalPhmm, LocalPhmm, PhmmError, PhmmNumber, PhmmState,
            SemiLocalPhmm,
            indexing::{Begin, DpIndex, End, FirstMatch, LastMatch, PhmmIndex, PhmmIndexable, SeqIndex},
            modules::{DomainModule, SemiLocalModule},
        },
    },
    data::{ByteIndexMap, cigar::Ciglet},
};
use std::ops::Range;

impl<T: PhmmNumber, const S: usize> DomainModule<T, S> {
    /// Lazily compute the score for skipping `inserted` residues at the
    /// beginning of the query. This should only be used for diagnostics or
    /// testing, otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    fn get_begin_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        if inserted.is_empty() {
            self.start_to_end
        } else {
            // Special casing needed in case insert_to_insert is infinite,
            // causing a NAN to appear when multiplied by 0
            let insert_to_insert = if inserted.len() > 1 {
                T::cast_from(inserted.len() - 1) * self.insert_to_insert
            } else {
                T::ZERO
            };

            self.start_to_insert
                + self.insert_to_end
                + (inserted
                    .iter()
                    .map(|x| self.background_emission[mapping.to_index(*x)])
                    .fold(T::ZERO, |acc, elem| acc + elem)
                    + insert_to_insert)
        }
    }

    /// Lazily compute the score for skipping `inserted` residues at the end of
    /// the query. This should only be used for diagnostics or testing,
    /// otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    fn get_end_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        if inserted.is_empty() {
            self.start_to_end
        } else {
            // Special casing needed in case insert_to_insert is infinite,
            // causing a NAN to appear when multiplied by 0
            let insert_to_insert = if inserted.len() > 1 {
                T::cast_from(inserted.len() - 1) * self.insert_to_insert
            } else {
                T::ZERO
            };

            self.start_to_insert
                + self.insert_to_end
                + (inserted
                    .iter()
                    .rev()
                    .map(|x| self.background_emission[mapping.to_index(*x)])
                    .fold(T::ZERO, |acc, elem| acc + elem)
                    + insert_to_insert)
        }
    }
}

impl<T: PhmmNumber, const S: usize> DomainPhmm<T, S> {
    /// Lazily compute the score for skipping `inserted` residues at the
    /// beginning of the query. This should only be used for diagnostics or
    /// testing, otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    #[inline]
    fn get_begin_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        self.begin().get_begin_score(inserted, mapping)
    }

    /// Lazily compute the score for skipping `inserted` residues at the end of
    /// the query. This should only be used for diagnostics or testing,
    /// otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    #[inline]
    fn get_end_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        self.end().get_end_score(inserted, mapping)
    }
}

impl<T: PhmmNumber, const S: usize> LocalPhmm<T, S> {
    /// Lazily compute the score for skipping `inserted` residues at the
    /// beginning of the query. This should only be used for diagnostics or
    /// testing, otherwise [`PrecomputedLocalModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedLocalModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedLocalModule
    fn get_begin_internal_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        self.begin().internal_params.get_begin_score(inserted, mapping)
    }

    /// Lazily compute the score for skipping to `index` while entering the
    /// [`CorePhmm`]. This should only be used for diagnostics or testing,
    /// otherwise [`PrecomputedLocalModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedLocalModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedLocalModule
    fn get_begin_external_score(&self, index: impl PhmmIndex) -> T {
        self.begin().external_params.get_score(index)
    }

    /// Lazily compute the score for exiting the [`CorePhmm`] from `index`. This
    /// should only be used for diagnostics or testing, otherwise
    /// [`PrecomputedLocalModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedLocalModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedLocalModule
    fn get_end_external_score(&self, index: impl PhmmIndex) -> T {
        self.end().external_params.get_score(index)
    }

    /// Lazily compute the score for skipping `inserted` residues at the end of
    /// the query. This should only be used for diagnostics or testing,
    /// otherwise [`PrecomputedLocalModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedLocalModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedLocalModule
    fn get_end_internal_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        self.end().internal_params.get_end_score(inserted, mapping)
    }
}

/// A wrapper around a pHMM parameter of type `T` which also contains metadata
/// about that parameter.
pub struct PhmmParam<T> {
    /// The parameter value, in negative log space.
    pub param: T,
    /// Information about the parameter and where it occurred in the model.
    pub kind:  PhmmParamKind<T>,
}

/// The location of a module in a pHMM.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum ModuleLocation {
    /// The module is placed at the beginning of the pHMM.
    Begin,
    /// The module is placed at the end of the pHMM.
    End,
}

/// Metadata about a pHMM parameter, used in [`PhmmParam`].
pub enum PhmmParamKind<T> {
    /// An emission parameter from a match state.
    EmissionMatch {
        /// The layer index of the match state.
        layer_idx:   DpIndex,
        /// The residue index that was emitted.
        residue_idx: usize,
    },
    /// An emission parameter from an insert state.
    EmissionInsert {
        /// The layer index of the insert state.
        layer_idx:   DpIndex,
        /// The residue index that was emitted.
        residue_idx: usize,
    },
    /// A transition parameter between two states in the [`CorePhmm`].
    Transition {
        /// The index of the layer that is being transitioned out of.
        starting_layer: DpIndex,
        /// The state that is being transitioned out of.
        starting_state: PhmmState,
        /// The state that is being transitioned into.
        ending_state:   PhmmState,
    },
    /// The (composite) parameter emitted by a [`LocalModule`] at the beginning
    /// or end of a [`LocalPhmm`].
    ///
    /// [`LocalModule`]: super::modules::LocalModule
    LocalModule {
        /// The (composite) parameter within the module, associated with
        /// skipping residues at the beginning of the query.
        internal_param: T,
        /// The skipped residues at the beginning of the query.
        skipped:        Vec<u8>,
        /// The parameter connecting the module to the [`CorePhmm`], associated
        /// with skipping to a layer in the pHMM.
        external_param: T,
        /// The layer in the [`CorePhmm`] which is entered.
        to_layer:       DpIndex,
        /// The location of the module.
        loc:            ModuleLocation,
    },
    /// The (composite) parameter emitted by a [`DomainModule`] at the beginning
    /// or end of a [`DomainPhmm`].
    ///
    /// [`DomainModule`]: super::modules::DomainModule
    DomainModule {
        /// The skipped residues at the beginning of the query.
        skipped: Vec<u8>,
        /// The location of the module.
        loc:     ModuleLocation,
    },
    /// The parameter emitted by a [`SemiLocalModule`] at the beginning or end
    /// of a [`SemiLocalPhmm`].
    ///
    /// [`SemiLocalModule`]: super::modules::SemiLocalModule
    SemiLocalModule {
        /// The layer in the [`CorePhmm`] which is entered/exited.
        core_layer: DpIndex,
        /// The location of the module.
        loc:        ModuleLocation,
    },
}

impl<T: PhmmNumber> PhmmParam<T> {
    /// Constructs a [`PhmmParam`] with kind [`EmissionMatch`].
    ///
    /// [`EmissionMatch`]: PhmmParamKind::EmissionMatch
    fn new_emission_match(param: T, layer_idx: impl PhmmIndex, residue_idx: usize, phmm: &impl PhmmIndexable) -> Self {
        Self {
            param,
            kind: PhmmParamKind::EmissionMatch {
                layer_idx: phmm.to_dp_index(layer_idx),
                residue_idx,
            },
        }
    }

    /// Constructs a [`PhmmParam`] with kind [`EmissionInsert`].
    ///
    /// [`EmissionInsert`]: PhmmParamKind::EmissionInsert
    fn new_emission_insert(param: T, layer_idx: impl PhmmIndex, residue_idx: usize, phmm: &impl PhmmIndexable) -> Self {
        Self {
            param,
            kind: PhmmParamKind::EmissionInsert {
                layer_idx: phmm.to_dp_index(layer_idx),
                residue_idx,
            },
        }
    }

    /// Constructs a [`PhmmParam`] with kind [`Transition`].
    ///
    /// [`Transition`]: PhmmParamKind::Transition
    fn new_transition(
        param: T, starting_layer: impl PhmmIndex, starting_state: PhmmState, ending_state: PhmmState,
        phmm: &impl PhmmIndexable,
    ) -> Self {
        Self {
            param,
            kind: PhmmParamKind::Transition {
                starting_layer: phmm.to_dp_index(starting_layer),
                starting_state,
                ending_state,
            },
        }
    }

    /// Constructs a [`PhmmParam`] with kind [`LocalModule`].
    ///
    /// [`LocalModule`]: PhmmParamKind::LocalModule
    fn new_local_module(
        loc: ModuleLocation, internal_param: T, skipped: Vec<u8>, external_param: T, to_layer: impl PhmmIndex,
        phmm: &impl PhmmIndexable,
    ) -> Self {
        Self {
            param: internal_param + external_param,
            kind:  PhmmParamKind::LocalModule {
                internal_param,
                skipped,
                external_param,
                to_layer: phmm.to_dp_index(to_layer),
                loc,
            },
        }
    }

    /// Constructs a [`PhmmParam`] with kind [`DomainModule`].
    ///
    /// [`DomainModule`]: PhmmParamKind::DomainModule
    fn new_domain_module(loc: ModuleLocation, param: T, skipped: Vec<u8>) -> Self {
        Self {
            param,
            kind: PhmmParamKind::DomainModule { skipped, loc },
        }
    }

    /// Constructs a [`PhmmParam`] with kind [`SemiLocalModule`].
    ///
    /// [`SemiLocalModule`]: PhmmParamKind::SemiLocalModule
    fn new_semilocal_module(loc: ModuleLocation, param: T, to_layer: impl PhmmIndex, phmm: &impl PhmmIndexable) -> Self {
        Self {
            param,
            kind: PhmmParamKind::SemiLocalModule {
                core_layer: phmm.to_dp_index(to_layer),
                loc,
            },
        }
    }
}

/// Given a parameter encountered in the pHMM, calls the closure `f` on it and
/// updates the score.
#[inline]
fn call_f<T, F>(f: &mut F, score: &mut T, param: PhmmParam<T>)
where
    T: PhmmNumber,
    F: FnMut(PhmmParam<T>), {
    *score += param.param;
    f(param);
}

/// Given a pHMM alignment specified by `ciglets`, call a closure on the
/// parameters encountered as the alignment is traversed through the portion
/// corresponding to the [`CorePhmm`].
///
/// The return value is the state from which the [`CorePhmm`] is exited (which
/// is the last one before the END state if there is no early exit).
///
/// `ref_range` specifies the range of the pHMM's underlying "reference" that
/// should be consumed while scoring. For global alignment or domain alignment,
/// this is the full length of the pHMM (as can be found with [`seq_len`]). For
/// other forms of alignment, this may be a smaller range.
///
/// `seq_in_alignment` should contain exactly the portion of the sequence
/// matched by the model in the given alignment.
///
/// `score` is the starting score to tally into (provided since order of
/// operations may matter with floating point arithmetic).
///
/// `f` is the closure which is to be called on each parameter, of type
/// [`PhmmParam`] (which is a helper struct containing any relevant information
/// about the parameter). The closure may mutate state as needed.
///
/// This function does not include the following parameters:
///
/// - Transitions from the BEGIN state
/// - Transitions into the END state
/// - Transitions from a starting or ending module into a match state
///
/// [`seq_len`]: PhmmIndexable::seq_len
fn visit_params_core<T, const S: usize, F>(
    core: &CorePhmm<T, S>, mapping: &ByteIndexMap<S>, seq_in_alignment: &[u8], ref_range: Range<usize>, ciglets: &[Ciglet],
    score: &mut T, f: &mut F,
) -> Result<PhmmState, PhmmError>
where
    T: PhmmNumber,
    F: FnMut(PhmmParam<T>), {
    use PhmmState::*;

    // TODO: Improve ergonomics later...
    let mut op_iter = ciglets.iter().flat_map(|ciglet| std::iter::repeat_n(ciglet.op, ciglet.inc));

    let Some(op) = op_iter.next() else {
        if ref_range.is_empty() {
            // If the CIGAR string is empty, then either the path entered the
            // BEGIN state and then exited, or entered the END state and then
            // exited. In either case, a match state was passed through.
            return Ok(Match);
        }
        return Err(PhmmError::FullModelNotUsed);
    };

    // Index into the query
    let mut i = 0;

    // Initialize current score, starting state, and starting layer index
    let (mut state, mut j) = match PhmmState::from_op(op)? {
        Delete => {
            // Cannot enter into a delete state except from BEGIN.
            if ref_range.start != 0 {
                return Err(PhmmError::InvalidPath);
            }
            // After entering the delete state, we are now in layer 1
            (Delete, 1)
        }
        Match => {
            let x_idx = mapping.to_index(seq_in_alignment[i]);
            // Need to look at the previous layer to get transitions into this
            // layer
            let layer_idx = SeqIndex(ref_range.start).prev_index(core);
            call_f(
                f,
                score,
                PhmmParam::new_emission_match(core.get_layer(layer_idx).emission_match[x_idx], layer_idx, x_idx, core),
            );
            i += 1;
            (Match, ref_range.start + 1)
        }
        Insert => {
            // Cannot enter into a insert state except from BEGIN.
            if ref_range.start != 0 {
                return Err(PhmmError::InvalidPath);
            }
            // Insertion before reference position n has emission parameters and
            // transition parameters stored in layer n
            let x_idx = mapping.to_index(seq_in_alignment[i]);
            call_f(
                f,
                score,
                PhmmParam::new_emission_insert(core.get_layer(Begin).emission_insert[x_idx], Begin, x_idx, core),
            );
            i += 1;
            // After entering the delete state, we are still in begin layer
            (Insert, 0)
        }
    };

    // Main loop: first score transition out of current layer index j, then
    // score emission in new state if applicable
    for op in op_iter {
        let layer = core.get_layer(DpIndex(j));

        match PhmmState::from_op(op)? {
            Match => {
                let x_idx = mapping.to_index(seq_in_alignment[i]);
                // Must perform the two additions separately, since floating
                // point addition is not associative
                call_f(
                    f,
                    score,
                    PhmmParam::new_transition(layer.transition[(state, Match)], DpIndex(j), state, Match, core),
                );
                call_f(
                    f,
                    score,
                    PhmmParam::new_emission_match(layer.emission_match[x_idx], DpIndex(j), x_idx, core),
                );
                state = Match;
                i += 1;
                j += 1;
            }
            Insert => {
                let x_idx = mapping.to_index(seq_in_alignment[i]);
                // Must perform the two additions separately, since floating
                // point addition is not associative
                call_f(
                    f,
                    score,
                    PhmmParam::new_transition(layer.transition[(state, Insert)], DpIndex(j), state, Insert, core),
                );
                call_f(
                    f,
                    score,
                    PhmmParam::new_emission_insert(layer.emission_insert[x_idx], DpIndex(j), x_idx, core),
                );
                state = Insert;
                i += 1;
            }
            Delete => {
                call_f(
                    f,
                    score,
                    PhmmParam::new_transition(layer.transition[(state, Delete)], DpIndex(j), state, Delete, core),
                );
                state = Delete;
                j += 1;
            }
        }
    }

    // i and j were incremented 1 past their last value
    if i != seq_in_alignment.len() {
        return Err(PhmmError::FullSeqNotUsed);
    }
    if j != ref_range.end {
        return Err(PhmmError::FullModelNotUsed);
    }

    Ok(state)
}

/// Finds information about the optimal path at the beginning of a [`LocalPhmm`]
/// or [`SemiLocalPhmm`], assuming that the `ref_range` starts at 0.
///
/// There may be ambiguity as to the path taken at the beginning of such a pHMM.
/// It may pass through the Begin state, or it may skip directly to the first
/// match state. This function determines which has the minimal score.
///
/// The function accepts the current score so far (which may be 0 or the
/// internal parameter for [`LocalPhmm`]). Due to floating point error, this is
/// required to determine the optimal path.
///
/// The return values are the parameter for transitioning from the module to the
/// pHMM, the index of the layer that is entered, the parameter for
/// transitioning from the Begin state to the next state (if the Begin state is
/// passed through), and the state that is transitioned to after the module and
/// Begin state.
fn resolve_ambiguous_start<T: PhmmNumber, const S: usize>(
    begin_module: &SemiLocalModule<T>, core: &CorePhmm<T, S>, score: T, ciglets: &[Ciglet],
) -> Result<(T, DpIndex, Option<T>, PhmmState), PhmmError> {
    use PhmmState::*;

    let first_op = ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?;
    let first_state = PhmmState::from_op(first_op)?;
    let external_begin_param_through_begin = begin_module.get_score(Begin);
    let transition_from_begin = core.get_layer(Begin).transition[(Match, first_state)];
    let score_through_begin = score + external_begin_param_through_begin + transition_from_begin;

    let (external_begin_param, to_layer, transition_from_begin) = match first_state {
        Delete | Insert => (
            external_begin_param_through_begin,
            core.to_dp_index(Begin),
            Some(transition_from_begin),
        ),
        Match => {
            // Handle ambiguity: it is unclear whether we pass through the BEGIN
            // state of the core pHMM or not
            let external_begin_param_skip_begin = begin_module.get_score(FirstMatch);
            let score_skipping_begin = score + external_begin_param_skip_begin;

            if score_through_begin <= score_skipping_begin {
                (
                    external_begin_param_through_begin,
                    core.to_dp_index(Begin),
                    Some(transition_from_begin),
                )
            } else {
                (external_begin_param_skip_begin, core.to_dp_index(FirstMatch), None)
            }
        }
    };

    Ok((external_begin_param, to_layer, transition_from_begin, first_state))
}

/// Finds information about the optimal path at the end of a [`LocalPhmm`] or
/// [`SemiLocalPhmm`], assuming that the `ref_range` ends at the length of the
/// [`CorePhmm`].
///
/// There may be ambiguity as to the path taken at the end of such a pHMM. It
/// may pass through the End state, or it may exit directly from the last match
/// state. This function determines which has the minimal score.
///
/// The function accepts the current score so far, and the `internal_end_param`
/// for [`LocalPhmm`]. Due to floating point error, this is required to
/// determine the optimal path.
///
/// <div class="warning important">
///
/// **Important**
///
/// `internal_end_param` must be passed for a [`LocalPhmm`], and should be
/// `None` for [`SemiLocalPhmm`].
///
/// </div>
///
/// The return values are the parameter for transitioning from pHMM to the
/// module, the index of the layer that is exited, and the parameter for
/// transitioning into the End state (if the End state is passed through).
fn resolve_ambiguous_end<T: PhmmNumber, const S: usize>(
    end_module: &SemiLocalModule<T>, core: &CorePhmm<T, S>, score: T, final_state: PhmmState, internal_end_param: Option<T>,
) -> (T, DpIndex, Option<T>) {
    use PhmmState::*;

    let transition_to_end = core.get_layer(LastMatch).transition[(final_state, Match)];
    let external_end_param_through_end = end_module.get_score(End);
    let end_param_through_end = if let Some(internal_end_param) = internal_end_param {
        internal_end_param + external_end_param_through_end
    } else {
        external_end_param_through_end
    };
    let score_through_end = score + transition_to_end + end_param_through_end;

    let (transition_to_end, external_end_param, from_layer) = match final_state {
        Delete | Insert => (Some(transition_to_end), external_end_param_through_end, core.to_dp_index(End)),
        Match => {
            // Handle ambiguity: it is unclear whether we pass through
            // the END state of the core pHMM or not
            let external_end_param_skip_end = end_module.get_score(LastMatch);
            let end_param_skip_end = if let Some(internal_end_param) = internal_end_param {
                internal_end_param + external_end_param_skip_end
            } else {
                external_end_param_skip_end
            };
            let score_skipping_end = score + end_param_skip_end;

            if score_through_end <= score_skipping_end {
                (Some(transition_to_end), external_end_param_through_end, core.to_dp_index(End))
            } else {
                (None, external_end_param_skip_end, core.to_dp_index(LastMatch))
            }
        }
    };

    (external_end_param, from_layer, transition_to_end)
}

impl<T: PhmmNumber, const S: usize> GlobalPhmm<T, S> {
    /// Performs an action on each parameter that is encountered as a particular
    /// alignment is taken through a pHMM.
    ///
    /// The action to perform is specified by closure `f`, which takes a
    /// [`PhmmParam`] (which is a helper struct containing any relevant
    /// information about the parameter). The closure may mutate state as
    /// needed.
    ///
    /// The score is also summed and returned. This is designed to give the
    /// exact same score as [`viterbi`] when the best `alignment` is passed,
    /// performing all arithmetic operations in the same order so as not to
    /// change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are `M`, `=`, `X`, `I`, and `D`.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    ///
    /// [`viterbi`]: GlobalPhmm::viterbi
    pub fn visit_params<Q, A, F>(&self, seq: Q, alignment: A, f: F) -> Result<T, PhmmError>
    where
        Q: AsRef<[u8]>,
        A: AsRef<[Ciglet]>,
        F: FnMut(PhmmParam<T>), {
        self.visit_params_helper(seq.as_ref(), alignment.as_ref(), f)
    }

    /// See [`visit_params`]. This is a helper function to reduce
    /// monomorphization.
    ///
    /// [`visit_params`]: GlobalPhmm::visit_params
    fn visit_params_helper<F>(&self, seq: &[u8], ciglets: &[Ciglet], mut f: F) -> Result<T, PhmmError>
    where
        F: FnMut(PhmmParam<T>), {
        use PhmmState::*;

        let mut score = T::ZERO;

        let first_op = ciglets.iter().next().ok_or(PhmmError::FullModelNotUsed)?.op;
        let first_state = PhmmState::from_op(first_op)?;

        // Get transition from BEGIN state the actual first state

        call_f(
            &mut f,
            &mut score,
            PhmmParam::new_transition(
                self.get_layer(Begin).transition[(Match, first_state)],
                Begin,
                Match,
                first_state,
                self,
            ),
        );

        // Add score from first state after START until final state before END
        let final_state = visit_params_core(
            self.core(),
            self.mapping(),
            seq,
            0..self.seq_len(),
            ciglets,
            &mut score,
            &mut f,
        )?;

        // Add transition into END state
        call_f(
            &mut f,
            &mut score,
            PhmmParam::new_transition(
                self.get_layer(LastMatch).transition[(final_state, Match)],
                LastMatch,
                final_state,
                Match,
                self,
            ),
        );

        Ok(score)
    }
}

impl<T: PhmmNumber, const S: usize> LocalPhmm<T, S> {
    /// Performs an action on each parameter that is encountered as a particular
    /// alignment is taken through a pHMM.
    ///
    /// The action to perform is specified by closure `f`, which takes a
    /// [`PhmmParam`] (which is a helper struct containing any relevant
    /// information about the parameter). The closure may mutate state as
    /// needed.
    ///
    /// The score is also summed and returned. This is designed to give the
    /// exact same score as [`viterbi`] when the best `alignment` is passed,
    /// performing all arithmetic operations in the same order so as not to
    /// change the floating point error.
    ///
    /// There may be ambiguity as to the path taken at the beginning and end.
    /// This is resolved by choosing the path with the minimal score.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are `M`, `=`, `X`, `I`, and `D`. If `path`
    /// does not correspond to a valid path through a [`LocalPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    ///
    /// [`viterbi`]: LocalPhmm::viterbi
    #[allow(clippy::too_many_lines)]
    pub fn visit_params<Q, A, F>(&self, seq: Q, alignment: A, ref_range: Range<usize>, f: F) -> Result<T, PhmmError>
    where
        Q: AsRef<[u8]>,
        A: AsRef<[Ciglet]>,
        F: FnMut(PhmmParam<T>), {
        self.visit_params_helper(seq.as_ref(), alignment.as_ref(), ref_range, f)
    }

    /// See [`visit_params`]. This is a helper function to reduce
    /// monomorphization.
    ///
    /// [`visit_params`]: LocalPhmm::visit_params
    #[allow(clippy::too_many_lines)]
    fn visit_params_helper<F>(
        &self, seq: &[u8], mut ciglets: &[Ciglet], ref_range: Range<usize>, mut f: F,
    ) -> Result<T, PhmmError>
    where
        F: FnMut(PhmmParam<T>), {
        use PhmmState::*;

        let mut score = T::ZERO;

        // Split query between begin module, core pHMM, and end module
        let (begin_seq, seq, end_seq) = {
            let begin_inserted = ciglets.next_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
            let end_inserted = ciglets.next_back_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
            if ciglets.is_empty() {
                return Ok(self.visit_params_empty_alignment(seq, f));
            }

            let (seq, end_seq) = seq.split_at(seq.len() - end_inserted);
            let (begin_seq, seq) = seq.split_at(begin_inserted);
            (begin_seq, seq, end_seq)
        };

        // Add contribution from internal parameters of begin module
        let internal_begin_param = self.get_begin_internal_score(begin_seq, self.mapping());

        // Add the contribution of the external parameters into the core pHMM
        // and the transitions out of the BEGIN state
        if ref_range.start == 0 {
            let (external_begin_param, to_layer, transition_from_begin, first_state) =
                resolve_ambiguous_start(&self.begin().external_params, self.core(), score, ciglets)?;

            call_f(
                &mut f,
                &mut score,
                PhmmParam::new_local_module(
                    ModuleLocation::Begin,
                    internal_begin_param,
                    begin_seq.to_vec(),
                    external_begin_param,
                    to_layer,
                    self,
                ),
            );

            if let Some(transition_from_begin) = transition_from_begin {
                call_f(
                    &mut f,
                    &mut score,
                    PhmmParam::new_transition(transition_from_begin, Begin, Match, first_state, self),
                );
            }
        } else {
            // We must enter into a match state if not going through BEGIN
            if PhmmState::from_op(ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?)? != Match {
                return Err(PhmmError::InvalidPath);
            }

            call_f(
                &mut f,
                &mut score,
                PhmmParam::new_local_module(
                    ModuleLocation::Begin,
                    internal_begin_param,
                    begin_seq.to_vec(),
                    self.get_begin_external_score(SeqIndex(ref_range.start)),
                    SeqIndex(ref_range.start),
                    self,
                ),
            );
        }

        // Add the contribution from the core pHMM excluding transitions into
        // the END state, and obtain the final state before the END state
        let final_state = visit_params_core(
            self.core(),
            self.mapping(),
            seq,
            ref_range.clone(),
            ciglets,
            &mut score,
            &mut f,
        )?;

        // Compute the contribution from the internal parameters of the end
        // module, which involves processing any soft clipping at the end of the
        // alignment
        let internal_end_param = self.get_end_internal_score(end_seq, self.mapping());

        if ref_range.end == self.seq_len() {
            let (external_end_param, from_layer, transition_to_end) = resolve_ambiguous_end(
                &self.end().external_params,
                self.core(),
                score,
                final_state,
                Some(internal_end_param),
            );

            if let Some(transition_to_end) = transition_to_end {
                call_f(
                    &mut f,
                    &mut score,
                    PhmmParam::new_transition(transition_to_end, LastMatch, final_state, Match, self),
                );
            }

            call_f(
                &mut f,
                &mut score,
                PhmmParam::new_local_module(
                    ModuleLocation::End,
                    internal_end_param,
                    end_seq.to_vec(),
                    external_end_param,
                    from_layer,
                    self,
                ),
            );
        } else {
            // We must exit from a match state if not going through END
            if final_state != Match {
                return Err(PhmmError::InvalidPath);
            }

            // Subtract 1 since range is end-exclusive
            call_f(
                &mut f,
                &mut score,
                PhmmParam::new_local_module(
                    ModuleLocation::End,
                    internal_end_param,
                    end_seq.to_vec(),
                    self.get_end_external_score(SeqIndex(ref_range.end - 1)),
                    SeqIndex(ref_range.end - 1),
                    self,
                ),
            );
        }

        Ok(score)
    }

    /// Performs an action on each parameter that is encountered as an empty
    /// alignment (all soft clipping) is taken through a pHMM.
    ///
    /// The action to perform is specified by closure `f`, which takes a
    /// [`PhmmParam`] (which is a helper struct containing any relevant
    /// information about the parameter). The closure may mutate state as
    /// needed.
    ///
    /// The score is also summed and returned. This is designed to give the
    /// exact same score as [`viterbi`] when the best alignment is indeed empty,
    /// performing all arithmetic operations in the same order so as not to
    /// change the floating point error.
    ///
    /// There is ambiguity as to which query bases are consumed by the module at
    /// the beginning of the pHMM and the model at the end. This is resolved by
    /// choosing the grouping with the minimal score. If there is a tie,
    /// residues are included in the module at the beginning.
    ///
    /// [`viterbi`]: LocalPhmm::viterbi
    fn visit_params_empty_alignment<F>(&self, seq: &[u8], mut f: F) -> T
    where
        F: FnMut(PhmmParam<T>), {
        let mut best_i = 0;
        let mut best_state = self.to_dp_index(Begin);
        let mut best_score = T::INFINITY;

        for i in 0..=seq.len() {
            let (inserted_begin, inserted_end) = seq.split_at(i);
            for through_state in [self.to_dp_index(Begin), self.to_dp_index(End)] {
                let begin_score = self.begin().internal_params.get_begin_score(inserted_begin, self.mapping())
                    + self.get_begin_external_score(through_state);
                let end_score =
                    self.get_end_internal_score(inserted_end, self.mapping()) + self.get_end_external_score(through_state);
                let score = begin_score + end_score;

                if score < best_score {
                    best_i = i;
                    best_state = through_state;
                    best_score = score;
                }
            }
        }

        let (inserted_begin, inserted_end) = seq.split_at(best_i);
        f(PhmmParam::new_local_module(
            ModuleLocation::Begin,
            self.begin().internal_params.get_begin_score(inserted_begin, self.mapping()),
            inserted_begin.to_vec(),
            self.get_begin_external_score(best_state),
            best_state,
            self,
        ));
        f(PhmmParam::new_local_module(
            ModuleLocation::End,
            self.get_end_internal_score(inserted_end, self.mapping()),
            inserted_end.to_vec(),
            self.get_end_external_score(best_state),
            best_state,
            self,
        ));
        best_score
    }
}

impl<T: PhmmNumber, const S: usize> DomainPhmm<T, S> {
    /// Performs an action on each parameter that is encountered as a particular
    /// alignment is taken through a pHMM.
    ///
    /// The action to perform is specified by closure `f`, which takes a
    /// [`PhmmParam`] (which is a helper struct containing any relevant
    /// information about the parameter). The closure may mutate state as
    /// needed.
    ///
    /// The score is also summed and returned. This is designed to give the
    /// exact same score as [`viterbi`] when the best `alignment` is passed,
    /// performing all arithmetic operations in the same order so as not to
    /// change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are `M`, `=`, `X`, `I`, and `D`. If `path`
    /// does not correspond to a valid path through a [`DomainPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    ///
    /// [`viterbi`]: DomainPhmm::viterbi
    pub fn visit_params<Q, A, F>(&self, seq: Q, alignment: A, f: F) -> Result<T, PhmmError>
    where
        Q: AsRef<[u8]>,
        A: AsRef<[Ciglet]>,
        F: FnMut(PhmmParam<T>), {
        self.visit_params_helper(seq.as_ref(), alignment.as_ref(), f)
    }

    /// See [`visit_params`]. This is a helper function to reduce
    /// monomorphization.
    ///
    /// [`visit_params`]: DomainPhmm::visit_params
    fn visit_params_helper<F>(&self, seq: &[u8], mut ciglets: &[Ciglet], mut f: F) -> Result<T, PhmmError>
    where
        F: FnMut(PhmmParam<T>), {
        use PhmmState::*;

        let mut score = T::ZERO;

        // Split query between begin module, core pHMM, and end module
        let (begin_seq, seq, end_seq) = {
            let begin_inserted = ciglets.next_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
            let end_inserted = ciglets.next_back_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);

            let (seq, end_seq) = seq.split_at(seq.len() - end_inserted);
            let (begin_seq, seq) = seq.split_at(begin_inserted);
            (begin_seq, seq, end_seq)
        };

        // Add contribution from parameters of begin module
        call_f(
            &mut f,
            &mut score,
            PhmmParam::new_domain_module(
                ModuleLocation::Begin,
                self.get_begin_score(begin_seq, self.mapping()),
                begin_seq.to_vec(),
            ),
        );

        // Add the transitions out of the BEGIN state
        let first_op = ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?;
        let first_state = PhmmState::from_op(first_op)?;
        call_f(
            &mut f,
            &mut score,
            PhmmParam::new_transition(
                self.get_layer(Begin).transition[(Match, first_state)],
                Begin,
                Match,
                first_state,
                self,
            ),
        );

        // Add the contribution from the core pHMM excluding transitions into
        // the END state, and obtain the final state before the END state
        let final_state = visit_params_core(
            self.core(),
            self.mapping(),
            seq,
            0..self.seq_len(),
            ciglets,
            &mut score,
            &mut f,
        )?;

        // Add the transition into the END state
        call_f(
            &mut f,
            &mut score,
            PhmmParam::new_transition(
                self.get_layer(LastMatch).transition[(final_state, Match)],
                LastMatch,
                final_state,
                Match,
                self,
            ),
        );

        // Add the contribution from parameters of the end module
        call_f(
            &mut f,
            &mut score,
            PhmmParam::new_domain_module(
                ModuleLocation::End,
                self.get_end_score(end_seq, self.mapping()),
                end_seq.to_vec(),
            ),
        );

        Ok(score)
    }
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<T, S> {
    /// Performs an action on each parameter that is encountered as a particular
    /// alignment is taken through a pHMM.
    ///
    /// The action to perform is specified by closure `f`, which takes a
    /// [`PhmmParam`] (which is a helper struct containing any relevant
    /// information about the parameter). The closure may mutate state as
    /// needed.
    ///
    /// The score is also summed and returned. This is designed to give the
    /// exact same score as [`viterbi`] when the best `alignment` is passed,
    /// performing all arithmetic operations in the same order so as not to
    /// change the floating point error.
    ///
    /// There may be ambiguity as to the path taken at the beginning and end.
    /// This is resolved by choosing the path with the minimal score.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are `M`, `=`, `X`, `I`, and `D`. If `path`
    /// does not correspond to a valid path through a [`SemiLocalPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    ///
    /// [`viterbi`]: SemiLocalPhmm::viterbi
    pub fn visit_params<Q, A, F>(&self, seq: Q, alignment: A, ref_range: Range<usize>, f: F) -> Result<T, PhmmError>
    where
        Q: AsRef<[u8]>,
        A: AsRef<[Ciglet]>,
        F: FnMut(PhmmParam<T>), {
        self.visit_params_helper(seq.as_ref(), alignment.as_ref(), ref_range, f)
    }

    /// See [`visit_params`]. This is a helper function to reduce
    /// monomorphization.
    ///
    /// [`visit_params`]: SemiLocalPhmm::visit_params
    fn visit_params_helper<F>(
        &self, seq: &[u8], ciglets: &[Ciglet], ref_range: Range<usize>, mut f: F,
    ) -> Result<T, PhmmError>
    where
        F: FnMut(PhmmParam<T>), {
        use PhmmState::*;

        let mut score = T::ZERO;

        if ciglets.peek_op().is_none() {
            return Ok(self.visit_params_empty_alignment(f));
        }

        // Add the contribution of the external parameters into the core pHMM
        // and the transitions out of the BEGIN state
        if ref_range.start == 0 {
            let (begin_param, to_layer, transition_from_begin, first_state) =
                resolve_ambiguous_start(self.begin(), self.core(), score, ciglets)?;

            call_f(
                &mut f,
                &mut score,
                PhmmParam::new_semilocal_module(ModuleLocation::Begin, begin_param, to_layer, self),
            );

            if let Some(transition_from_begin) = transition_from_begin {
                call_f(
                    &mut f,
                    &mut score,
                    PhmmParam::new_transition(transition_from_begin, Begin, Match, first_state, self),
                );
            }
        } else {
            // We must enter into a match state if not going through BEGIN
            if PhmmState::from_op(ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?)? != Match {
                return Err(PhmmError::InvalidPath);
            }

            call_f(
                &mut f,
                &mut score,
                PhmmParam::new_semilocal_module(
                    ModuleLocation::Begin,
                    self.get_begin_score(SeqIndex(ref_range.start)),
                    SeqIndex(ref_range.start),
                    self,
                ),
            );
        }

        // Add the contribution from the core pHMM excluding transitions into
        // the END state, and obtain the final state before the END state
        let final_state = visit_params_core(
            self.core(),
            self.mapping(),
            seq,
            ref_range.clone(),
            ciglets,
            &mut score,
            &mut f,
        )?;

        if ref_range.end == self.seq_len() {
            let (end_param, from_layer, transition_to_end) =
                resolve_ambiguous_end(self.end(), self.core(), score, final_state, None);

            if let Some(transition_to_end) = transition_to_end {
                call_f(
                    &mut f,
                    &mut score,
                    PhmmParam::new_transition(transition_to_end, LastMatch, final_state, Match, self),
                );
            }

            call_f(
                &mut f,
                &mut score,
                PhmmParam::new_semilocal_module(ModuleLocation::End, end_param, from_layer, self),
            );
        } else {
            // We must exit from a match state if not going through END
            if final_state != Match {
                return Err(PhmmError::InvalidPath);
            }

            // Subtract 1 since range is end-exclusive
            call_f(
                &mut f,
                &mut score,
                PhmmParam::new_semilocal_module(
                    ModuleLocation::End,
                    self.get_end_score(SeqIndex(ref_range.end - 1)),
                    SeqIndex(ref_range.end - 1),
                    self,
                ),
            );
        }

        Ok(score)
    }

    /// Performs an action on each parameter that is encountered as an empty
    /// alignment is taken through a pHMM. This means that the original sequence
    /// was empty too!
    ///
    /// The action to perform is specified by closure `f`, which takes a
    /// [`PhmmParam`] (which is a helper struct containing any relevant
    /// information about the parameter). The closure may mutate state as
    /// needed.
    ///
    /// The score is also summed and returned. This is designed to give the
    /// exact same score as [`viterbi`] when the best alignment is indeed empty,
    /// performing all arithmetic operations in the same order so as not to
    /// change the floating point error.
    ///
    /// There is ambiguity as to whether to pass through the Begin or End state.
    /// This is resolved by choosing the state with the minimal score. If there
    /// is a tie, the Begin state is used.
    ///
    /// [`viterbi`]: SemiLocalPhmm::viterbi
    fn visit_params_empty_alignment<F>(&self, mut f: F) -> T
    where
        F: FnMut(PhmmParam<T>), {
        let score_through_begin = self.get_begin_score(Begin) + self.get_end_score(Begin);
        let score_through_end = self.get_begin_score(End) + self.get_end_score(End);

        let (score, through_state) = if score_through_begin <= score_through_end {
            (score_through_begin, self.to_dp_index(Begin))
        } else {
            (score_through_end, self.to_dp_index(End))
        };

        f(PhmmParam::new_semilocal_module(
            ModuleLocation::Begin,
            self.get_begin_score(through_state),
            through_state,
            self,
        ));
        f(PhmmParam::new_semilocal_module(
            ModuleLocation::End,
            self.get_end_score(through_state),
            through_state,
            self,
        ));
        score
    }
}
