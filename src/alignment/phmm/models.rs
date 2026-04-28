use crate::{
    alignment::phmm::{
        InvalidModelError, PhmmError, PhmmNumber, PhmmState,
        indexing::{PhmmIndex, PhmmIndexRange, PhmmIndexable},
        modules::{DomainModule, LocalModule, SemiLocalModule},
    },
    data::mappings::ByteIndexMap,
};
use std::{
    ops::{Index, IndexMut},
    slice::GetDisjointMutError,
};

/// The transition probabilities for a layer of the pHMM.
///
/// The parameters are stored in log space (see
/// [here](crate::alignment::phmm#log-space-parameters) for more details). They
/// are arranged in the format:
/// ```text
/// [[m->m, d->m, i->m]
///  [m->d, d->d, i->d]
///  [m->i, d->i, i->i]]
/// ```
/// This means the first element is the transition probabilities into the match
/// state, the second element is the probabilities into the delete state, and
/// the third element is the probabilities into the insert state.
///
/// Note that this is a different layout than
/// [SAM](https://tr.soe.ucsc.edu/sites/default/files/technical-reports/UCSC-CRL-96-22.pdf).
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// access this struct.
///
/// </div>
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) struct TransitionParams<T>(pub(crate) [[T; 3]; 3]);

impl<T: Copy> TransitionParams<T> {
    /// Retrieves an array containing the three parameters for exiting `state`
    /// and moving to the next layer.
    #[inline]
    #[must_use]
    pub fn exiting_params(&self, state: PhmmState) -> [T; 3] {
        [
            self[(state, PhmmState::Match)],
            self[(state, PhmmState::Delete)],
            self[(state, PhmmState::Insert)],
        ]
    }

    /// Retrieves an array containing the three parameters for entering `state`
    /// from the previous layer.
    #[inline]
    #[must_use]
    pub fn entering_params(&self, state: PhmmState) -> [T; 3] {
        self[state]
    }

    /// Retrieves the inner array from the [`TransitionParams`]. The layout of
    /// this array is subject to change. Consider indexing directly into the
    /// [`TransitionParams`] instead.
    #[inline]
    #[must_use]
    #[cfg_attr(feature = "dev-phmm-regression", visibility::make(pub))]
    pub fn as_array(&self) -> &[[T; 3]; 3] {
        &self.0
    }

    /// Constructs a [`TransitionParams`] from an array. The interpretation of
    /// this array is subject to change. Consider using [`Default`] and then
    /// mutating each entry.
    #[inline]
    #[must_use]
    #[cfg_attr(feature = "dev-phmm-regression", visibility::make(pub))]
    pub fn from_array(arr: [[T; 3]; 3]) -> Self {
        Self(arr)
    }
}

impl<T: PhmmNumber> Default for TransitionParams<T> {
    /// Initializes the transition parameters so that all transitions are
    /// probability zero.
    #[inline]
    fn default() -> Self {
        // Zero probability is infinity in negative log space
        Self([[T::INFINITY; 3]; 3])
    }
}

impl<T> Index<(PhmmState, PhmmState)> for TransitionParams<T> {
    type Output = T;

    /// Retrieves the transition parameter corresponding to moving from state
    /// `index.0` to state `index.1`.
    #[inline]
    fn index(&self, index: (PhmmState, PhmmState)) -> &Self::Output {
        &self.0[usize::from(index.1)][usize::from(index.0)]
    }
}

impl<T> IndexMut<(PhmmState, PhmmState)> for TransitionParams<T> {
    /// Retrieves a mutable reference to the transition parameter corresponding
    /// to moving from state `index.0` to state `index.1`.
    #[inline]
    fn index_mut(&mut self, index: (PhmmState, PhmmState)) -> &mut Self::Output {
        &mut self.0[usize::from(index.1)][usize::from(index.0)]
    }
}

impl<T> Index<PhmmState> for TransitionParams<T> {
    type Output = [T; 3];

    /// Retrieves the transition parameters for moving into state `index` from
    /// the match, delete, and insert states respectively.
    #[inline]
    fn index(&self, index: PhmmState) -> &Self::Output {
        &self.0[usize::from(index)]
    }
}

impl<T> IndexMut<PhmmState> for TransitionParams<T> {
    /// Retrieves a mutable reference to the transition parameters for moving
    /// into state `index` from the match, delete, and insert states
    /// respectively.
    #[inline]
    fn index_mut(&mut self, index: PhmmState) -> &mut Self::Output {
        &mut self.0[usize::from(index)]
    }
}

/// A set of emission probabilities.
///
/// The parameters are stored in log space (see
/// [here](crate::alignment::phmm#log-space-parameters) for more details).
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct EmissionParams<T, const S: usize>([T; S]);

impl<T, const S: usize> EmissionParams<T, S> {
    /// Returns a [`EmissionParams`] object from an array.
    ///
    /// The elements in the array should be in the same order as the keys in the
    /// [`ByteIndexMap`] used by the pHMM.
    #[inline]
    #[must_use]
    pub fn from_array(arr: [T; S]) -> Self {
        Self(arr)
    }

    /// Returns the parameters as a slice.
    ///
    /// The order of the elements matches the order of the keys in the
    /// [`ByteIndexMap`] used by the pHMM.
    #[inline]
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.0
    }

    /// Returns a reference to the underlying array storing the parameters.
    ///
    /// The order of the elements matches the order of the keys in the
    /// [`ByteIndexMap`] used by the pHMM.
    #[inline]
    #[must_use]
    pub fn as_array(&self) -> &[T; S] {
        &self.0
    }

    /// Retrieves an iterator over the parameters.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.0.iter()
    }
}

impl<T: PhmmNumber, const S: usize> EmissionParams<T, S> {
    /// Generates the emission parameters for a uniform distribution.
    #[inline]
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn uniform() -> Self {
        Self([T::from_prob(1.0f64 / (S as f64)); S])
    }
}

impl<T: PhmmNumber, const S: usize> Default for EmissionParams<T, S> {
    /// Initializes the emission parameters so that all residues are probability
    /// zero.
    #[inline]
    fn default() -> Self {
        // Zero probability is infinity in negative log space
        Self([T::INFINITY; S])
    }
}

impl<'a, T, const S: usize> IntoIterator for &'a EmissionParams<T, S> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<T, const S: usize> Index<usize> for EmissionParams<T, S> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

/// The parameters for a single layer of the pHMM.
///
/// The parameters are stored in log space (see
/// [here](crate::alignment::phmm#log-space-parameters) for more details). Also
/// see [`TransitionParams`], [`EmissionParams`], and [`CorePhmm`].
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// access this struct.
///
/// </div>
#[derive(Clone, Eq, PartialEq, Debug)]
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) struct LayerParams<T, const S: usize> {
    pub transition:      TransitionParams<T>,
    pub emission_match:  EmissionParams<T, S>,
    pub emission_insert: EmissionParams<T, S>,
}

impl<T: PhmmNumber, const S: usize> Default for LayerParams<T, S> {
    #[inline]
    fn default() -> Self {
        Self {
            transition:      TransitionParams::default(),
            emission_match:  EmissionParams::default(),
            emission_insert: EmissionParams::default(),
        }
    }
}

/// The core profile hidden Markov model (pHMM) used by [`GlobalPhmm`],
/// [`LocalPhmm`], [`DomainPhmm`], and [`SemiLocalPhmm`].
///
/// This includes the layers of the pHMM without any modules at the beginning or
/// end. This struct guarantees that at least two layers are present
/// (corresponding to a reference length of one).
///
/// Internally, each element stores:
///
/// - The transition probabilities into the insert state within the current
///   layer, as well as the insert state's emission probabilities
/// - The transition probabilities from the current layer's states into the next
///   layer's match and delete states, as well as the next layer's emission
///   probabilities for the match state
///
/// The first element's match state is the BEGIN state, while the last element's
/// transition probabilities reflect the probabilities entering the END state.
///
/// Note that this is different from
/// [SAM](https://tr.soe.ucsc.edu/sites/default/files/technical-reports/UCSC-CRL-96-22.pdf),
/// which instead has:
///
/// - The emission probabilities for the insert and match states are stored in
///   that layer
/// - The transition probabilities from the previous layer to the current layer
///   are stored
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// access this struct.
///
/// </div>
#[derive(Clone, Eq, PartialEq, Debug)]
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) struct CorePhmm<T, const S: usize>(Vec<LayerParams<T, S>>);

impl<T, const S: usize> CorePhmm<T, S> {
    /// Create a new [`CorePhmm`] from a `Vec` of the parameters.
    ///
    /// ## Errors
    ///
    /// At least two layers are required, otherwise [`InvalidModel`] with cause
    /// [`TooFewLayers`] is returned.
    ///
    /// [`InvalidModel`]: PhmmError::InvalidModel
    /// [`TooFewLayers`]:
    ///     crate::alignment::phmm::errors::InvalidModelError::TooFewLayers
    #[inline]
    #[allow(dead_code)]
    pub fn new(layers: Vec<LayerParams<T, S>>) -> Result<Self, PhmmError> {
        if layers.len() >= 2 {
            Ok(CorePhmm(layers))
        } else {
            Err(InvalidModelError::TooFewLayers(2).into())
        }
    }

    /// Creates a new [`CorePhmm`] from a `Vec` of the parameters, without
    /// checking that the number of layers is valid.
    ///
    /// ## Validity
    ///
    /// The length of `layers` must be at least 2, corresponding to a reference
    /// length of at least 1.
    #[inline]
    #[must_use]
    pub(crate) fn new_unchecked(layers: Vec<LayerParams<T, S>>) -> Self {
        CorePhmm(layers)
    }
}

/// An implementation of a profile hidden Markov model (pHMM) for global
/// alignment (aligning a full sequence to a full model).
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct GlobalPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    mapping: &'static ByteIndexMap<S>,
    /// The model parameters
    core:    CorePhmm<T, S>,
}

impl<T, const S: usize> GlobalPhmm<T, S> {
    /// Creates a new [`GlobalPhmm`] from the specified mapping and
    /// [`CorePhmm`].
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    #[inline]
    #[must_use]
    #[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
    pub(crate) fn new(mapping: &'static ByteIndexMap<S>, core: CorePhmm<T, S>) -> GlobalPhmm<T, S> {
        Self { mapping, core }
    }

    /// Returns a reference to the [`ByteIndexMap`] used by the global pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

/// An implementation of a profile hidden Markov model (pHMM) for local
/// alignment (aligning a subsequence to a submodel).
///
/// This is created from a [`GlobalPhmm`] using [`into_local_phmm`]. Two
/// [`LocalModule`] modules are added to either end which can match arbitrarily
/// many bases at the beginning or end of the sequence, and can skip arbitrarily
/// many states at the beginning or end of the pHMM.
///
/// [`into_local_phmm`]: GlobalPhmm::into_local_phmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct LocalPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    begin:   LocalModule<T, S>,
    /// The module for handling any bases after the core model.
    end:     LocalModule<T, S>,
}

impl<T, const S: usize> LocalPhmm<T, S> {
    /// Creates a new [`LocalPhmm`] from the specified parts.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    #[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
    pub(crate) fn new(
        mapping: &'static ByteIndexMap<S>, core: CorePhmm<T, S>, begin: LocalModule<T, S>, end: LocalModule<T, S>,
    ) -> LocalPhmm<T, S> {
        Self {
            mapping,
            core,
            begin,
            end,
        }
    }

    /// Returns a reference to the [`ByteIndexMap`] used by the local pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

/// An implementation of a profile hidden Markov model (pHMM) for domain
/// alignment (aligning a subsequence to a full model).
///
/// This is created from a [`GlobalPhmm`] using [`into_domain_phmm`]. Two
/// [`DomainModule`] modules are added to either end which can match arbitrarily
/// many bases at the beginning or end of the sequence.
///
/// [`into_domain_phmm`]: GlobalPhmm::into_domain_phmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct DomainPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    begin:   DomainModule<T, S>,
    /// The module for handling any bases after the core model.
    end:     DomainModule<T, S>,
}

impl<T, const S: usize> DomainPhmm<T, S> {
    /// Creates a new [`DomainPhmm`] from the specified parts.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    #[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
    pub(crate) fn new(
        mapping: &'static ByteIndexMap<S>, core: CorePhmm<T, S>, begin: DomainModule<T, S>, end: DomainModule<T, S>,
    ) -> Self {
        Self {
            mapping,
            core,
            begin,
            end,
        }
    }

    /// Returns a reference to the [`ByteIndexMap`] used by the domain pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

/// An implementation of a profile hidden Markov model (pHMM) for semilocal
/// alignment (aligning a full sequence to a submodel).
///
/// This is created from a [`GlobalPhmm`] using [`into_semilocal_phmm`]. Two
/// [`SemiLocalModule`] modules are added to either end which can skip
/// arbitrarily many states at the beginning or end of the pHMM.
///
/// [`into_semilocal_phmm`]: GlobalPhmm::into_semilocal_phmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SemiLocalPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    begin:   SemiLocalModule<T>,
    /// The module for handling any bases after the core model.
    end:     SemiLocalModule<T>,
}

impl<T, const S: usize> SemiLocalPhmm<T, S> {
    /// Creates a new [`SemiLocalPhmm`] from the specified parts.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    #[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
    pub(crate) fn new(
        mapping: &'static ByteIndexMap<S>, core: CorePhmm<T, S>, begin: SemiLocalModule<T>, end: SemiLocalModule<T>,
    ) -> Self {
        Self {
            mapping,
            core,
            begin,
            end,
        }
    }

    /// Returns a reference to the [`ByteIndexMap`] used by the semilocal pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<T, S> {
    /// Gets the score for transitioning into a given [`PhmmIndex`] from the
    /// [`SemiLocalModule`] at the beginning of the pHMM.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    pub(crate) fn get_begin_score(&self, index: impl PhmmIndex) -> T {
        self.begin.get_score(index)
    }

    /// Gets the score for transitioning out of a given [`PhmmIndex`] into the
    /// [`SemiLocalModule`] at the end of the pHMM.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    pub(crate) fn get_end_score(&self, index: impl PhmmIndex) -> T {
        self.end.get_score(index)
    }
}

/// Options for how to construct a [`LocalPhmm`] from a [`GlobalPhmm`].
#[non_exhaustive]
pub enum LocalConfig<T, const S: usize> {
    /// Does not penalize any transitions in the [`LocalModule`], instead solely
    /// penalizing the emissions for any insertions.
    NoPenalty { background_emission: EmissionParams<T, S> },
    /// Use a custom [`LocalModule`] for the begin and end.
    Custom {
        begin: LocalModule<T, S>,
        end:   LocalModule<T, S>,
    },
}

/// Options for how to construct a [`DomainPhmm`] from a [`GlobalPhmm`].
#[non_exhaustive]
pub enum DomainConfig<T, const S: usize> {
    /// Does not penalize any transitions in the [`DomainModule`], instead
    /// solely penalizing the emissions for any insertions.
    NoPenalty { background_emission: EmissionParams<T, S> },
    /// Use a custom [`DomainModule`] for the begin and end.
    Custom {
        begin: DomainModule<T, S>,
        end:   DomainModule<T, S>,
    },
}

/// Options for how to construct a [`SemiLocalPhmm`] from a [`GlobalPhmm`].
#[non_exhaustive]
pub enum SemiLocalConfig<T> {
    /// Does not penalize any transitions in the [`SemiLocalModule`].
    NoPenalty,
    /// Use a custom [`SemiLocalModule`] for the begin and end.
    Custom {
        begin: SemiLocalModule<T>,
        end:   SemiLocalModule<T>,
    },
}

impl<T: PhmmNumber, const S: usize> GlobalPhmm<T, S> {
    /// Creates a [`LocalPhmm`] from a [`GlobalPhmm`].
    ///
    /// The method to use for defining the local alignment behavior is specified
    /// with `config`.
    #[inline]
    #[must_use]
    pub fn into_local_phmm(self, config: LocalConfig<T, S>) -> LocalPhmm<T, S> {
        let (begin, end) = match config {
            LocalConfig::NoPenalty { background_emission } => (
                LocalModule::no_penalty(&self.core, background_emission.clone()),
                LocalModule::no_penalty(&self.core, background_emission),
            ),
            LocalConfig::Custom { begin, end } => (begin, end),
        };
        LocalPhmm {
            mapping: self.mapping,
            core: self.core,
            begin,
            end,
        }
    }

    /// Creates a [`DomainPhmm`] from a [`GlobalPhmm`].
    ///
    /// The method to use for defining the local alignment behavior is specified
    /// with `config`.
    #[inline]
    #[must_use]
    pub fn into_domain_phmm(self, config: DomainConfig<T, S>) -> DomainPhmm<T, S> {
        let (begin, end) = match config {
            DomainConfig::NoPenalty { background_emission } => (
                DomainModule::no_penalty(background_emission.clone()),
                DomainModule::no_penalty(background_emission),
            ),
            DomainConfig::Custom { begin, end } => (begin, end),
        };
        DomainPhmm {
            mapping: self.mapping,
            core: self.core,
            begin,
            end,
        }
    }

    /// Converts a [`GlobalPhmm`] into a [`SemiLocalPhmm`] using the provided
    /// `config`.
    #[inline]
    #[must_use]
    pub fn into_semilocal_phmm(self, config: SemiLocalConfig<T>) -> SemiLocalPhmm<T, S> {
        let (begin, end) = match config {
            SemiLocalConfig::NoPenalty => (
                SemiLocalModule::no_penalty(&self.core),
                SemiLocalModule::no_penalty(&self.core),
            ),
            SemiLocalConfig::Custom { begin, end } => (begin, end),
        };
        SemiLocalPhmm {
            mapping: self.mapping,
            core: self.core,
            begin,
            end,
        }
    }
}

/// A trait providing read-only access to the modules at the beginning and end
/// of a pHMM.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// use this trait.
///
/// </div>
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) trait GetModule {
    /// The type of the module at the beginning of the pHMM.
    type Begin;
    /// The type of the module at the end of the pHMM.
    type End;

    /// Returns a reference to the module at the start of the pHMM.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    #[must_use]
    fn begin(&self) -> &Self::Begin;

    /// Returns a reference to the module at the end of the pHMM.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    #[must_use]
    fn end(&self) -> &Self::End;
}

/// A trait providing mutable access to the modules at the beginning and end of
/// a pHMM.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// use this trait.
///
/// </div>
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) trait GetModuleMut: GetModule {
    /// Returns a mutable reference to the module at the start of the pHMM.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    #[must_use]
    #[allow(dead_code)]
    fn begin_mut(&mut self) -> &mut Self::Begin;

    /// Returns a mutable reference to the module at the end of the pHMM.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
    /// to use this method.
    ///
    /// </div>
    #[must_use]
    #[allow(dead_code)]
    fn end_mut(&mut self) -> &mut Self::End;
}

// This is a separate trait from GetLayer in order to prevent core from being
// called on a CorePhmm, which is an easy way to have infinite recursion in an
// implementation.

/// A trait providing read-only access to the [`CorePhmm`] within a larger pHMM.
pub(crate) trait GetCore<T, const S: usize> {
    /// Returns a reference to the [`CorePhmm`] holding the core parameters.
    #[must_use]
    fn core(&self) -> &CorePhmm<T, S>;
}

// This is a separate trait from GetLayerMut in order to prevent core_mut from
// being called on a CorePhmm, which is an easy way to have infinite recursion
// in an implementation.

/// A trait providing read-only access to the [`CorePhmm`] within a larger pHMM.
pub(crate) trait GetCoreMut<T, const S: usize> {
    /// Returns a mutable reference to the [`CorePhmm`] holding the core
    /// parameters.
    #[must_use]
    #[allow(dead_code)]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S>;
}

/// A trait providing read-only accessors to the layers of a pHMM.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml`
/// to use this trait.
///
/// </div>
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) trait GetLayer<T, const S: usize>: PhmmIndexable {
    /// Retrieves a slice of the layers contained within the core pHMM.
    ///
    /// This slice will be at least 2 in length.
    #[must_use]
    fn layers(&self) -> &[LayerParams<T, S>];

    /// Returns the last layer, as well as all previous layers.
    ///
    /// This is an infallible version of `model.layers().split_last()`.
    #[must_use]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]);

    /// Returns a reference to the parameters for the layer containing the BEGIN
    /// state.
    ///
    /// This is an infallible version of `model.get_layer(Begin)`.
    #[inline]
    #[must_use]
    fn begin_layer(&self) -> &LayerParams<T, S> {
        &self.layers()[0]
    }

    /// Returns a reference to the parameters for the layer containing the first
    /// match state.
    ///
    /// This is an infallible version of `model.get_layer(FirstMatch)`.
    #[inline]
    #[must_use]
    fn first_match(&self) -> &LayerParams<T, S> {
        &self.layers()[1]
    }

    /// Returns a reference to the parameters for the layer containing the last
    /// match state.
    ///
    /// This is an infallible version of `model.get_layer(LastMatch)`.
    #[inline]
    #[must_use]
    fn last_match(&self) -> &LayerParams<T, S> {
        &self.layers()[self.layers().len() - 1]
    }

    /// Gets a layer from within the core pHMM.
    ///
    /// This returns `None` if the index is out of bounds or [`End`] (since
    /// there is no layer corresponding to the END state).
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    #[inline]
    #[must_use]
    fn get_layer(&self, j: impl PhmmIndex) -> Option<&LayerParams<T, S>> {
        self.layers().get(self.get_dp_index(j))
    }

    /// Gets a range of layers from within the core pHMM.
    ///
    /// If any of the indices are out of bounds, this will return `None`.
    /// Particularly, if the range is end-inclusive and ends with `End` (e.g.,
    /// `..=End`), this will return `None`.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layers(&self, range: impl PhmmIndexRange) -> Option<&[LayerParams<T, S>]> {
        self.layers().get(self.get_dp_range(range))
    }
}

/// A trait providing mutable accessors to the layers of a pHMM.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// use this trait.
///
/// </div>
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) trait GetLayerMut<T, const S: usize>: GetLayer<T, S> {
    /// Retrieves a mutable slice of the layers contained within the core pHMM.
    #[must_use]
    #[allow(dead_code)]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>];

    /// Returns a mutable reference to the vector of layer parameters stored in
    /// the core pHMM.
    #[must_use]
    #[allow(dead_code)]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>>;

    /// Returns a reference to the parameters for the layer containing the BEGIN
    /// state.
    ///
    /// This is an infallible version of `model.get_layer_mut(Begin)`.
    #[inline]
    #[must_use]
    fn begin_layer_mut(&mut self) -> &mut LayerParams<T, S> {
        &mut self.layers_mut()[0]
    }

    /// Returns a reference to the parameters for the layer containing the first
    /// match state.
    ///
    /// This is an infallible version of `model.get_layer_mut(FirstMatch)`.
    #[inline]
    #[must_use]
    fn first_match_mut(&mut self) -> &mut LayerParams<T, S> {
        &mut self.layers_mut()[1]
    }

    /// Returns a reference to the parameters for the layer containing the last
    /// match state.
    ///
    /// This is an infallible version of `model.get_layer_mut(LastMatch)`.
    #[inline]
    #[must_use]
    fn last_match_mut(&mut self) -> &mut LayerParams<T, S> {
        let idx = self.layers().len() - 1;
        &mut self.layers_mut()[idx]
    }

    /// Gets a mutable reference to a layer from within the core pHMM.
    ///
    /// This returns `None` if the index is out of bounds or [`End`] (since
    /// there is no layer corresponding to the END state).
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layer_mut(&mut self, j: impl PhmmIndex) -> Option<&mut LayerParams<T, S>> {
        let idx = self.get_dp_index(j);
        self.layers_mut().get_mut(idx)
    }

    /// Get a range of mutable layers from within the core pHMM.
    ///
    /// If any of the indices are out of bounds, this will return `None`.
    /// Particularly, if the range is end-inclusive and ends with `End` (e.g.,
    /// `..=End`), this will return `None`.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layers_mut(&mut self, range: impl PhmmIndexRange) -> Option<&mut [LayerParams<T, S>]> {
        let range = self.get_dp_range(range);
        self.layers_mut().get_mut(range)
    }

    /// Gets mutable references to two distinct layers within the core pHMM.
    ///
    /// ## Errors
    ///
    /// - [`IndexOutOfBounds`] if either index is out of bounds or [`End`]
    ///   (since there is no layer corresponding to the END state)
    /// - [`OverlappingIndices`] if `j1` and `j2` are the same index
    ///
    /// [`IndexOutOfBounds`]: GetDisjointMutError::IndexOutOfBounds
    /// [`End`]: crate::alignment::phmm::indexing::End
    /// [`OverlappingIndices`]: GetDisjointMutError::OverlappingIndices
    #[inline]
    #[allow(dead_code)]
    fn get_two_layers_mut(
        &mut self, j1: impl PhmmIndex, j2: impl PhmmIndex,
    ) -> Result<(&mut LayerParams<T, S>, &mut LayerParams<T, S>), GetDisjointMutError> {
        let j1 = self.get_dp_index(j1);
        let j2 = self.get_dp_index(j2);

        let [l1, l2] = self.layers_mut().get_disjoint_mut([j1, j2])?;

        Ok((l1, l2))
    }
}

/// A trait providing simultaneous mutable accessors to the parts of a pHMM.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// use this trait.
///
/// </div>
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) trait GetPartsMut<T, const S: usize>: GetModule + GetLayer<T, S> {
    /// Returns simultaneous mutable references to the core pHMM, the begin
    /// module, and the end module.
    ///
    /// Calling the individual mutable accessors and simultaneously using them
    /// is not allowed by the borrow checker, hence this function.
    #[must_use]
    fn parts_mut(&mut self) -> (&mut CorePhmm<T, S>, &mut Self::Begin, &mut Self::End);
}

impl<T, const S: usize> GetModule for LocalPhmm<T, S> {
    type Begin = LocalModule<T, S>;
    type End = LocalModule<T, S>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        &self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        &self.end
    }
}

impl<T, const S: usize> GetModule for DomainPhmm<T, S> {
    type Begin = DomainModule<T, S>;
    type End = DomainModule<T, S>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        &self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        &self.end
    }
}

impl<T, const S: usize> GetModule for SemiLocalPhmm<T, S> {
    type Begin = SemiLocalModule<T>;
    type End = SemiLocalModule<T>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        &self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        &self.end
    }
}

impl<T, const S: usize> GetModuleMut for LocalPhmm<T, S> {
    #[inline]
    fn begin_mut(&mut self) -> &mut Self::Begin {
        &mut self.begin
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        &mut self.end
    }
}

impl<T, const S: usize> GetModuleMut for DomainPhmm<T, S> {
    #[inline]
    fn begin_mut(&mut self) -> &mut Self::Begin {
        &mut self.begin
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        &mut self.end
    }
}

impl<T, const S: usize> GetModuleMut for SemiLocalPhmm<T, S> {
    #[inline]
    fn begin_mut(&mut self) -> &mut Self::Begin {
        &mut self.begin
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        &mut self.end
    }
}

impl<T, const S: usize> GetLayer<T, S> for CorePhmm<T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.0.as_slice()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.0.as_slice().split_last().expect("A CorePhmm has at least two layers")
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for CorePhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.0.as_mut_slice()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        &mut self.0
    }
}

impl<T, const S: usize> GetCore<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetCoreMut<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        &mut self.core
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core.layers_mut()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core.layers_mut_vec()
    }
}

impl<T, const S: usize> GetCore<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetCoreMut<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        &mut self.core
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core.layers_mut()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core.layers_mut_vec()
    }
}

impl<T, const S: usize> GetCore<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetCoreMut<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        &mut self.core
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core.layers_mut()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core.layers_mut_vec()
    }
}

impl<T, const S: usize> GetCore<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    #[inline]
    fn split_last_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_last_layer()
    }
}

impl<T, const S: usize> GetCoreMut<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        &mut self.core
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core.layers_mut()
    }

    #[inline]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core.layers_mut_vec()
    }
}

impl<T, const S: usize> GetPartsMut<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn parts_mut(&mut self) -> (&mut CorePhmm<T, S>, &mut Self::Begin, &mut Self::End) {
        (&mut self.core, &mut self.begin, &mut self.end)
    }
}

impl<T, const S: usize> GetPartsMut<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn parts_mut(&mut self) -> (&mut CorePhmm<T, S>, &mut Self::Begin, &mut Self::End) {
        (&mut self.core, &mut self.begin, &mut self.end)
    }
}

impl<T, const S: usize> GetPartsMut<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn parts_mut(&mut self) -> (&mut CorePhmm<T, S>, &mut Self::Begin, &mut Self::End) {
        (&mut self.core, &mut self.begin, &mut self.end)
    }
}
