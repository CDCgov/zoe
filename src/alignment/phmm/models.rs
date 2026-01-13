use crate::{
    alignment::phmm::{
        PhmmError, PhmmNumber, PhmmState,
        indexing::{LastMatch, PhmmIndex, PhmmIndexRange, PhmmIndexable},
        modules::{DomainModule, LocalModule, SemiLocalModule},
    },
    data::ByteIndexMap,
};
use std::ops::{Index, IndexMut};

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
    /// At least two layers are required, otherwise [`PhmmError::TooFewLayers`]
    /// is returned.
    #[inline]
    #[allow(dead_code)]
    pub fn new(layers: Vec<LayerParams<T, S>>) -> Result<Self, PhmmError> {
        if layers.len() >= 2 {
            Ok(CorePhmm(layers))
        } else {
            Err(PhmmError::TooFewLayers(2))
        }
    }

    /// Creates a new [`CorePhmm`] from a `Vec` of the parameters, without
    /// checking that the number of layers is valid (the corresponding reference
    /// length should be at least 1).
    #[inline]
    #[must_use]
    pub(crate) fn new_unchecked(layers: Vec<LayerParams<T, S>>) -> Self {
        CorePhmm(layers)
    }

    /// Returns a reference to the layer parameters stored in the core pHMM.
    #[inline]
    #[must_use]
    pub fn layers(&self) -> &[LayerParams<T, S>] {
        self.0.as_slice()
    }

    /// Returns a mutable reference to the layer parameters stored in the core
    /// pHMM.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    pub fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.0.as_mut_slice()
    }

    /// Returns a mutable reference to the vector of layer parameters stored in
    /// the core pHMM.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    pub fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        &mut self.0
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

/// A trait providing accessors to the modules at the beginning and end of a
/// pHMM.
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
    type Begin;
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

impl<T, const S: usize> GetModule for LocalPhmm<T, S> {
    type Begin = LocalModule<T, S>;
    type End = LocalModule<T, S>;

    #[inline]
    fn begin(&self) -> &Self::Begin {
        &self.begin
    }

    #[inline]
    fn begin_mut(&mut self) -> &mut Self::Begin {
        &mut self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        &self.end
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        &mut self.end
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
    fn begin_mut(&mut self) -> &mut Self::Begin {
        &mut self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        &self.end
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        &mut self.end
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
    fn begin_mut(&mut self) -> &mut Self::Begin {
        &mut self.begin
    }

    #[inline]
    fn end(&self) -> &Self::End {
        &self.end
    }

    #[inline]
    fn end_mut(&mut self) -> &mut Self::End {
        &mut self.end
    }
}

/// A trait providing accessors to the layers of a pHMM.
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
    /// Returns a reference to the [`CorePhmm`] holding the core parameters.
    #[must_use]
    fn core(&self) -> &CorePhmm<T, S>;

    /// Returns a mutable reference to the [`CorePhmm`] holding the core
    /// parameters.
    #[must_use]
    #[allow(dead_code)]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S>;

    /// Retrieves a slice of the layers contained within the core pHMM.
    #[inline]
    #[must_use]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.core().layers()
    }

    /// Retrieves a mutable slice of the layers contained within the core pHMM.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>] {
        self.core_mut().layers_mut()
    }

    /// Returns a mutable reference to the vector of layer parameters stored in
    /// the core pHMM.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn layers_mut_vec(&mut self) -> &mut Vec<LayerParams<T, S>> {
        self.core_mut().layers_mut_vec()
    }

    /// Get a layer from within the core pHMM.
    ///
    /// Although there is no actual layer for the `End` state, for readability
    /// we let `End` be synonymous with `LastMatch` since `End` emphasizes it is
    /// the last layer.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layer(&self, j: impl PhmmIndex) -> &LayerParams<T, S> {
        if j.is_end() {
            self.get_layer(LastMatch)
        } else {
            &self.layers()[self.get_dp_index(j)]
        }
    }

    /// Get a mutable reference to a layer from within the core pHMM.
    ///
    /// Although there is no actual layer for the `End` state, for readability
    /// we let `End` be synonymous with `LastMatch` since `End` emphasizes it is
    /// the last layer.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layer_mut(&mut self, j: impl PhmmIndex) -> &mut LayerParams<T, S> {
        if j.is_end() {
            self.get_layer_mut(LastMatch)
        } else {
            let idx = self.get_dp_index(j);
            &mut self.layers_mut()[idx]
        }
    }

    /// Get a range of layers from within the core pHMM.
    ///
    /// ## Panics
    ///
    /// If any of the indices are out of bounds, this will panic. Particularly,
    /// if the range is end-inclusive and ends with `End` (e.g., `..=End`), this
    /// will panic. This is different behavior than [`get_layer`].
    ///
    /// [`get_layer`]: GetLayer::get_layer
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layers(&self, range: impl PhmmIndexRange) -> &[LayerParams<T, S>] {
        &self.layers()[self.get_dp_range(range)]
    }

    /// Get a range of mutable layers from within the core pHMM.
    ///
    /// ## Panics
    ///
    /// If any of the indices are out of bounds, this will panic. Particularly,
    /// if the range is end-inclusive and ends with `End` (e.g., `..=End`), this
    /// will panic. This is different behavior than [`get_layer_mut`].
    ///
    /// [`get_layer_mut`]: GetLayer::get_layer_mut
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_layers_mut(&mut self, range: impl PhmmIndexRange) -> &mut [LayerParams<T, S>] {
        let range = self.get_dp_range(range);
        &mut self.layers_mut()[range]
    }

    /// Gets mutable references to two distinct layers within the core pHMM.
    ///
    /// Although there is no actual layer for the `End` state, for readability
    /// we let `End` be synonymous with `LastMatch` since `End` emphasizes it is
    /// the last layer.
    ///
    /// ## Panics
    ///
    /// If the requested indices correspond to the same dynamic programming
    /// index or are out of bounds, this will panic.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    fn get_two_layers_mut(
        &mut self, j1: impl PhmmIndex, j2: impl PhmmIndex,
    ) -> (&mut LayerParams<T, S>, &mut LayerParams<T, S>) {
        let j1 = if j1.is_end() {
            self.get_dp_index(LastMatch)
        } else {
            self.get_dp_index(j1)
        };
        let j2 = if j2.is_end() {
            self.get_dp_index(LastMatch)
        } else {
            self.get_dp_index(j2)
        };
        let (j1, j2) = (std::cmp::min(j1, j2), std::cmp::max(j1, j2));
        // Split into [0, j1] and (j1, len). This is the same as [0, j1+1) and
        // [j1+1, len). This will fail if both j1 and j2 are out of bounds
        let (s1, s2) = self.layers_mut().split_at_mut(j1 + 1);
        // This will fail if j2 is out of bounds or equal to j1
        (&mut s1[j1], &mut s2[j2 - (j1 + 1)])
    }
}

impl<T, const S: usize> GetLayer<T, S> for CorePhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self
    }

    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        self
    }
}

impl<T, const S: usize> GetLayer<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }

    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        &mut self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }

    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        &mut self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }

    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        &mut self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }

    #[inline]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S> {
        &mut self.core
    }
}
