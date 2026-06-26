use crate::{
    alignment::phmm::{
        PhmmNumber,
        indexing::{GetCore, GetCoreMut, GetLayer, GetLayerMut, GetModule, GetModuleMut, GetPartsMut, PhmmIndex},
        modules::{DomainModule, LocalModule, SemiLocalModule},
    },
    data::mappings::ByteIndexMap,
};

#[cfg(not(feature = "alignment-diagnostics"))]
#[doc(auto_cfg(hide(feature = "alignment-diagnostics")))]
mod components;
#[cfg(not(feature = "alignment-diagnostics"))]
#[doc(auto_cfg(hide(feature = "alignment-diagnostics")))]
pub use components::EmissionParams;
#[cfg(not(feature = "alignment-diagnostics"))]
#[doc(auto_cfg(hide(feature = "alignment-diagnostics")))]
pub(crate) use components::*;

#[cfg(feature = "alignment-diagnostics")]
mod components;
#[cfg(feature = "alignment-diagnostics")]
pub use components::*;

/// An implementation of a profile hidden Markov model (pHMM) for global
/// alignment (aligning a full sequence to a full model).
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct GlobalPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub(crate) mapping: &'static ByteIndexMap<S>,
    /// The model parameters
    pub(crate) core:    CorePhmm<T, S>,
}

impl<T, const S: usize> GlobalPhmm<T, S> {
    /// Creates a new [`GlobalPhmm`] from the specified mapping and
    /// [`CorePhmm`].
    #[inline]
    #[must_use]
    #[cfg(feature = "alignment-diagnostics")]
    pub fn new(mapping: &'static ByteIndexMap<S>, core: CorePhmm<T, S>) -> GlobalPhmm<T, S> {
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
    pub(crate) mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub(crate) core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    pub(crate) begin:   LocalModule<T, S>,
    /// The module for handling any bases after the core model.
    pub(crate) end:     LocalModule<T, S>,
}

impl<T, const S: usize> LocalPhmm<T, S> {
    /// Creates a new [`LocalPhmm`] from the specified parts.
    #[inline]
    #[must_use]
    #[cfg(feature = "alignment-diagnostics")]
    pub fn new(
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
    pub(crate) mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub(crate) core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    pub(crate) begin:   DomainModule<T, S>,
    /// The module for handling any bases after the core model.
    pub(crate) end:     DomainModule<T, S>,
}

impl<T, const S: usize> DomainPhmm<T, S> {
    /// Creates a new [`DomainPhmm`] from the specified parts.
    #[inline]
    #[must_use]
    #[cfg(feature = "alignment-diagnostics")]
    pub fn new(
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
    pub(crate) mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub(crate) core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    pub(crate) begin:   SemiLocalModule<T>,
    /// The module for handling any bases after the core model.
    pub(crate) end:     SemiLocalModule<T>,
}

impl<T, const S: usize> SemiLocalPhmm<T, S> {
    /// Creates a new [`SemiLocalPhmm`] from the specified parts.
    #[inline]
    #[must_use]
    #[cfg(feature = "alignment-diagnostics")]
    pub fn new(
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
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
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
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
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
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
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
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.core().split_first_layer()
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
