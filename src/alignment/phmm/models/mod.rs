//! The definitions of the pHMMs and their components.

use crate::{
    alignment::phmm::{
        PhmmNumber,
        at_least_two::VecAtLeast2,
        components::{CorePhmm, EmissionParams, LayerParams},
        indexing::{AlnIndex, AlnIndexable, GetCore, GetLayer, GetLayerMut, GetMapping, GetModule},
        modules::{DomainModule, LocalModule, SemiLocalModule},
    },
    data::mappings::ByteIndexMap,
};
use std::{
    error::Error,
    fmt::{Debug, Display},
};

pub mod at_least_two;
pub mod components;

#[cfg(feature = "fuzzing")]
mod float_compare;

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
    pub fn from_parts(mapping: &'static ByteIndexMap<S>, core: CorePhmm<T, S>) -> GlobalPhmm<T, S> {
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
    ///
    /// ## Errors
    ///
    /// [`IncompatibleModuleError`] is returned if the length of `begin` or
    /// `end` doesn't match the length of `core`.
    #[inline]
    pub fn from_parts(
        mapping: &'static ByteIndexMap<S>, core: CorePhmm<T, S>, begin: LocalModule<T, S>, end: LocalModule<T, S>,
    ) -> Result<LocalPhmm<T, S>, IncompatibleModuleError> {
        if core.seq_len() != begin.semilocal_params.seq_len() || core.seq_len() != end.semilocal_params.seq_len() {
            return Err(IncompatibleModuleError);
        }

        Ok(Self {
            mapping,
            core,
            begin,
            end,
        })
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
    pub fn from_parts(
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
    ///
    /// ## Errors
    ///
    /// [`IncompatibleModuleError`] is returned if the length of `begin` or
    /// `end` doesn't match the length of `core`.
    #[inline]
    pub fn from_parts(
        mapping: &'static ByteIndexMap<S>, core: CorePhmm<T, S>, begin: SemiLocalModule<T>, end: SemiLocalModule<T>,
    ) -> Result<Self, IncompatibleModuleError> {
        if core.seq_len() != begin.seq_len() || core.seq_len() != end.seq_len() {
            return Err(IncompatibleModuleError);
        }

        Ok(Self {
            mapping,
            core,
            begin,
            end,
        })
    }

    /// Returns a reference to the [`ByteIndexMap`] used by the semilocal pHMM.
    #[inline]
    #[must_use]
    pub fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<T, S> {
    /// Gets the score for transitioning into a given [`AlnIndex`] from the
    /// [`SemiLocalModule`] at the beginning of the pHMM.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    pub(crate) fn get_begin_score(&self, index: impl AlnIndex) -> T {
        self.begin.get_score(index)
    }

    /// Gets the score for transitioning out of a given [`AlnIndex`] into the
    /// [`SemiLocalModule`] at the end of the pHMM.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    pub(crate) fn get_end_score(&self, index: impl AlnIndex) -> T {
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
    /// Converts a [`GlobalPhmm`] into a [`LocalPhmm`] using the provided
    /// `config`.
    ///
    /// ## Errors
    ///
    /// [`IncompatibleModuleError`] is returned if the length of either module
    /// is incorrect (for [`LocalConfig::Custom`]).
    #[inline]
    pub fn into_local_phmm(self, config: LocalConfig<T, S>) -> Result<LocalPhmm<T, S>, IncompatibleModuleError> {
        let (begin, end) = match config {
            LocalConfig::NoPenalty { background_emission } => (
                LocalModule::no_penalty(&self.core, background_emission.clone()),
                LocalModule::no_penalty(&self.core, background_emission),
            ),
            LocalConfig::Custom { begin, end } => (begin, end),
        };

        LocalPhmm::from_parts(self.mapping, self.core, begin, end)
    }

    /// Converts a [`GlobalPhmm`] into a [`DomainPhmm`] using the provided
    /// `config`.
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
        DomainPhmm::from_parts(self.mapping, self.core, begin, end)
    }

    /// Converts a [`GlobalPhmm`] into a [`SemiLocalPhmm`] using the provided
    /// `config`.
    ///
    /// ## Errors
    ///
    /// [`IncompatibleModuleError`] is returned if the length of either module
    /// is incorrect (for [`SemiLocalConfig::Custom`]).
    #[inline]
    pub fn into_semilocal_phmm(self, config: SemiLocalConfig<T>) -> Result<SemiLocalPhmm<T, S>, IncompatibleModuleError> {
        let (begin, end) = match config {
            SemiLocalConfig::NoPenalty => (
                SemiLocalModule::no_penalty(&self.core),
                SemiLocalModule::no_penalty(&self.core),
            ),
            SemiLocalConfig::Custom { begin, end } => (begin, end),
        };

        SemiLocalPhmm::from_parts(self.mapping, self.core, begin, end)
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

impl<T, const S: usize> GetCore<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        self.core().layers()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for GlobalPhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>> {
        self.core.layers_mut()
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
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        self.core().layers()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for DomainPhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>> {
        self.core.layers_mut()
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
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        self.core().layers()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>> {
        self.core.layers_mut()
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
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        self.core().layers()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for LocalPhmm<T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>> {
        self.core.layers_mut()
    }
}

impl<T, const S: usize> GetMapping<S> for GlobalPhmm<T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for LocalPhmm<T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for SemiLocalPhmm<T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for DomainPhmm<T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for &GlobalPhmm<T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for &LocalPhmm<T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for &SemiLocalPhmm<T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

impl<T, const S: usize> GetMapping<S> for &DomainPhmm<T, S> {
    #[inline]
    fn mapping(&self) -> &'static ByteIndexMap<S> {
        self.mapping
    }
}

/// An error representing an incompatible module in a pHMM.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct IncompatibleModuleError;

impl Display for IncompatibleModuleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "The pHMM's alignment modules are incompatible with the core pHMM!")
    }
}

impl Error for IncompatibleModuleError {}

/// An enum representing errors that can happen when working with pHMMs.
#[derive(Copy, Clone, Eq, PartialEq)]
pub enum InvalidModelError {
    /// The model has no layers
    EmptyModel,
    /// Not enough layers were specified
    TooFewLayers(usize),
    /// Alignment modules were incompatible with core pHMM
    IncompatibleModule,
}

impl Display for InvalidModelError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            InvalidModelError::EmptyModel => write!(f, "No layers were present in the pHMM!"),
            InvalidModelError::TooFewLayers(num) => {
                write!(f, "Too few layers were specified for the pHMM! At least {num} are required")
            }
            InvalidModelError::IncompatibleModule => {
                write!(f, "{IncompatibleModuleError}")
            }
        }
    }
}

impl Debug for InvalidModelError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for InvalidModelError {}
