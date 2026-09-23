//! The definitions of the pHMMs and their components.

use crate::{
    alignment::phmm::{
        PhmmNumber,
        at_least_two::VecAtLeast2,
        components::{CorePhmm, EmissionParams, LayerParams},
        indexing::{AlnIndexable, GetCore, GetLayer, GetLayerMut, GetMapping, GetModule},
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
///
/// ## Parameters
///
/// - `'a`: The lifetime of the alphabet. This can be `'static` for compile-time
///   alphabets.
/// - `T`: The numeric type used to store the parameters.
/// - `S`: The size of the alphabet.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct GlobalPhmm<'a, T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub(crate) mapping: &'a ByteIndexMap<S>,
    /// The model parameters
    pub(crate) core:    CorePhmm<T, S>,
}

impl<'a, T, const S: usize> GlobalPhmm<'a, T, S> {
    /// Creates a new [`GlobalPhmm`] from the specified mapping and
    /// [`CorePhmm`].
    #[inline]
    #[must_use]
    pub fn from_parts(mapping: &'a ByteIndexMap<S>, core: CorePhmm<T, S>) -> Self {
        Self { mapping, core }
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
/// ## Parameters
///
/// - `'a`: The lifetime of the alphabet. This can be `'static` for compile-time
///   alphabets.
/// - `T`: The numeric type used to store the parameters.
/// - `S`: The size of the alphabet.
///
/// [`into_local_phmm`]: GlobalPhmm::into_local_phmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct LocalPhmm<'a, T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub(crate) mapping: &'a ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub(crate) core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    pub(crate) begin:   LocalModule<T, S>,
    /// The module for handling any bases after the core model.
    pub(crate) end:     LocalModule<T, S>,
}

impl<'a, T, const S: usize> LocalPhmm<'a, T, S> {
    /// Creates a new [`LocalPhmm`] from the specified parts.
    ///
    /// ## Errors
    ///
    /// [`IncompatibleModuleError`] is returned if the length of `begin` or
    /// `end` doesn't match the length of `core`.
    #[inline]
    pub fn from_parts(
        mapping: &'a ByteIndexMap<S>, core: CorePhmm<T, S>, begin: LocalModule<T, S>, end: LocalModule<T, S>,
    ) -> Result<Self, IncompatibleModuleError> {
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
}

/// An implementation of a profile hidden Markov model (pHMM) for domain
/// alignment (aligning a subsequence to a full model).
///
/// This is created from a [`GlobalPhmm`] using [`into_domain_phmm`]. Two
/// [`DomainModule`] modules are added to either end which can match arbitrarily
/// many bases at the beginning or end of the sequence.
///
/// ## Parameters
///
/// - `'a`: The lifetime of the alphabet. This can be `'static` for compile-time
///   alphabets.
/// - `T`: The numeric type used to store the parameters.
/// - `S`: The size of the alphabet.
///
/// [`into_domain_phmm`]: GlobalPhmm::into_domain_phmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct DomainPhmm<'a, T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub(crate) mapping: &'a ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub(crate) core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    pub(crate) begin:   DomainModule<T, S>,
    /// The module for handling any bases after the core model.
    pub(crate) end:     DomainModule<T, S>,
}

impl<'a, T, const S: usize> DomainPhmm<'a, T, S> {
    /// Creates a new [`DomainPhmm`] from the specified parts.
    #[inline]
    #[must_use]
    pub fn from_parts(
        mapping: &'a ByteIndexMap<S>, core: CorePhmm<T, S>, begin: DomainModule<T, S>, end: DomainModule<T, S>,
    ) -> Self {
        Self {
            mapping,
            core,
            begin,
            end,
        }
    }
}

/// An implementation of a profile hidden Markov model (pHMM) for semilocal
/// alignment (aligning a full sequence to a submodel).
///
/// This is created from a [`GlobalPhmm`] using [`into_semilocal_phmm`]. Two
/// [`SemiLocalModule`] modules are added to either end which can skip
/// arbitrarily many states at the beginning or end of the pHMM.
///
/// ## Parameters
///
/// - `'a`: The lifetime of the alphabet. This can be `'static` for compile-time
///   alphabets.
/// - `T`: The numeric type used to store the parameters.
/// - `S`: The size of the alphabet.
///
/// [`into_semilocal_phmm`]: GlobalPhmm::into_semilocal_phmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SemiLocalPhmm<'a, T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub(crate) mapping: &'a ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub(crate) core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model.
    pub(crate) begin:   SemiLocalModule<T>,
    /// The module for handling any bases after the core model.
    pub(crate) end:     SemiLocalModule<T>,
}

impl<'a, T, const S: usize> SemiLocalPhmm<'a, T, S> {
    /// Creates a new [`SemiLocalPhmm`] from the specified parts.
    ///
    /// ## Errors
    ///
    /// [`IncompatibleModuleError`] is returned if the length of `begin` or
    /// `end` doesn't match the length of `core`.
    #[inline]
    pub fn from_parts(
        mapping: &'a ByteIndexMap<S>, core: CorePhmm<T, S>, begin: SemiLocalModule<T>, end: SemiLocalModule<T>,
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

impl<'a, T: PhmmNumber, const S: usize> GlobalPhmm<'a, T, S> {
    /// Converts a [`GlobalPhmm`] into a [`LocalPhmm`] using the provided
    /// `config`.
    ///
    /// ## Errors
    ///
    /// [`IncompatibleModuleError`] is returned if the length of either module
    /// is incorrect (for [`LocalConfig::Custom`]).
    #[inline]
    pub fn into_local_phmm(self, config: LocalConfig<T, S>) -> Result<LocalPhmm<'a, T, S>, IncompatibleModuleError> {
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
    pub fn into_domain_phmm(self, config: DomainConfig<T, S>) -> DomainPhmm<'a, T, S> {
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
    pub fn into_semilocal_phmm(
        self, config: SemiLocalConfig<T>,
    ) -> Result<SemiLocalPhmm<'a, T, S>, IncompatibleModuleError> {
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

impl<T, const S: usize> GetModule for LocalPhmm<'_, T, S> {
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

impl<T, const S: usize> GetModule for DomainPhmm<'_, T, S> {
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

impl<T, const S: usize> GetModule for SemiLocalPhmm<'_, T, S> {
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

impl<T, const S: usize> GetCore<T, S> for GlobalPhmm<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for GlobalPhmm<'_, T, S> {
    #[inline]
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        self.core().layers()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for GlobalPhmm<'_, T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>> {
        self.core.layers_mut()
    }
}

impl<T, const S: usize> GetCore<T, S> for DomainPhmm<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for DomainPhmm<'_, T, S> {
    #[inline]
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        self.core().layers()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for DomainPhmm<'_, T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>> {
        self.core.layers_mut()
    }
}

impl<T, const S: usize> GetCore<T, S> for SemiLocalPhmm<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for SemiLocalPhmm<'_, T, S> {
    #[inline]
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        self.core().layers()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for SemiLocalPhmm<'_, T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>> {
        self.core.layers_mut()
    }
}

impl<T, const S: usize> GetCore<T, S> for LocalPhmm<'_, T, S> {
    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.core
    }
}

impl<T, const S: usize> GetLayer<T, S> for LocalPhmm<'_, T, S> {
    #[inline]
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        self.core().layers()
    }
}

impl<T, const S: usize> GetLayerMut<T, S> for LocalPhmm<'_, T, S> {
    #[inline]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>> {
        self.core.layers_mut()
    }
}

impl<'a, T, const S: usize> GetMapping<'a, S> for GlobalPhmm<'a, T, S> {
    #[inline]
    fn mapping(&self) -> &'a ByteIndexMap<S> {
        self.mapping
    }
}

impl<'a, T, const S: usize> GetMapping<'a, S> for LocalPhmm<'a, T, S> {
    #[inline]
    fn mapping(&self) -> &'a ByteIndexMap<S> {
        self.mapping
    }
}

impl<'a, T, const S: usize> GetMapping<'a, S> for SemiLocalPhmm<'a, T, S> {
    #[inline]
    fn mapping(&self) -> &'a ByteIndexMap<S> {
        self.mapping
    }
}

impl<'a, T, const S: usize> GetMapping<'a, S> for DomainPhmm<'a, T, S> {
    #[inline]
    fn mapping(&self) -> &'a ByteIndexMap<S> {
        self.mapping
    }
}

impl<'a, P, const S: usize> GetMapping<'a, S> for &P
where
    P: GetMapping<'a, S>,
{
    #[inline]
    fn mapping(&self) -> &'a ByteIndexMap<S> {
        (*self).mapping()
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
