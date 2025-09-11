//! Structs and algorithms for profile Hidden Markov Models (pHMMs).
//!
//! <div class="warning note">
//!
//! **Note**
//!
//! You must enable the *dev-phmm* feature in your `Cargo.toml` to use these
//! functions. They are in active development and are not complete.
//!
//! </div>

#![allow(dead_code)]
// TODO: Remove this when more pHMM stuff is added

use crate::{
    alignment::phmm::indexing::{LastMatch, PhmmIndex},
    data::ByteIndexMap,
    math::{CastAs, CastAsNumeric, CastFrom, CastFromNumeric, Float},
};
use std::ops::{Add, AddAssign, Index, IndexMut, Mul};

mod alignment_modes;
mod errors;
mod indexing;
mod sam_parser;
mod score_from_path;
mod state;
mod traits;
mod viterbi;

pub use alignment_modes::*;
pub use errors::*;
pub(crate) use indexing::*;
pub use sam_parser::*;
pub(crate) use state::*;
pub use viterbi::*;

/// A trait for numeric types compatible with pHMMs.
///
/// These numeric types are used for performing pHMM calculations in negative
/// log space.
pub trait PhmmNumber:
    Copy + Add<Output = Self> + Mul<Output = Self> + AddAssign + PartialOrd + CastAs + CastAsNumeric + CastFrom + CastFromNumeric
{
    /// Infinity, the negative log space score corresponding to probability zero
    const INFINITY: Self;
    /// Zero, the negative log space score corresponding to probability one
    const ZERO: Self;

    /// Converts a probability (a floating point value) into negative log space
    fn from_prob<T: Float>(prob: T) -> Self;

    /// Converts a negative log space score back to a probability
    fn to_prob<T: Float>(self) -> T;

    /// Returns the negative log space score as a float
    fn to_float<T: Float>(self) -> T;

    /// Computes the minimum of two negative log space scores
    #[must_use]
    fn min(self, other: Self) -> Self;
}

impl PhmmNumber for f32 {
    const INFINITY: Self = f32::INFINITY;
    const ZERO: Self = 0.0;

    #[inline]
    fn from_prob<T: Float>(prob: T) -> Self {
        // Increase precision by converting to f32 last
        let param = (-prob.ln()).cast_as::<f32>();
        if param.is_nan() { Self::INFINITY } else { param }
    }

    #[inline]
    fn to_prob<T: Float>(self) -> T {
        // Increase precision by converting from f32 first
        (-T::cast_from(self)).exp()
    }

    #[inline]
    fn to_float<T: Float>(self) -> T {
        T::cast_from(self)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        self.min(other)
    }
}

impl PhmmNumber for f64 {
    const INFINITY: Self = f64::INFINITY;
    const ZERO: Self = 0.0;

    #[inline]
    fn from_prob<T: Float>(prob: T) -> Self {
        // Increase precision by converting to f64 first
        let param = -prob.cast_as::<f64>().ln();
        if param.is_nan() { Self::INFINITY } else { param }
    }

    #[inline]
    fn to_prob<T: Float>(self) -> T {
        // Increase precision by converting from f64 last
        T::cast_from((-self).exp())
    }

    #[inline]
    fn to_float<T: Float>(self) -> T {
        T::cast_from(self)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        self.min(other)
    }
}

/// The transition probabilities for a layer of the pHMM.
///
/// The parameters are converted to log space with $-\operatorname{ln}(\cdot)$.
/// They are arranged in the same format as SAM: `[[dd, md, id], [dm, mm, im],
/// [di, mi, ii]]`. This means the first element is the transition probabilities
/// into the delete state, the second element is the probabilities into the
/// match state, and the third element is the probabilities into the insert
/// state.
///
/// All parameters reflect the probability of transitioning from the
/// previous layer into the current layer, except for transitions into the
/// insert state, which are transitions within the same layer.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct TransitionParams<T>(pub(crate) [[T; 3]; 3]);

impl<T: PhmmNumber> Default for TransitionParams<T> {
    /// Initialize the transition parameters so that all transitions are
    /// probability zero
    #[inline]
    fn default() -> Self {
        // Zero probability is infinity in negative log space
        Self([[T::INFINITY; 3]; 3])
    }
}

impl<T> Index<(PhmmState, PhmmState)> for TransitionParams<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: (PhmmState, PhmmState)) -> &Self::Output {
        &self.0[usize::from(index.1)][usize::from(index.0)]
    }
}

impl<T> IndexMut<(PhmmState, PhmmState)> for TransitionParams<T> {
    #[inline]
    fn index_mut(&mut self, index: (PhmmState, PhmmState)) -> &mut Self::Output {
        &mut self.0[usize::from(index.1)][usize::from(index.0)]
    }
}

impl<T> Index<PhmmState> for TransitionParams<T> {
    type Output = [T; 3];

    #[inline]
    fn index(&self, index: PhmmState) -> &Self::Output {
        &self.0[usize::from(index)]
    }
}

impl<T> IndexMut<PhmmState> for TransitionParams<T> {
    #[inline]
    fn index_mut(&mut self, index: PhmmState) -> &mut Self::Output {
        &mut self.0[usize::from(index)]
    }
}

/// A set of emission probabilities, converted to log space with
/// $-\operatorname{ln}(\cdot)$.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct EmissionParams<T, const S: usize>(pub(crate) [T; S]);

impl<T: PhmmNumber, const S: usize> Default for EmissionParams<T, S> {
    /// Initialize the emission parameters so that all residues are probability
    /// zero
    #[inline]
    fn default() -> Self {
        // Zero probability is infinity in negative log space
        Self([T::INFINITY; S])
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

impl<T, const S: usize> Index<usize> for EmissionParams<T, S> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

/// The parameters for a single layer of the pHMM.
///
/// See [`TransitionParams`] and [`EmissionParams`] for more details.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct LayerParams<T, const S: usize> {
    pub(crate) transition:      TransitionParams<T>,
    pub(crate) emission_match:  EmissionParams<T, S>,
    pub(crate) emission_insert: EmissionParams<T, S>,
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
/// [`LocalPhmm`], [`DomainPhmm`], and other PHMMs.
///
/// This includes the layers of the pHMM without any modules at the beginning or
/// end.
///
/// All probabilities are converted to log space with
/// $-\operatorname{ln}(\cdot)$. This struct guarantees that at least two layers
/// are present (corresponding to a reference length of one).
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct CorePhmm<T, const S: usize>(
    /// Each element stores:
    /// * The transition probabilities into the insert state within the current
    ///   layer, as well as the insert state's emission probabilities
    /// * The transition probabilities from the current layer's states into the next
    ///   layer's match and delete states, as well as the next layer's emission
    ///   probabilities for the match state
    ///
    /// Note that this is different from SAM. The first element's match state is the
    /// BEGIN state, while the last element's transition probabilities reflect the
    /// probabilities entering the END state.
    Vec<LayerParams<T, S>>,
);

impl<T, const S: usize> CorePhmm<T, S> {
    /// Create a new [`CorePhmm`] from a `Vec` of the parameters.
    ///
    /// ## Errors
    ///
    /// At least two layers are required, otherwise [`PhmmError::TooFewLayers`]
    /// is returned.
    pub(crate) fn new(layers: Vec<LayerParams<T, S>>) -> Result<Self, PhmmError> {
        if layers.len() >= 2 {
            Ok(CorePhmm(layers))
        } else {
            Err(PhmmError::TooFewLayers(2))
        }
    }

    /// Creates a new [`CorePhmm`] from a `Vec` of the parameters, without
    /// checking that the number of layers is valid (at least 2).
    pub(crate) fn new_unchecked(layers: Vec<LayerParams<T, S>>) -> Self {
        CorePhmm(layers)
    }

    /// Get a layer from within the core pHMM.
    ///
    /// Although there is no actual layer for the `End` state, for readability
    /// we let `End` be synonymous with `LastMatch` since `End` emphasizes it is
    /// the last layer.
    pub(crate) fn get_layer(&self, j: impl PhmmIndex) -> &LayerParams<T, S> {
        if j.is_end() {
            self.get_layer(LastMatch)
        } else {
            &self.0[self.get_dp_index(j)]
        }
    }

    /// Get a mutable reference to a layer from within the core pHMM.
    ///
    /// Although there is no actual layer for the `End` state, for readability
    /// we let `End` be synonymous with `LastMatch` since `End` emphasizes it is
    /// the last layer.
    pub(crate) fn get_layer_mut(&mut self, j: impl PhmmIndex) -> &mut LayerParams<T, S> {
        if j.is_end() {
            self.get_layer_mut(LastMatch)
        } else {
            let idx = self.get_dp_index(j);
            &mut self.0[idx]
        }
    }

    /// Returns the length of the "reference" represented by the pHMM
    #[inline]
    #[must_use]
    pub fn ref_length(&self) -> usize {
        // A pHMM of length n has n transitions between layers, which means n+1
        // layers. One is the begin, and one is the end. So n-1 layers
        // correspond to indices in the reference.
        self.0.len() - 1
    }
}

/// An implementation of a profile hidden Markov model (pHMM) for global
/// alignment (aligning a full sequence to a full model).
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct GlobalPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub mapping: &'static ByteIndexMap<S>,
    /// The model parameters
    pub core:    CorePhmm<T, S>,
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
#[derive(Debug)]
pub struct LocalPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model
    pub begin:   LocalModule<T, S>,
    /// The module for handling any bases after the core model
    pub end:     LocalModule<T, S>,
}

impl<T: PhmmNumber, const S: usize> LocalPhmm<T, S> {
    /// Gets the score incurred by skipping `inserted` bases from the beginning
    /// of the query, before entering the [`CorePhmm`].
    pub(crate) fn get_begin_internal_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        self.begin.internal_params.get_begin_score(inserted, mapping)
    }

    /// Gets the score incurred by skipping to `index` while entering the
    /// [`CorePhmm`].
    pub(crate) fn get_begin_external_score(&self, index: impl PhmmIndex) -> T {
        self.begin.external_params.get_score(index)
    }

    /// Gets the score incurred by exiting the [`CorePhmm`] from `index`.
    pub(crate) fn get_end_external_score(&self, index: impl PhmmIndex) -> T {
        self.end.external_params.get_score(index)
    }

    /// Gets the score incurred by skipping `inserted` bases from the end of the
    /// query, after exiting the [`CorePhmm`].
    pub(crate) fn get_end_internal_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        self.end.internal_params.get_end_score(inserted, mapping)
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
#[derive(Debug)]
pub struct DomainPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model
    pub begin:   DomainModule<T, S>,
    /// The module for handling any bases after the core model
    pub end:     DomainModule<T, S>,
}

impl<T: PhmmNumber, const S: usize> DomainPhmm<T, S> {
    /// Gets the score for transitioning into a given [`PhmmIndex`] from the
    /// [`DomainModule`] at the beginning of the pHMM.
    ///
    /// This is lazily computed. Use [`PrecomputedDomainModule`] for a more
    /// efficienct alternative.
    #[inline]
    pub(crate) fn get_begin_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        self.begin.get_begin_score(inserted, mapping)
    }

    /// Gets the score for transitioning out of a given [`PhmmIndex`] into the
    /// [`DomainModule`] at the end of the pHMM.
    ///
    /// This is lazily computed. Use [`PrecomputedDomainModule`] for a more
    /// efficienct alternative.
    #[inline]
    pub(crate) fn get_end_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        self.end.get_end_score(inserted, mapping)
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
#[derive(Debug)]
pub struct SemiLocalPhmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    pub mapping: &'static ByteIndexMap<S>,
    /// The core model containing the parameters.
    pub core:    CorePhmm<T, S>,
    /// The module for handling any bases before the core model
    pub begin:   SemiLocalModule<T>,
    /// The module for handling any bases after the core model
    pub end:     SemiLocalModule<T>,
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<T, S> {
    /// Gets the score for transitioning into a given [`PhmmIndex`] from the
    /// [`SemiLocalModule`] at the beginning of the pHMM.
    pub(crate) fn get_begin_score(&self, index: impl PhmmIndex) -> T {
        self.begin.get_score(index)
    }

    /// Gets the score for transitioning out of a given [`PhmmIndex`] into the
    /// [`SemiLocalModule`] at the end of the pHMM.
    pub(crate) fn get_end_score(&self, index: impl PhmmIndex) -> T {
        self.end.get_score(index)
    }
}

/// Options for how to construct a [`LocalPhmm`] from a [`GlobalPhmm`]
#[non_exhaustive]
pub enum LocalConfig<T, const S: usize> {
    /// Does not penalize any transitions in the [`LocalModule`], instead solely
    /// penalizing the emissions for any insertions
    NoPenalty { background_emission: EmissionParams<T, S> },
    /// Use a custom [`LocalModule`] for the begin and end.
    Custom {
        begin: LocalModule<T, S>,
        end:   LocalModule<T, S>,
    },
}

/// Options for how to construct a [`DomainPhmm`] from a [`GlobalPhmm`]
#[non_exhaustive]
pub enum DomainConfig<T, const S: usize> {
    /// Does not penalize any transitions in the [`DomainModule`], instead
    /// solely penalizing the emissions for any insertions
    NoPenalty { background_emission: EmissionParams<T, S> },
    /// Use a custom [`DomainModule`] for the begin and end.
    Custom {
        begin: DomainModule<T, S>,
        end:   DomainModule<T, S>,
    },
}

/// Options for how to construct a [`SemiLocalPhmm`] from a [`GlobalPhmm`]
#[non_exhaustive]
pub enum SemiLocalConfig<T> {
    /// Does not penalize any transitions in the [`SemiLocalModule`]
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
