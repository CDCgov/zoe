#![allow(dead_code)]
// TODO: Remove this when more pHMM stuff is added

use crate::{
    alignment::phmm::indexing::{LastMatch, PhmmIndex, PhmmIndexable},
    data::ByteIndexMap,
    math::Float,
};
use std::ops::{Index, IndexMut};

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

impl<T: Float> Default for TransitionParams<T> {
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

impl<T: Float, const S: usize> Default for EmissionParams<T, S> {
    /// Initialize the emission parameters so that all residues are probability
    /// zero
    #[inline]
    fn default() -> Self {
        // Zero probability is infinity in negative log space
        Self([T::INFINITY; S])
    }
}

// `f32` and `f64` both implement From<u16>
impl<T: Float + From<u16>, const S: usize> EmissionParams<T, S> {
    /// Generates the emission parameters for a uniform distribution.
    #[inline]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn uniform() -> Self {
        const { assert!(S < u16::MAX as usize) }
        Self([-(T::ONE / T::from(S as u16)).ln(); S])
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

impl<T: Float, const S: usize> Default for LayerParams<T, S> {
    #[inline]
    fn default() -> Self {
        Self {
            transition:      TransitionParams::default(),
            emission_match:  EmissionParams::default(),
            emission_insert: EmissionParams::default(),
        }
    }
}

// TODO: Fix doc comments (add local)
/// The core profile hidden Markov model (pHMM) used by [`GlobalPhmm`],
/// [`LocalPhmm`] and other PHMMs.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct CorePhmm<T, const S: usize>(
    /// The parameters for the pHMM.
    ///
    /// All probabilities are converted to log space with
    /// $-\operatorname{ln}(\cdot)$.
    ///
    /// Each element stores:
    /// * The transition probabilities into the insert state within the current
    ///   layer, as well as the insert state's emission probabilities
    /// * The transition probabilities from the current layer's states into the
    ///   next layer's match and delete states, as well as the next layer's
    ///   emission probabilities for the match state
    ///
    /// Note that this is different from SAM. The first element's match state is
    /// the BEGIN state, while the last element's transition probabilities
    /// reflect the probabilities entering the END state.
    pub(crate) Vec<LayerParams<T, S>>,
);

impl<T, const S: usize> CorePhmm<T, S> {
    /// Get a layer from within the core pHMM. Although there is no actual layer
    /// for the `End` state, for readability we let `End` be synonymous with
    /// `LastMatch` since `End` emphasizes it is the last layer.
    pub(crate) fn get_layer(&self, j: impl PhmmIndex) -> &LayerParams<T, S> {
        if j.is_end() {
            self.get_layer(LastMatch)
        } else {
            &self.0[self.get_dp_index(j)]
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
/// [`LocalModule`] modules are added to either end of the [`GlobalPhmm`], which
/// can match arbitrarily many bases at the beginning or end of the sequence,
/// and can skip arbitrarily many states at the beginning or end of the pHMM.
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

impl<T: Float, const S: usize> GlobalPhmm<T, S> {
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
}
