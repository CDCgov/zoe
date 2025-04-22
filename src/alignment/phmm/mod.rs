#![allow(dead_code)]
// TODO: Remove this when more pHMM stuff is added

use crate::{data::ByteIndexMap, math::Float};
use std::ops::{Index, IndexMut};

mod sam_parser;

pub use sam_parser::*;

/// An enum representing the three states within each layer of a pHMM. This is
/// used for readability when indexing.
#[derive(Clone, Copy, Debug)]
pub enum PhmmState {
    Delete,
    Match,
    Insert,
}

impl From<PhmmState> for usize {
    #[inline]
    fn from(value: PhmmState) -> Self {
        match value {
            PhmmState::Delete => 0,
            PhmmState::Match => 1,
            PhmmState::Insert => 2,
        }
    }
}

/// Stores three values, one associated to each pHMM state (delete, match, and
/// insert). This is used for readability, allowing indexing with `PhmmState`.
#[derive(Clone, Debug)]
pub struct InfoByPhmmState<T>([T; 3]);

impl<T> Index<PhmmState> for InfoByPhmmState<T> {
    type Output = T;

    fn index(&self, index: PhmmState) -> &Self::Output {
        &self.0[usize::from(index)]
    }
}

impl<T> IndexMut<PhmmState> for InfoByPhmmState<T> {
    fn index_mut(&mut self, index: PhmmState) -> &mut Self::Output {
        &mut self.0[usize::from(index)]
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
pub struct TransitionParams<T>([[T; 3]; 3]);

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

/// A set of emission probabilities, converted to log space with
/// $-\operatorname{ln}(\cdot)$.
pub struct EmissionParams<T, const S: usize>([T; S]);

impl<T: Float + From<u16>, const S: usize> EmissionParams<T, S> {
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
pub struct LayerParams<T, const S: usize> {
    transition:      TransitionParams<T>,
    emission_match:  EmissionParams<T, S>,
    emission_insert: EmissionParams<T, S>,
}

/// An implementation of a profile hidden Markov model (pHMM).
pub struct Phmm<T, const S: usize> {
    /// The mapping used when processing the bases. This will vary depending on
    /// the alphabet used.
    mapping: &'static ByteIndexMap<S>,

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
    params: Vec<LayerParams<T, S>>,
}
