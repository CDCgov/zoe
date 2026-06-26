use crate::alignment::phmm::{
    InvalidModelError, PhmmError, PhmmNumber,
    indexing::{GetLayer, GetLayerMut},
    state::PhmmState,
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
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct TransitionParams<T>(pub(crate) [[T; 3]; 3]);

impl<T: Copy> TransitionParams<T> {
    /// Retrieves an array containing the three parameters for exiting `state`
    /// and moving to the next layer.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn entering_params(&self, state: PhmmState) -> [T; 3] {
        self[state]
    }

    /// Retrieves the inner array from the [`TransitionParams`]. The layout of
    /// this array is subject to change. Consider indexing directly into the
    /// [`TransitionParams`] instead.
    #[inline]
    #[must_use]
    #[cfg(feature = "dev-phmm-regression")]
    pub fn as_array(&self) -> &[[T; 3]; 3] {
        &self.0
    }

    /// Constructs a [`TransitionParams`] from an array. The interpretation of
    /// this array is subject to change. Consider using [`Default`] and then
    /// mutating each entry.
    #[inline]
    #[must_use]
    #[cfg(feature = "dev-phmm-regression")]
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
    ///
    /// [`ByteIndexMap`]: crate::data::mappings::ByteIndexMap
    #[inline]
    #[must_use]
    pub fn from_array(arr: [T; S]) -> Self {
        Self(arr)
    }

    /// Returns the parameters as a slice.
    ///
    /// The order of the elements matches the order of the keys in the
    /// [`ByteIndexMap`] used by the pHMM.
    ///
    /// [`ByteIndexMap`]: crate::data::mappings::ByteIndexMap
    #[inline]
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.0
    }

    /// Returns a reference to the underlying array storing the parameters.
    ///
    /// The order of the elements matches the order of the keys in the
    /// [`ByteIndexMap`] used by the pHMM.
    ///
    /// [`ByteIndexMap`]: crate::data::mappings::ByteIndexMap
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
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct LayerParams<T, const S: usize> {
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
/// [`GlobalPhmm`]: crate::alignment::phmm::models::GlobalPhmm
/// [`DomainPhmm`]: crate::alignment::phmm::models::DomainPhmm
/// [`LocalPhmm`]: crate::alignment::phmm::models::LocalPhmm
/// [`SemiLocalPhmm`]: crate::alignment::phmm::models::SemiLocalPhmm
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct CorePhmm<T, const S: usize>(Vec<LayerParams<T, S>>);

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

impl<T, const S: usize> GetLayer<T, S> for CorePhmm<T, S> {
    #[inline]
    fn layers(&self) -> &[LayerParams<T, S>] {
        self.0.as_slice()
    }

    #[inline]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]) {
        self.0.as_slice().split_first().expect("A CorePhmm has at least two layers")
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
