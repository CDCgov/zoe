use crate::alignment::phmm::{
    CorePhmm, LayerParams,
    indexing::{PhmmIndex, PhmmIndexRange, PhmmIndexable},
};
use std::slice::GetDisjointMutError;

/// A trait providing read-only access to the modules at the beginning and end
/// of a pHMM.
pub trait GetModule {
    /// The type of the module at the beginning of the pHMM.
    type Begin;
    /// The type of the module at the end of the pHMM.
    type End;

    /// Returns a reference to the module at the start of the pHMM.
    #[must_use]
    fn begin(&self) -> &Self::Begin;

    /// Returns a reference to the module at the end of the pHMM.
    #[must_use]
    fn end(&self) -> &Self::End;
}

/// A trait providing mutable access to the modules at the beginning and end of
/// a pHMM.
#[allow(dead_code)]
pub trait GetModuleMut: GetModule {
    /// Returns a mutable reference to the module at the start of the pHMM.
    #[must_use]
    fn begin_mut(&mut self) -> &mut Self::Begin;

    /// Returns a mutable reference to the module at the end of the pHMM.
    #[must_use]
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
#[allow(dead_code)]
pub trait GetLayer<T, const S: usize>: PhmmIndexable {
    /// Retrieves a slice of the layers contained within the core pHMM.
    ///
    /// This slice will be at least 2 in length.
    #[must_use]
    fn layers(&self) -> &[LayerParams<T, S>];

    /// Returns the first layer, as well as all subsequent layers.
    ///
    /// This is an infallible version of `model.layers().split_first()`.
    #[must_use]
    fn split_first_layer(&self) -> (&LayerParams<T, S>, &[LayerParams<T, S>]);

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
    fn get_layers(&self, range: impl PhmmIndexRange) -> Option<&[LayerParams<T, S>]> {
        self.layers().get(self.get_dp_range(range))
    }
}

/// A trait providing mutable accessors to the layers of a pHMM.
#[allow(dead_code)]
pub trait GetLayerMut<T, const S: usize>: GetLayer<T, S> {
    /// Retrieves a mutable slice of the layers contained within the core pHMM.
    #[must_use]
    fn layers_mut(&mut self) -> &mut [LayerParams<T, S>];

    /// Returns a mutable reference to the vector of layer parameters stored in
    /// the core pHMM.
    #[must_use]
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
pub trait GetPartsMut<T, const S: usize>: GetModule + GetLayer<T, S> {
    /// Returns simultaneous mutable references to the core pHMM, the begin
    /// module, and the end module.
    ///
    /// Calling the individual mutable accessors and simultaneously using them
    /// is not allowed by the borrow checker, hence this function.
    #[must_use]
    fn parts_mut(&mut self) -> (&mut CorePhmm<T, S>, &mut Self::Begin, &mut Self::End);
}
