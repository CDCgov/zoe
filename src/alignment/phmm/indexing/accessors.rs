use crate::{
    alignment::phmm::{
        components::{CorePhmm, LayerParams},
        indexing::{Begin, FirstMatch, LastMatch, PhmmIndex, PhmmIndexRange, PhmmIndexable},
        nonempty_vec::NonEmptyVec,
    },
    data::ByteIndexMap,
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
pub trait GetCore<T, const S: usize> {
    /// Returns a reference to the [`CorePhmm`] holding the core parameters.
    #[must_use]
    fn core(&self) -> &CorePhmm<T, S>;
}

// This is a separate trait from GetLayerMut in order to prevent core_mut from
// being called on a CorePhmm, which is an easy way to have infinite recursion
// in an implementation.

/// A trait providing read-only access to the [`CorePhmm`] within a larger pHMM.
pub trait GetCoreMut<T, const S: usize> {
    /// Returns a mutable reference to the [`CorePhmm`] holding the core
    /// parameters.
    #[must_use]
    #[allow(dead_code)]
    fn core_mut(&mut self) -> &mut CorePhmm<T, S>;
}

/// A trait providing read-only accessors to the layers of a pHMM.
pub trait GetLayer<T, const S: usize>: PhmmIndexable {
    /// Retrieves a slice of the layers contained within the core pHMM.
    ///
    /// This slice will be at least 2 in length.
    #[must_use]
    fn layers(&self) -> &NonEmptyVec<LayerParams<T, S>>;

    /// Returns the layers, split at a particular index.
    ///
    /// If the index is out of bounds, `None` is returned.
    #[inline]
    #[must_use]
    #[allow(clippy::type_complexity)]
    fn split_layers_at(&self, j: impl PhmmIndex) -> Option<(&[LayerParams<T, S>], &[LayerParams<T, S>])> {
        self.layers().split_at_checked(self.get_dp_index(j))
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

    /// Returns a reference to the parameters for the specified layer which is
    /// guaranteed to exist in the pHMM, either [`Begin`], [`FirstResidue`], or
    /// [`LastResidue`].
    #[inline]
    #[must_use]
    fn layer(&self, idx: impl InfallibleLayerIdx) -> &LayerParams<T, S> {
        idx.layer(self.layers())
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

impl<P: GetLayer<T, S>, T, const S: usize> GetLayer<T, S> for &P {
    fn layers(&self) -> &NonEmptyVec<LayerParams<T, S>> {
        P::layers(self)
    }
}

/// A trait providing mutable accessors to the layers of a pHMM.
pub trait GetLayerMut<T, const S: usize>: GetLayer<T, S> {
    /// Retrieves a mutable slice of the layers contained within the core pHMM.
    #[must_use]
    fn layers_mut(&mut self) -> &mut NonEmptyVec<LayerParams<T, S>>;

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

    /// Returns a mutable reference to the parameters for the specified layer
    /// which is guaranteed to exist in the pHMM, either [`Begin`],
    /// [`FirstMatch`], or [`LastMatch`].
    #[inline]
    #[must_use]
    fn layer_mut(&mut self, idx: impl InfallibleLayerIdx) -> &mut LayerParams<T, S> {
        idx.layer_mut(self.layers_mut())
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

/// A trait providing access to the underlying alphabet of a pHMM.
pub trait GetMapping<const S: usize> {
    /// Returns a reference to the underlying alphabet of the pHMM.
    #[must_use]
    fn mapping(&self) -> &'static ByteIndexMap<S>;
}

/// A trait unifying [`Begin`], [`FirstResidue`], and [`LastResidue`], which are
/// layers of a pHMM that are guaranteed to exist.
pub trait InfallibleLayerIdx: Copy {
    /// Returns the specified layer from `layers`.
    fn layer<T, const S: usize>(self, layers: &NonEmptyVec<LayerParams<T, S>>) -> &LayerParams<T, S>;

    /// Returns a mutable reference to the specified layer from `layers`.
    fn layer_mut<T, const S: usize>(self, layers: &mut NonEmptyVec<LayerParams<T, S>>) -> &mut LayerParams<T, S>;
}

impl InfallibleLayerIdx for Begin {
    #[inline]
    fn layer<T, const S: usize>(self, layers: &NonEmptyVec<LayerParams<T, S>>) -> &LayerParams<T, S> {
        layers.first()
    }

    #[inline]
    fn layer_mut<T, const S: usize>(self, layers: &mut NonEmptyVec<LayerParams<T, S>>) -> &mut LayerParams<T, S> {
        layers.first_mut()
    }
}

impl InfallibleLayerIdx for FirstMatch {
    #[inline]
    fn layer<T, const S: usize>(self, layers: &NonEmptyVec<LayerParams<T, S>>) -> &LayerParams<T, S> {
        &layers[1]
    }

    #[inline]
    fn layer_mut<T, const S: usize>(self, layers: &mut NonEmptyVec<LayerParams<T, S>>) -> &mut LayerParams<T, S> {
        &mut layers[1]
    }
}

impl InfallibleLayerIdx for LastMatch {
    #[inline]
    fn layer<T, const S: usize>(self, layers: &NonEmptyVec<LayerParams<T, S>>) -> &LayerParams<T, S> {
        layers.last()
    }

    #[inline]
    fn layer_mut<T, const S: usize>(self, layers: &mut NonEmptyVec<LayerParams<T, S>>) -> &mut LayerParams<T, S> {
        layers.last_mut()
    }
}
