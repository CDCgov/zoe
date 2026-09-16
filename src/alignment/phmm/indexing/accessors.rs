use crate::{
    alignment::phmm::{
        at_least_two::VecAtLeast2,
        components::{CorePhmm, LayerParams},
        indexing::{AlnIndex, AlnIndexRange, AlnIndexable, Begin, FirstResidue, IndexRangeInner, LastResidue},
    },
    data::ByteIndexMap,
};

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

/// A trait providing read-only accessors to the layers of a pHMM.
pub trait GetLayer<T, const S: usize>: AlnIndexable {
    /// Retrieves a vector of the layers contained within the core pHMM.
    ///
    /// This vector will be at least 2 in length.
    #[must_use]
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>>;

    /// Returns the layers, split at a particular index.
    ///
    /// If the index is out of bounds, `None` is returned.
    #[inline]
    #[must_use]
    #[allow(clippy::type_complexity)]
    fn split_layers_at(&self, phmm_idx: impl AlnIndex) -> Option<(&[LayerParams<T, S>], &[LayerParams<T, S>])> {
        self.layers().split_at_checked(phmm_idx.to_dp_index(self).0)
    }

    /// Gets a layer from within the core pHMM.
    ///
    /// This returns `None` if the index is out of bounds or [`End`] (since
    /// there is no layer corresponding to the END state).
    ///
    /// [`End`]: crate::alignment::phmm::indexing::End
    #[inline]
    #[must_use]
    fn get_layer(&self, phmm_idx: impl AlnIndex) -> Option<&LayerParams<T, S>> {
        self.layers().get(phmm_idx.to_dp_index(self).0)
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
    fn get_layers(&self, range: impl AlnIndexRange) -> Option<&[LayerParams<T, S>]> {
        let range = range.to_dp_range(self);
        self.layers().get(range.into_inner())
    }
}

impl<P: GetLayer<T, S>, T, const S: usize> GetLayer<T, S> for &P {
    fn layers(&self) -> &VecAtLeast2<LayerParams<T, S>> {
        P::layers(self)
    }
}

/// A trait providing mutable accessors to the layers of a pHMM.
pub(crate) trait GetLayerMut<T, const S: usize>: GetLayer<T, S> {
    /// Retrieves a mutable vector of the layers contained within the core pHMM.
    #[must_use]
    fn layers_mut(&mut self) -> &mut VecAtLeast2<LayerParams<T, S>>;

    /// Returns a mutable reference to the parameters for the specified layer
    /// which is guaranteed to exist in the pHMM, either [`Begin`],
    /// [`FirstResidue`], or [`LastResidue`].
    #[inline]
    #[must_use]
    fn layer_mut(&mut self, idx: impl InfallibleLayerIdx) -> &mut LayerParams<T, S> {
        idx.layer_mut(self.layers_mut())
    }
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
    fn layer<T, const S: usize>(self, layers: &VecAtLeast2<LayerParams<T, S>>) -> &LayerParams<T, S>;

    /// Returns a mutable reference to the specified layer from `layers`.
    fn layer_mut<T, const S: usize>(self, layers: &mut VecAtLeast2<LayerParams<T, S>>) -> &mut LayerParams<T, S>;
}

impl InfallibleLayerIdx for Begin {
    #[inline]
    fn layer<T, const S: usize>(self, layers: &VecAtLeast2<LayerParams<T, S>>) -> &LayerParams<T, S> {
        layers.first()
    }

    #[inline]
    fn layer_mut<T, const S: usize>(self, layers: &mut VecAtLeast2<LayerParams<T, S>>) -> &mut LayerParams<T, S> {
        layers.first_mut()
    }
}

impl InfallibleLayerIdx for FirstResidue {
    #[inline]
    fn layer<T, const S: usize>(self, layers: &VecAtLeast2<LayerParams<T, S>>) -> &LayerParams<T, S> {
        layers.second()
    }

    #[inline]
    fn layer_mut<T, const S: usize>(self, layers: &mut VecAtLeast2<LayerParams<T, S>>) -> &mut LayerParams<T, S> {
        layers.second_mut()
    }
}

impl InfallibleLayerIdx for LastResidue {
    #[inline]
    fn layer<T, const S: usize>(self, layers: &VecAtLeast2<LayerParams<T, S>>) -> &LayerParams<T, S> {
        layers.last()
    }

    #[inline]
    fn layer_mut<T, const S: usize>(self, layers: &mut VecAtLeast2<LayerParams<T, S>>) -> &mut LayerParams<T, S> {
        layers.last_mut()
    }
}
