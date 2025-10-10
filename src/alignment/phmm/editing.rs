//! Functions for editing a pHMM after it has been created/loaded.

use crate::alignment::phmm::{
    CorePhmm, DomainPhmm, GetLayer, GetModule, GlobalPhmm, LocalPhmm, PhmmNumber, PhmmState, SemiLocalPhmm,
    indexing::{PhmmIndex, PhmmIndexRange, PhmmIndexable},
};

/// A trait facilitating editing of a pHMM by removing layers.
#[allow(dead_code)]
pub trait RemoveLayer {
    /// Removes a single layer from the pHMM along with the transitions from the
    /// previous layer.
    ///
    /// Consider three consecutive layers `0,1,2` of a pHMM. Removing layer 1
    /// with this function will erase the emission parameters for the states
    /// corresponding to layer 1, as well as the transition probabilities from 0
    /// to 1. The transition probabilities from 1 to 2 are now used for 0 to 2.
    ///
    /// This method may result in a model that is no longer strictly
    /// probabilistic (i.e., the probabilities out of a given state may not sum
    /// to 1).
    ///
    /// ## Panics
    ///
    /// The layer must not be the BEGIN layer (dynamic programming index 0)
    /// since there are no transitions into this layer.
    fn remove_layer_and_pre_transition(&mut self, layer: impl PhmmIndex);

    /// Removes a single layer from the pHMM along with the transitions into the
    /// next layer.
    ///
    /// Consider three consecutive layers `0,1,2` of a pHMM. Removing layer 1
    /// with this function will erase the emission parameters for the states
    /// corresponding to layer 1, as well as the transition probabilities from 1
    /// to 2. The transition probabilities from 0 to 1 are now used for 0 to 2.
    ///
    /// If the BEGIN layer is removed with this function, then the layer
    /// following it becomes the new BEGIN layer (and hence will not have a
    /// match state that consumes any residues).
    ///
    /// This method may result in a model that is no longer strictly
    /// probabilistic (i.e., the probabilities out of a given state may not sum
    /// to 1).
    ///
    /// ## Panics
    ///
    /// The layer must be in the range of the pHMM.
    fn remove_layer_and_post_transition(&mut self, layer: impl PhmmIndex);

    /// Removes a range of layers from the pHMM along with the transitions into
    /// those layer.
    ///
    /// This is equivalent to calling [`remove_layer_and_pre_transition`] on
    /// each index in the range. For example, consider five consecutive layers
    /// `0,1,2,3,4` of a pHMM. Removing layers `1..=3` with this function will
    /// erase the emission parameters for the states corresponding to layer
    /// those layers, as well as the transition probabilities from 0 to 1, 1 to
    /// 2, and 2 to 3. The transition probabilities from 3 to 4 are now used for
    /// 0 to 4.
    ///
    /// This method may result in a model that is no longer strictly
    /// probabilistic (i.e., the probabilities out of a given state may not sum
    /// to 1).
    ///
    /// ## Panics
    ///
    /// The range must not include the BEGIN layer (dynamic programming index 0)
    /// since there are no transitions into this layer.
    ///
    /// [`remove_layer_and_pre_transition`]:
    ///     RemoveLayer::remove_layer_and_pre_transition
    fn remove_layers_and_pre_transitions<R: PhmmIndexRange>(&mut self, range: R);

    /// Removes a range of layers from the pHMM along with the transitions out
    /// of those layers.
    ///
    /// This is equivalent to calling [`remove_layer_and_post_transition`] on
    /// each index in the range. For example, consider five consecutive layers
    /// `0,1,2,3,4` of a pHMM. Removing layers `1..=3` with this function will
    /// erase the emission parameters for the states corresponding to layer
    /// those layers, as well as the transition probabilities from 1 to 2, 2 to
    /// 3, and 3 to 4. The transition probabilities from 0 to 1 are now used for
    /// 0 to 4.
    ///
    /// If the BEGIN layer is included in the range, then the first layer after
    /// the range becomes the new BEGIN layer (and hence will not have a match
    /// state that consumes any residues).
    ///
    /// This method may result in a model that is no longer strictly
    /// probabilistic (i.e., the probabilities out of a given state may not sum
    /// to 1).
    ///
    /// ## Panics
    ///
    /// The layers must be in the range of the pHMM.
    ///
    /// [`remove_layer_and_post_transition`]:
    ///     RemoveLayer::remove_layer_and_post_transition
    fn remove_layers_and_post_transitions<R: PhmmIndexRange>(&mut self, range: R);
}

impl<T: PhmmNumber, const S: usize> RemoveLayer for CorePhmm<T, S> {
    fn remove_layer_and_pre_transition(&mut self, layer: impl PhmmIndex) {
        // Zoe stores the transitions from the current layer into the next, so
        // we must shift these transitions back before deleting so that the
        let prev_layer = layer.prev_index(self);
        let (layer1, layer2) = self.get_two_layers_mut(prev_layer, layer);
        std::mem::swap(
            &mut layer1.transition[PhmmState::Match],
            &mut layer2.transition[PhmmState::Match],
        );
        std::mem::swap(
            &mut layer1.transition[PhmmState::Delete],
            &mut layer2.transition[PhmmState::Delete],
        );
        std::mem::swap(&mut layer1.emission_match, &mut layer2.emission_match);
        let layer_idx = self.get_dp_index(layer);
        self.layers_mut_vec().remove(layer_idx);
    }

    fn remove_layer_and_post_transition(&mut self, layer: impl PhmmIndex) {
        if self.get_dp_index(layer) > 0 {
            // We are not deleting the begin state, so the previous layer holds
            // the match emission for this layer which must be deleted, and the
            // this layer holds the match emission for the next layer which must
            // be saved
            let prev_layer = layer.prev_index(self);
            let (layer1, layer2) = self.get_two_layers_mut(prev_layer, layer);
            std::mem::swap(&mut layer1.emission_match, &mut layer2.emission_match);
        }

        let layer_idx = self.get_dp_index(layer);
        self.layers_mut_vec().remove(layer_idx);
    }

    fn remove_layers_and_pre_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        let range = self.to_dp_range(range);
        if range.is_empty() {
            return;
        }

        let swap1_idx = range.start.prev_index(self);
        let swap2_idx = range.end.prev_index(self);
        let (layer1, layer2) = self.get_two_layers_mut(swap1_idx, swap2_idx);
        std::mem::swap(
            &mut layer1.transition[PhmmState::Match],
            &mut layer2.transition[PhmmState::Match],
        );
        std::mem::swap(
            &mut layer1.transition[PhmmState::Delete],
            &mut layer2.transition[PhmmState::Delete],
        );
        std::mem::swap(&mut layer1.emission_match, &mut layer2.emission_match);
        let layer_idxs = self.get_dp_range(range);
        self.layers_mut_vec().drain(layer_idxs);
    }

    fn remove_layers_and_post_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        // Ensure BEGIN is not modified by converting to sequence range
        let range = self.to_seq_range(range);
        if range.is_empty() {
            return;
        }

        if self.get_dp_index(range.start) > 0 {
            // We are not deleting the begin state, so the layer before
            // range.start holds the match emission for range.start which must
            // be deleted, and the layer before range.end holds the match
            // emission for range.end which is not deleted and hence must be
            // saved
            let swap1_idx = range.start.prev_index(self);
            let swap2_idx = range.end.prev_index(self);
            let (layer1, layer2) = self.get_two_layers_mut(swap1_idx, swap2_idx);
            std::mem::swap(&mut layer1.emission_match, &mut layer2.emission_match);
        }

        let layer_idxs = self.get_dp_range(range);
        self.layers_mut_vec().drain(layer_idxs);
    }
}

impl<T: PhmmNumber, const S: usize> RemoveLayer for GlobalPhmm<T, S> {
    #[inline]
    fn remove_layer_and_pre_transition(&mut self, layer: impl PhmmIndex) {
        self.core_mut().remove_layer_and_pre_transition(layer);
    }

    #[inline]
    fn remove_layer_and_post_transition(&mut self, layer: impl PhmmIndex) {
        self.core_mut().remove_layer_and_post_transition(layer);
    }

    #[inline]
    fn remove_layers_and_pre_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        self.core_mut().remove_layers_and_pre_transitions(range);
    }

    #[inline]
    fn remove_layers_and_post_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        self.core_mut().remove_layers_and_post_transitions(range);
    }
}

impl<T: PhmmNumber, const S: usize> RemoveLayer for LocalPhmm<T, S> {
    #[inline]
    fn remove_layer_and_pre_transition(&mut self, layer: impl PhmmIndex) {
        let layer_idx = self.get_dp_index(layer);
        self.begin_mut().external_params.0.remove(layer_idx);
        self.end_mut().external_params.0.remove(layer_idx);
        self.core_mut().remove_layer_and_pre_transition(layer);
    }

    #[inline]
    fn remove_layer_and_post_transition(&mut self, layer: impl PhmmIndex) {
        let layer_idx = self.get_dp_index(layer);
        self.begin_mut().external_params.0.remove(layer_idx);
        self.end_mut().external_params.0.remove(layer_idx);
        self.core_mut().remove_layer_and_post_transition(layer);
    }

    #[inline]
    fn remove_layers_and_pre_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        let range_idxs = self.get_dp_range(range.clone());
        self.begin_mut().external_params.0.drain(range_idxs.clone());
        self.end_mut().external_params.0.drain(range_idxs);
        self.core_mut().remove_layers_and_pre_transitions(range);
    }

    #[inline]
    fn remove_layers_and_post_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        let range_idxs = self.get_dp_range(range.clone());
        self.begin_mut().external_params.0.drain(range_idxs.clone());
        self.end_mut().external_params.0.drain(range_idxs);
        self.core_mut().remove_layers_and_post_transitions(range);
    }
}

impl<T: PhmmNumber, const S: usize> RemoveLayer for SemiLocalPhmm<T, S> {
    #[inline]
    fn remove_layer_and_pre_transition(&mut self, layer: impl PhmmIndex) {
        let layer_idx = self.get_dp_index(layer);
        self.begin_mut().0.remove(layer_idx);
        self.end_mut().0.remove(layer_idx);
        self.core_mut().remove_layer_and_pre_transition(layer);
    }

    #[inline]
    fn remove_layer_and_post_transition(&mut self, layer: impl PhmmIndex) {
        let layer_idx = self.get_dp_index(layer);
        self.begin_mut().0.remove(layer_idx);
        self.end_mut().0.remove(layer_idx);
        self.core_mut().remove_layer_and_post_transition(layer);
    }

    #[inline]
    fn remove_layers_and_pre_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        let range_idxs = self.get_dp_range(range.clone());
        self.begin_mut().0.drain(range_idxs.clone());
        self.end_mut().0.drain(range_idxs);
        self.core_mut().remove_layers_and_pre_transitions(range);
    }

    #[inline]
    fn remove_layers_and_post_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        let range_idxs = self.get_dp_range(range.clone());
        self.begin_mut().0.drain(range_idxs.clone());
        self.end_mut().0.drain(range_idxs);
        self.core_mut().remove_layers_and_post_transitions(range);
    }
}

impl<T: PhmmNumber, const S: usize> RemoveLayer for DomainPhmm<T, S> {
    #[inline]
    fn remove_layer_and_pre_transition(&mut self, layer: impl PhmmIndex) {
        self.core_mut().remove_layer_and_pre_transition(layer);
    }

    #[inline]
    fn remove_layer_and_post_transition(&mut self, layer: impl PhmmIndex) {
        self.core_mut().remove_layer_and_post_transition(layer);
    }

    #[inline]
    fn remove_layers_and_pre_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        self.core_mut().remove_layers_and_pre_transitions(range);
    }

    #[inline]
    fn remove_layers_and_post_transitions<R: PhmmIndexRange>(&mut self, range: R) {
        self.core_mut().remove_layers_and_post_transitions(range);
    }
}
