//! Implementations of the forward algorithm on pHMMs for sequence alignment.

use crate::alignment::phmm::{
    PhmmNumber,
    components::LayerParams,
    state::PhmmState::{Delete, Insert, Match},
    sum_probs,
};

mod global;

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the parameter for the next insert state and the state.
#[inline]
#[must_use]
fn update_insert<T: PhmmNumber, const S: usize>(
    layer: &LayerParams<T, S>, x_idx: usize, match_val: T, delete_val: T, insert_val: T,
) -> T {
    layer.emission_insert[x_idx]
        + sum_probs([
            match_val + layer.transition[(Match, Insert)],
            delete_val + layer.transition[(Delete, Insert)],
            insert_val + layer.transition[(Insert, Insert)],
        ])
}

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the parameter for the next delete state.
#[inline]
#[must_use]
fn update_delete<T: PhmmNumber, const S: usize>(layer: &LayerParams<T, S>, match_val: T, delete_val: T, insert_val: T) -> T {
    sum_probs([
        match_val + layer.transition[(Match, Delete)],
        delete_val + layer.transition[(Delete, Delete)],
        insert_val + layer.transition[(Insert, Delete)],
    ])
}
