use crate::alignment::phmm::{LayerParams, PhmmNumber, PhmmState, PhmmStateOrEnter, best_state, indexing::SeqIndex};

mod domain;
mod global;
mod local;
mod semilocal;

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the best score for the next insert state and the state from which
/// it was reached.
#[inline]
#[must_use]
fn update_insert<T: PhmmNumber, const S: usize>(
    layer: &LayerParams<T, S>, x_idx: usize, mut match_val: T, mut delete_val: T, mut insert_val: T,
) -> (PhmmState, T) {
    use crate::alignment::phmm::PhmmState::*;

    match_val += layer.transition[(Match, Insert)];
    delete_val += layer.transition[(Delete, Insert)];
    insert_val += layer.transition[(Insert, Insert)];

    let (state, best) = best_state(match_val, delete_val, insert_val);
    (state, best + layer.emission_insert[x_idx])
}

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the best score for the next delete state and the state from which
/// it was reached.
#[inline]
#[must_use]
fn update_delete<T: PhmmNumber, const S: usize>(
    layer: &LayerParams<T, S>, mut match_val: T, mut delete_val: T, mut insert_val: T,
) -> (PhmmState, T) {
    use crate::alignment::phmm::PhmmState::*;

    match_val += layer.transition[(Match, Delete)];
    delete_val += layer.transition[(Delete, Delete)];
    insert_val += layer.transition[(Insert, Delete)];

    best_state(match_val, delete_val, insert_val)
}

/// A traceback matrix for Viterbi alignment.
///
/// The data is arranged with the query residues as the columns and the pHMM
/// layers as the rows, and then row-major order is used.
///
/// The number of columns is one more than the query length (to allow the
/// condition of no bases in the query being matched to be represented). The
/// number of rows is one more than the number of pHMM match states (i.e., the
/// number of match states along with the BEGIN state).
struct ViterbiTraceback<T> {
    data: Vec<T>,
    cols: usize,
}

impl<T: Clone> ViterbiTraceback<T> {
    /// Creates a new [`ViterbiTraceback`] from the given dimensions.
    ///
    /// `query_dim` should be the length of the query plus one, since the
    /// traceback stores a column for [`NoBases`] as well as all [`SeqIndex`]
    /// values. `phmm_dim` should be the [`seq_len`] of the pHMM plus one, since
    /// the traceback stores a row for [`Begin`] as well as all [`SeqIndex`]
    /// values (but not for [`End`]).
    ///
    /// [`NoBases`]: super::indexing::NoBases
    /// [`seq_len`]: super::indexing::PhmmIndexable::seq_len
    /// [`Begin`]: super::indexing::Begin
    /// [`End`]: super::indexing::End
    #[inline]
    #[must_use]
    pub fn new(default: T, query_dim: usize, phmm_dim: usize) -> Self {
        Self {
            data: vec![default; query_dim * phmm_dim],
            cols: query_dim,
        }
    }
}

impl<T> ViterbiTraceback<T> {
    /// Retrieves the state stored in the traceback given the dynamic
    /// programming index.
    ///
    /// `i` is the dynamic programming query index, so that `0` corresponds to
    /// matching no residues, `1` corresponds to the first residue in the query,
    /// and so on.
    ///
    /// `j` is the dynamic programming layer index, so that `0` corresponds to
    /// the BEGIN state, `1` corresponds to the first layer, and so on.
    #[inline]
    #[must_use]
    pub fn get(&self, i: usize, j: usize) -> &T {
        &self.data[self.cols * j + i]
    }

    /// Retrieves a mutable reference to the state stored in the traceback given
    /// the dynamic programming index.
    ///
    /// `i` is the dynamic programming query index, so that `0` corresponds to
    /// matching no residues, `1` corresponds to the first residue in the query,
    /// and so on.
    ///
    /// `j` is the dynamic programming layer index, so that `0` corresponds to
    /// the BEGIN state, `1` corresponds to the first layer, and so on.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        &mut self.data[self.cols * j + i]
    }

    /// Retrieves a mutable reference to a row stored in the traceback given the
    /// dynamic programming index of the pHMM layer.
    ///
    /// `0` corresponds to the BEGIN state, `1` corresponds to the first layer,
    /// and so on.
    #[inline]
    #[must_use]
    #[allow(dead_code)]
    pub fn get_rows_mut(&mut self, j: usize) -> (&mut [T], &mut [T]) {
        let start = self.cols * j;
        self.data[start..2 * self.cols].split_at_mut(start + self.cols)
    }
}

/// The location from which the alignment exits the [`CorePhmm`] for local or
/// semilocal alignment.
///
/// [`CorePhmm`]: super::CorePhmm
#[derive(Clone, Copy)]
enum ExitLocation {
    /// The alignment exited the core pHMM through the BEGIN state.
    Begin,
    /// The alignment exited the core pHMM through a match state (not BEGIN or
    /// END). The index of the match state is stored.
    Match(SeqIndex),
    /// The alignment exited the core pHMM through the END state. The transition
    /// used to reach the END state is stored.
    End(PhmmStateOrEnter),
}
