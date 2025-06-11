use crate::{
    alignment::{
        Alignment,
        phmm::{
            CorePhmm, DpIndex, EnumArray, LayerParams, PhmmError, PhmmIndex, PhmmNumber, PhmmState, PhmmStateArray,
            PhmmStateEnum, PhmmStateOrEnter, QueryIndex, SeqIndex,
        },
    },
    data::ByteIndexMap,
};

mod domain;
mod global;
mod local;
mod semilocal;

pub use domain::*;
pub use global::*;
pub use local::*;
pub use semilocal::*;

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the best score for the next match state and the state from which
/// it was reached.
///
/// `vals` may also contain any additional transitions, such as the "enter"
/// score for the case of local alignment.
#[inline]
#[must_use]
fn update_match<T: PhmmNumber, const S: usize, E: PhmmStateEnum, const N: usize>(
    layer: &LayerParams<T, S>, x_idx: usize, mut vals: EnumArray<T, E, N>,
) -> (E, T) {
    use crate::alignment::phmm::state::PhmmState::*;

    vals[Delete.into()] += layer.transition[(Delete, Match)];
    vals[Match.into()] += layer.transition[(Match, Match)];
    vals[Insert.into()] += layer.transition[(Insert, Match)];

    let (state, best) = vals.locate_min();
    (state, best + layer.emission_match[x_idx])
}

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the best score for the next insert state and the state from which
/// it was reached.
#[inline]
#[must_use]
fn update_insert<T: PhmmNumber, const S: usize>(
    layer: &LayerParams<T, S>, x_idx: usize, mut vals: PhmmStateArray<T>,
) -> (PhmmState, T) {
    use crate::alignment::phmm::PhmmState::*;

    vals[Delete] += layer.transition[(Delete, Insert)];
    vals[Match] += layer.transition[(Match, Insert)];
    vals[Insert] += layer.transition[(Insert, Insert)];

    let (state, best) = vals.locate_min();
    (state, best + layer.emission_insert[x_idx])
}

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the best score for the next delete state and the state from which
/// it was reached.
#[inline]
#[must_use]
fn update_delete<T: PhmmNumber, const S: usize>(layer: &LayerParams<T, S>, mut vals: PhmmStateArray<T>) -> (PhmmState, T) {
    use crate::alignment::phmm::PhmmState::*;

    vals[Delete] += layer.transition[(Delete, Delete)];
    vals[Match] += layer.transition[(Match, Delete)];
    vals[Insert] += layer.transition[(Insert, Delete)];

    vals.locate_min()
}

/// A traceback matrix for Viterbi alignment.
///
/// The data is arranged with the query residues as the columns and the PHMM
/// layers as the rows, and then row-major order is used.
///
/// The number of columns is one more than the query length (to allow the
/// condition of no bases in the query being matched to be represented). The
/// number of rows is one more than the number of PHMM match states (i.e., the
/// number of match states along with the BEGIN state).
pub struct ViterbiTraceback<T> {
    data: Vec<T>,
    cols: usize,
}

impl<T: Clone> ViterbiTraceback<T> {
    /// Creates a new [`ViterbiTraceback`] from the given query length and the
    /// layer length.
    ///
    /// The layer length is the number of match states excluding BEGIN and END.
    #[inline]
    #[must_use]
    fn new(default: T, query_len: usize, layers_len: usize) -> Self {
        Self {
            data: vec![default; (query_len + 1) * (layers_len + 1)],
            cols: query_len + 1,
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
    pub fn get_rows_mut(&mut self, j: usize) -> (&mut [T], &mut [T]) {
        let start = self.cols * j;
        self.data[start..2 * self.cols].split_at_mut(start + self.cols)
    }
}

/// A strategy for how to properly run the Viterbi algorithm for different
/// alignment modes. By specifying various components of the algorithm via the
/// required trait methods, this trait provides an implementation of the Viterbi
/// algorithm.
pub(crate) trait ViterbiStrategy<'a, T: PhmmNumber + 'a, const S: usize>: Sized {
    /// The type used in the traceback matrix, either [`PhmmState`] or
    /// [`PhmmStateOrEnter`].
    ///
    /// [`PhmmStateOrEnter`]: super::PhmmStateOrEnter
    type TracebackState: From<PhmmState> + Copy;

    /// The type used for tracking the best score and end of the alignment.
    type BestScore: Default + BestScore<T, S, Strategy<'a> = Self>;

    /// The default value to initialize the traceback matrix with.
    const TRACEBACK_DEFAULT: PhmmStateArray<Self::TracebackState>;

    /// Retrieves the underlying core pHMM.
    fn core(&self) -> &CorePhmm<T, S>;

    /// Retrieves the byte mapping used to represent the alphabet.
    fn mapping(&self) -> &ByteIndexMap<S>;

    /// Retrieves the query.
    fn query(&self) -> &[u8];

    /// Performs any initial checks on the pHMM
    fn initial_check(&self) -> Result<(), PhmmError>;

    /// Fills `v_m` with the initialization required for the first column (j=0).
    fn fill_vm(&self, v_m: &mut [T]);

    /// Updates the score for entering the match state.
    fn update_match_score(
        &self, layer: &LayerParams<T, S>, x_idx: usize, cur_vals: PhmmStateArray<T>, i: impl QueryIndex, j: impl PhmmIndex,
    ) -> (Self::TracebackState, T);

    /// Performs the traceback.
    fn perform_traceback(
        self, best_score: Self::BestScore, traceback: ViterbiTraceback<PhmmStateArray<Self::TracebackState>>,
    ) -> Alignment<T>;

    /// Obtains the best scoring alignment along with its score via the Viterbi
    /// algorithm.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found, an error is given. An
    /// error is also returned if the model has no layers.
    #[allow(clippy::range_plus_one)]
    fn viterbi(self) -> Result<Alignment<T>, PhmmError> {
        use crate::alignment::phmm::PhmmState::*;

        self.initial_check()?;

        let [layers @ .., end] = self.core().0.as_slice() else {
            return Err(PhmmError::EmptyModel);
        };

        let query_dim = self.query().len() + 1;
        let layer_dim = layers.len() + 1;

        let mut v_m = vec![T::INFINITY; query_dim];
        self.fill_vm(&mut v_m);
        let mut v_i = vec![T::INFINITY; query_dim];
        let mut v_d = vec![T::INFINITY; query_dim];

        let mut j = 0;

        let mut traceback = ViterbiTraceback {
            data: vec![Self::TRACEBACK_DEFAULT; query_dim * layer_dim],
            cols: query_dim,
        };

        let mut best_score = Self::BestScore::default();

        for layer in layers {
            let mut cur_m = v_m[0];

            let start = query_dim * j;
            let (traceback_curr_row, traceback_next_row) =
                traceback.data[start..start + 2 * query_dim].split_at_mut(query_dim);

            let traceback_curr_row = &mut traceback_curr_row[0..self.query().len() + 1];
            let traceback_next_row = &mut traceback_next_row[0..self.query().len() + 1];

            for (i, x_idx) in self.query().iter().map(|x| self.mapping().to_index(*x)).enumerate() {
                let cur_vals = PhmmStateArray::new([v_d[i], cur_m, v_i[i]]);
                best_score.update(&self, cur_vals, DpIndex(i), DpIndex(j));

                let (state_m, match_score) = self.update_match_score(layer, x_idx, cur_vals, DpIndex(i), DpIndex(j));
                traceback_next_row[i + 1][Match] = state_m;
                cur_m = std::mem::replace(&mut v_m[i + 1], match_score);

                let (state_i, insert_score) = update_insert(layer, x_idx, cur_vals);
                traceback_curr_row[i + 1][Insert] = state_i.into();
                v_i[i + 1] = insert_score;

                let (state_d, delete_score) = update_delete(layer, cur_vals);
                traceback_next_row[i][Delete] = state_d.into();
                v_d[i] = delete_score;
            }

            let i = self.query().len();
            let cur_vals = PhmmStateArray::new([v_d[i], cur_m, v_i[i]]);
            best_score.update_seq_end(&self, cur_vals, DpIndex(j));

            let (state_d, delete_score) = update_delete(layer, cur_vals);
            traceback_next_row[i][Delete] = state_d.into();
            v_d[i] = delete_score;

            j += 1;
            // No version of alignment can enter a match state after BEGIN
            // without consuming a base
            v_m[0] = T::INFINITY;
        }

        // Calculate last column of insertion table and previous states for END
        // state
        let start = query_dim * j;
        let traceback_curr_row = &mut traceback.data[start..start + self.query().len() + 1];
        for (i, x_idx) in self.query().iter().map(|x| self.mapping().to_index(*x)).enumerate() {
            let cur_vals = PhmmStateArray::new([v_d[i], v_m[i], v_i[i]]);
            best_score.update_last_layer(&self, end, cur_vals, DpIndex(i));

            let (state_i, insert_score) = update_insert(end, x_idx, cur_vals);
            traceback_curr_row[i + 1][Insert] = state_i.into();
            v_i[i + 1] = insert_score;
        }

        let i = self.query().len();
        let cur_vals = PhmmStateArray::new([v_d[i], v_m[i], v_i[i]]);
        best_score.update_seq_end_last_layer(&self, end, cur_vals);

        // This is a necessary check, otherwise the traceback may panic
        if best_score.score() == T::INFINITY {
            return Err(PhmmError::NoAlignmentFound);
        }

        Ok(self.perform_traceback(best_score, traceback))
    }
}

/// Specifications for how to track the best score for the Viterbi algorithm.
///
/// This is designed to encapsulate behavior between the various types of
/// alignment. Four hooks are provided to allow the best score to update at
/// various points in the Viterbi algorithm.
pub(crate) trait BestScore<T, const S: usize> {
    /// The [`ViterbiStrategy`] struct this [`BestScore`] type corresponds with.
    type Strategy<'a>
    where
        T: 'a;

    /// Retrieves the current best score.
    fn score(&self) -> T;

    /// Updates the best score before the end of the model and sequence.
    ///
    /// ## Arguments
    ///
    /// * `strategy`: the [`ViterbiStrategy`] struct this [`BestScore`] type
    ///   corresponds with
    /// * `vals`: the current values in the delete, match, and insert states
    /// * `i`: The current number of bases consumed
    /// * `j`: The current model layer
    fn update(&mut self, _strategy: &Self::Strategy<'_>, _vals: PhmmStateArray<T>, _i: impl QueryIndex, _j: impl PhmmIndex) {
    }

    /// Updates the best score before the end of the model but at the end of the
    /// sequence.
    ///
    /// ## Arguments
    ///
    /// * `strategy`: the [`ViterbiStrategy`] struct this [`BestScore`] type
    ///   corresponds with
    /// * `vals`: the current values in the delete, match, and insert states
    /// * `j`: The current model layer
    fn update_seq_end(&mut self, _strategy: &Self::Strategy<'_>, _vals: PhmmStateArray<T>, _j: impl PhmmIndex) {}

    /// Updates the best score before the end of the sequence but at the end of
    /// the model.
    ///
    /// ## Arguments
    ///
    /// * `strategy`: the [`ViterbiStrategy`] struct this [`BestScore`] type
    ///   corresponds with
    /// * `layer`: the last layer, including transition probabilities into the
    ///   END state
    /// * `vals`: the values in the delete, match, and insert states of the
    ///   layer before the END state. This function must handle early exits out
    ///   of that state as well as exits through the END state.
    /// * `i`: The current number of bases consumed
    fn update_last_layer(
        &mut self, _strategy: &Self::Strategy<'_>, _layer: &LayerParams<T, S>, _vals: PhmmStateArray<T>, _i: impl QueryIndex,
    ) {
    }

    /// Update the best score at the end of the model and sequence.
    ///
    /// ## Arguments
    ///
    /// * `strategy`: the [`ViterbiStrategy`] struct this [`BestScore`] type
    ///   corresponds with
    /// * `layer`: the last layer, including transition probabilities into the
    ///   END state
    /// * `vals`: the current values in the delete, match, and insert states
    fn update_seq_end_last_layer(
        &mut self, _strategy: &Self::Strategy<'_>, _layer: &LayerParams<T, S>, _vals: PhmmStateArray<T>,
    ) {
    }
}

/// The location from which the alignment exits the [`CorePhmm`] for local or
/// semilocal alignment.
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
