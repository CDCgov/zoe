use super::ViterbiTraceback;
use crate::alignment::{
    Alignment, AlignmentStates,
    phmm::{
        GetLayer, GlobalPhmm, LayerParams, PhmmBacktrackFlags, PhmmError, PhmmNumber,
        PhmmState::{self, Delete, Insert, Match},
        PhmmTracebackState, best_state,
        indexing::PhmmIndexable,
        viterbi::{update_delete, update_insert},
    },
};

/// A tracker for the best score for the global Viterbi algorithm.
struct GlobalBestScore<T> {
    /// The state from which the END state was reached
    state: PhmmState,
    /// The score in the END state
    score: T,
}

impl<T: PhmmNumber> GlobalBestScore<T> {
    /// Updates the best state and score from the score values at the end of the
    /// sequence and last layer.
    ///
    /// ## Validity
    ///
    /// NaN values should not be passed within the layers or the score values,
    /// and may result in inconsistent behavior.
    fn update_seq_end_last_layer<const S: usize>(
        &mut self, layer: &LayerParams<T, S>, mut match_val: T, mut delete_val: T, mut insert_val: T,
    ) {
        use crate::alignment::phmm::PhmmState::*;

        match_val += layer.transition[(Match, Match)];
        delete_val += layer.transition[(Delete, Match)];
        insert_val += layer.transition[(Insert, Match)];

        (self.state, self.score) = best_state(match_val, delete_val, insert_val);
    }
}

impl<T: PhmmNumber> Default for GlobalBestScore<T> {
    #[inline]
    fn default() -> Self {
        Self {
            state: PhmmState::Match,
            score: T::INFINITY,
        }
    }
}

impl<T: PhmmNumber, const S: usize> GlobalPhmm<T, S> {
    /// Obtain the best scoring global alignment along with its score via the
    /// Viterbi algorithm.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found, an error is given. An
    /// error is also returned if the model has no layers.
    pub fn viterbi<Q: AsRef<[u8]>>(&self, seq: Q) -> Result<Alignment<T>, PhmmError> {
        let seq = seq.as_ref();

        let [layers @ .., end] = self.core().layers() else {
            return Err(PhmmError::EmptyModel);
        };

        let query_dim = seq.len() + 1;
        // This is equivalent to self.seq_len()+1, but may help with bounds
        // checking
        let phmm_dim = layers.len() + 1;

        let mut v_m = vec![T::INFINITY; query_dim];
        v_m[0] = T::ZERO;
        let mut v_i = vec![T::INFINITY; query_dim];
        let mut v_d = vec![T::INFINITY; query_dim];

        let mut j = 0;

        let mut traceback = ViterbiTraceback::new(PhmmBacktrackFlags::new(), query_dim, phmm_dim);

        for layer in layers {
            let mut cur_m = v_m[0];

            let start = query_dim * j;
            let (traceback_curr_row, traceback_next_row) =
                traceback.data[start..start + 2 * query_dim].split_at_mut(query_dim);

            let traceback_curr_row = &mut traceback_curr_row[0..query_dim];
            let traceback_next_row = &mut traceback_next_row[0..query_dim];

            for (i, x_idx) in seq.iter().map(|x| self.mapping().to_index(*x)).enumerate() {
                let match_val = cur_m;
                let delete_val = v_d[i];
                let insert_val = v_i[i];

                let (state_m, match_score) = {
                    let (state, best) = best_state(
                        match_val + layer.transition[(Match, Match)],
                        delete_val + layer.transition[(Delete, Match)],
                        insert_val + layer.transition[(Insert, Match)],
                    );
                    (state, best + layer.emission_match[x_idx])
                };
                traceback_next_row[i + 1].set_match(state_m);
                cur_m = std::mem::replace(&mut v_m[i + 1], match_score);

                let (state_i, insert_score) = update_insert(layer, x_idx, match_val, delete_val, insert_val);
                traceback_curr_row[i + 1].set_insert(state_i);
                v_i[i + 1] = insert_score;

                let (state_d, delete_score) = update_delete(layer, match_val, delete_val, insert_val);
                traceback_next_row[i].set_delete(state_d);
                v_d[i] = delete_score;
            }

            let i = seq.len();
            let (state_d, delete_score) = update_delete(layer, cur_m, v_d[i], v_i[i]);
            traceback_next_row[i].set_delete(state_d);
            v_d[i] = delete_score;

            j += 1;
            // Cannot enter a match state after BEGIN without consuming a base
            v_m[0] = T::INFINITY;
        }

        // Calculate last column of insertion table and previous states for END
        // state
        let start = query_dim * j;
        let traceback_curr_row = &mut traceback.data[start..start + query_dim];
        for (i, x_idx) in seq.iter().map(|x| self.mapping().to_index(*x)).enumerate() {
            let (state_i, insert_score) = update_insert(end, x_idx, v_m[i], v_d[i], v_i[i]);
            traceback_curr_row[i + 1].set_insert(state_i);
            v_i[i + 1] = insert_score;
        }

        let i = seq.len();
        let mut best_score = GlobalBestScore::default();
        best_score.update_seq_end_last_layer(end, v_m[i], v_d[i], v_i[i]);

        // This is a necessary check, otherwise the traceback may panic
        if best_score.score == T::INFINITY {
            return Err(PhmmError::NoAlignmentFound);
        }

        let mut cursor = traceback.data.len() - 1;
        let GlobalBestScore { state, score } = best_score;
        let mut state = PhmmTracebackState::from(state);

        let mut states = AlignmentStates::new();

        while cursor > 0 {
            let next_state = traceback.data[cursor].get_prev_state(state);

            if state.is_match() {
                states.add_state(b'M');
                // i -= 1; j -= 1;
                cursor -= traceback.cols + 1;
            } else if state.is_delete() {
                states.add_state(b'D');
                // j -= 1;
                cursor -= traceback.cols;
            } else {
                states.add_state(b'I');
                // i -= 1;
                cursor -= 1;
            }

            state = next_state;
        }

        states.make_reverse();

        Ok(Alignment {
            score,
            ref_range: 0..self.seq_len(),
            query_range: 0..seq.len(),
            states,
            ref_len: self.seq_len(),
            query_len: seq.len(),
        })
    }
}
