use crate::alignment::{
    Alignment, AlignmentStates,
    phmm::{
        DomainPhmm, GetLayer, GetModule, LayerParams, PhmmBacktrackFlags, PhmmError, PhmmNumber,
        PhmmState::{self, Delete, Insert, Match},
        PhmmTracebackState, best_state,
        indexing::{DpIndex, LastBase, LastMatch, PhmmIndexable, QueryIndex, QueryIndexable},
        modules::PrecomputedDomainModule,
        viterbi::{ViterbiTraceback, update_delete, update_insert},
    },
};
use std::ops::Bound::{Excluded, Included};

/// A tracker for the best score for the domain Viterbi algorithm.
struct DomainBestScore<T> {
    /// The best score found at the end of the model
    score: T,
    /// The last query index consumed by the [`CorePhmm`]
    ///
    /// [`CorePhmm`]: crate::alignment::phmm::CorePhmm
    i:     DpIndex,
    /// The state from which the END state was reached
    state: PhmmState,
}

impl<T: PhmmNumber> DomainBestScore<T> {
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn update_last_layer<const S: usize>(
        &mut self, layer: &LayerParams<T, S>, mut match_val: T, mut delete_val: T, mut insert_val: T, i: impl QueryIndex,
        seq: &[u8], end: &PrecomputedDomainModule<T, S>,
    ) {
        use crate::alignment::phmm::PhmmState::*;

        match_val += layer.transition[(Match, Match)];
        delete_val += layer.transition[(Delete, Match)];
        insert_val += layer.transition[(Insert, Match)];

        let (state, mut score) = best_state(match_val, delete_val, insert_val);
        score += end.get_score(i);

        if score < self.score {
            self.score = score;
            self.i = seq.to_dp_index(i);
            self.state = state;
        }
    }
}

impl<T: PhmmNumber> Default for DomainBestScore<T> {
    #[inline]
    fn default() -> Self {
        Self {
            score: T::INFINITY,
            i:     DpIndex(0),
            state: PhmmState::Match,
        }
    }
}

impl<T: PhmmNumber, const S: usize> DomainPhmm<T, S> {
    /// Computes the best scoring domain alignment along with its score via the
    /// Viterbi algorithm.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found, an error is given. An
    /// error is also returned if the model has no layers.
    #[allow(clippy::too_many_lines)]
    pub fn viterbi<Q: AsRef<[u8]>>(&self, seq: Q) -> Result<Alignment<T>, PhmmError> {
        let seq = seq.as_ref();

        let begin_mod = self.begin().precompute_begin_mod(seq, self.mapping());
        let end_mod = self.end().precompute_end_mod(seq, self.mapping());

        // Check for compatibility between the pHMM and the start/end modules
        if begin_mod.seq_len() != seq.len() || end_mod.seq_len() != seq.len() {
            return Err(PhmmError::IncompatibleModule);
        }

        let [layers @ .., end] = self.core().layers() else {
            return Err(PhmmError::EmptyModel);
        };

        let query_dim = seq.len() + 1;
        // This is equivalent to self.seq_len()+1, but may help with bounds
        // checking
        let phmm_dim = layers.len() + 1;

        let mut v_m = vec![T::INFINITY; query_dim];
        for (i, value) in v_m.iter_mut().enumerate() {
            *value = begin_mod.get_score(DpIndex(i));
        }
        let mut v_i = vec![T::INFINITY; query_dim];
        let mut v_d = vec![T::INFINITY; query_dim];

        let mut j = 0;

        let mut traceback = ViterbiTraceback::new(PhmmBacktrackFlags::new(), query_dim, phmm_dim);

        let mut best_score: DomainBestScore<T> = DomainBestScore::default();

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
            // No version of alignment can enter a match state after BEGIN
            // without consuming a base
            v_m[0] = T::INFINITY;
        }

        // Calculate last column of insertion table and previous states for END
        // state
        let start = query_dim * j;
        let traceback_curr_row = &mut traceback.data[start..start + query_dim];
        for (i, x_idx) in seq.iter().map(|x| self.mapping().to_index(*x)).enumerate() {
            best_score.update_last_layer(end, v_m[i], v_d[i], v_i[i], DpIndex(i), seq, &end_mod);

            let (state_i, insert_score) = update_insert(end, x_idx, v_m[i], v_d[i], v_i[i]);
            traceback_curr_row[i + 1].set_insert(state_i);
            v_i[i + 1] = insert_score;
        }

        let i = seq.len();
        best_score.update_last_layer(end, v_m[i], v_d[i], v_i[i], LastBase, seq, &end_mod);

        // This is a necessary check, otherwise the traceback may panic
        if best_score.score == T::INFINITY {
            return Err(PhmmError::NoAlignmentFound);
        }

        let DomainBestScore {
            score,
            i: DpIndex(end_i),
            state,
        } = best_score;
        let end_j = self.get_dp_index(LastMatch);

        let mut state = PhmmTracebackState::from(state);
        let mut i = end_i;
        let mut j = end_j;

        let mut states = AlignmentStates::new();
        states.soft_clip(seq.len() - i);

        // Continue traceback until BEGIN state is reached
        while j > 0 || !state.is_match() {
            let next_state = traceback.get(i, j).get_prev_state(state);

            if state.is_match() {
                states.add_state(b'M');
                i -= 1;
                j -= 1;
            } else if state.is_delete() {
                states.add_state(b'D');
                j -= 1;
            } else {
                states.add_state(b'I');
                i -= 1;
            }

            state = next_state;
        }

        states.soft_clip(i);

        states.make_reverse();

        // start is excluded since we subtracted 1 at the end of the loop, but
        // end is included because that was where we started
        let start_j = Excluded(DpIndex(j));
        let end_j = Included(DpIndex(end_j));
        let start_i = Excluded(DpIndex(i));
        let end_i = Included(DpIndex(end_i));

        Ok(Alignment {
            score,
            ref_range: self.get_seq_range((start_j, end_j)),
            query_range: seq.get_seq_range(start_i, end_i),
            states,
            ref_len: self.seq_len(),
            query_len: seq.len(),
        })
    }
}
