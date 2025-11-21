use crate::alignment::{
    Alignment, AlignmentStates,
    phmm::{
        GetLayer, GetModule, LayerParams, PhmmBacktrackFlags, PhmmError, PhmmNumber,
        PhmmState::{self, Delete, Insert, Match},
        PhmmStateOrEnter, PhmmTracebackState, SemiLocalPhmm, best_state, best_state_or_enter,
        indexing::{Begin, DpIndex, End, LastMatch, NoBases, PhmmIndex, PhmmIndexable, QueryIndexable},
        viterbi::{ExitLocation, ViterbiTraceback, update_delete, update_insert},
    },
};
use std::ops::Bound::{Excluded, Included};

/// A tracker for the best score for the semilocal Viterbi algorithm.
struct SemiLocalBestScore<T> {
    /// The best score found at the end of the model
    score: T,
    /// The location from which the alignment exits the [`CorePhmm`]
    ///
    /// [`CorePhmm`]: crate::alignment::phmm::CorePhmm
    loc:   ExitLocation,
}

impl<T: PhmmNumber> SemiLocalBestScore<T> {
    #[inline]
    fn update_seq_end<const S: usize>(&mut self, match_val: T, j: impl PhmmIndex, phmm: &SemiLocalPhmm<T, S>) {
        let score = match_val + phmm.end().get_score(j);
        if score < self.score {
            self.score = score;
            self.loc = match phmm.to_seq_index(j) {
                Some(loc) => ExitLocation::Match(loc),
                None => ExitLocation::Begin,
            }
        }
    }

    #[inline]
    fn update_seq_end_last_layer<const S: usize>(
        &mut self, layer: &LayerParams<T, S>, mut match_val: T, mut delete_val: T, mut insert_val: T,
        phmm: &SemiLocalPhmm<T, S>, seq: &[u8],
    ) {
        use crate::alignment::phmm::PhmmState::*;

        // Option 1: Early exit from this layer
        self.update_seq_end(match_val, LastMatch, phmm);

        // Option 2: Go through END state
        match_val += layer.transition[(Match, Match)];
        delete_val += layer.transition[(Delete, Match)];
        insert_val += layer.transition[(Insert, Match)];

        let (state, mut score) = if seq.is_empty() {
            let enter_val = phmm.begin().get_score(End);
            best_state_or_enter(match_val, delete_val, insert_val, enter_val)
        } else {
            let (state, score) = best_state(match_val, delete_val, insert_val);
            (PhmmStateOrEnter::from(state), score)
        };

        score += phmm.end().get_score(End);

        if score < self.score {
            self.score = score;
            self.loc = ExitLocation::End(state);
        }
    }
}

impl<T: PhmmNumber> Default for SemiLocalBestScore<T> {
    #[inline]
    fn default() -> Self {
        Self {
            score: T::INFINITY,
            loc:   ExitLocation::End(PhmmStateOrEnter::Match),
        }
    }
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<T, S> {
    /// Computes the best scoring local alignment along with its score via the
    /// Viterbi algorithm.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found, an error is given. An
    /// error is also returned if the model has no layers.
    #[allow(clippy::too_many_lines)]
    pub fn viterbi<Q: AsRef<[u8]>>(&self, seq: Q) -> Result<Alignment<T>, PhmmError> {
        let seq = seq.as_ref();

        // Check for compatibility between the pHMM and the start/end modules
        if self.begin().num_pseudomatch() != self.num_pseudomatch() || self.end().num_pseudomatch() != self.num_pseudomatch()
        {
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
        v_m[0] = self.begin().get_score(Begin);
        let mut v_i = vec![T::INFINITY; query_dim];
        let mut v_d = vec![T::INFINITY; query_dim];

        let mut j = 0;

        let mut traceback = ViterbiTraceback::new(PhmmBacktrackFlags::new(), query_dim, phmm_dim);

        let mut best_score = SemiLocalBestScore::default();

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
                    // We are updating the match score of the layer after j after consuming
                    // one more base after i. To reach this cell via entering, one must
                    // consume all bases up to i in the begin module, then the (i+1)st is
                    // consumed in this match state. The emission parameter is added within
                    // `update_match`.
                    if i == seq.get_dp_index(NoBases) {
                        let next_layer_idx = DpIndex(j).next_index(self);

                        let (state, best) = best_state_or_enter(
                            match_val + layer.transition[(Match, Match)],
                            delete_val + layer.transition[(Delete, Match)],
                            insert_val + layer.transition[(Insert, Match)],
                            self.begin().get_score(next_layer_idx),
                        );
                        (state, best + layer.emission_match[x_idx])
                    } else {
                        let (state, best) = best_state(
                            match_val + layer.transition[(Match, Match)],
                            delete_val + layer.transition[(Delete, Match)],
                            insert_val + layer.transition[(Insert, Match)],
                        );
                        (PhmmStateOrEnter::from(state), best + layer.emission_match[x_idx])
                    }
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
            best_score.update_seq_end(cur_m, DpIndex(j), self);

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
            let (state_i, insert_score) = update_insert(end, x_idx, v_m[i], v_d[i], v_i[i]);
            traceback_curr_row[i + 1].set_insert(state_i);
            v_i[i + 1] = insert_score;
        }

        let i = seq.len();
        best_score.update_seq_end_last_layer(end, v_m[i], v_d[i], v_i[i], self, seq);

        // This is a necessary check, otherwise the traceback may panic
        if best_score.score == T::INFINITY {
            return Err(PhmmError::NoAlignmentFound);
        }

        let SemiLocalBestScore { score, loc } = best_score;
        let (state, end_j) = match loc {
            ExitLocation::Begin => {
                let mut states = AlignmentStates::new();
                states.soft_clip(seq.len());
                return Ok(Alignment {
                    score,
                    ref_range: 0..0,
                    query_range: 0..0,
                    states,
                    ref_len: self.seq_len(),
                    query_len: seq.len(),
                });
            }
            ExitLocation::Match(j) => (Match, self.get_dp_index(j)),
            ExitLocation::End(ptr) => {
                let Some(ptr) = PhmmState::get_from(ptr) else {
                    let mut states = AlignmentStates::new();
                    states.soft_clip(seq.len());
                    return Ok(Alignment {
                        score,
                        ref_range: self.get_seq_range(End..End),
                        query_range: 0..0,
                        states,
                        ref_len: self.seq_len(),
                        query_len: seq.len(),
                    });
                };
                (ptr, self.get_dp_index(LastMatch))
            }
        };
        let end_i = seq.len();

        let mut state = PhmmTracebackState::from(state);
        let mut i = end_i;
        let mut j = end_j;

        let mut states = AlignmentStates::new();
        states.soft_clip(seq.len() - i);

        // Continue traceback until BEGIN state is reached, or until Enter is
        // encountered (see break statement)
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

            if next_state.is_enter() {
                break;
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
