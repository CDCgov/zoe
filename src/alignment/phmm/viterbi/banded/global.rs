use crate::alignment::{
    Alignment, AlignmentStates,
    phmm::{
        GlobalPhmm, PhmmNumber,
        components::LayerParams,
        indexing::{AlnIndexable, GetLayer},
        state::{
            PhmmState::{self, Delete, Insert, Match},
            PhmmTracebackState, best_state,
        },
        viterbi::{ViterbiError, update_delete, update_insert},
    },
};

use super::BandedViterbiTraceback;

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

/// Returns the half-open query-column range in the band for DP row `j`.
#[inline]
fn band_bounds(j: usize, band_width: usize, query_dim: usize) -> (usize, usize) {
    let start_col = j.saturating_sub(band_width);
    let end_col = j.saturating_add(band_width).saturating_add(1).min(query_dim);
    (start_col, end_col)
}

impl<T: PhmmNumber, const S: usize> GlobalPhmm<T, S> {
    /// Computes the best scoring global alignment contained within a diagonal
    /// band using the Viterbi algorithm.
    ///
    /// `band_width` is the maximum difference between the number of pHMM
    /// layers and query residues consumed by an alignment path. Consequently,
    /// a band narrower than the difference between the model and query lengths
    /// cannot contain a global alignment.
    ///
    /// This is exact within the band, but may differ from [`Self::viterbi`] if
    /// the unrestricted best path leaves the band.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found within the band, an
    /// error is given.
    pub fn viterbi_banded<Q: AsRef<[u8]>>(&self, seq: Q, band_width: usize) -> Result<Alignment<T>, ViterbiError> {
        let seq = seq.as_ref();
        let (end, layers) = self.layers().split_last();

        let query_dim = seq.len() + 1;
        // This is equivalent to self.seq_len() + 1.
        let phmm_dim = layers.len() + 1;

        let mut v_m = vec![T::INFINITY; query_dim];
        v_m[0] = T::ZERO;
        let mut v_i = vec![T::INFINITY; query_dim];
        let mut v_d = vec![T::INFINITY; query_dim];

        let mut traceback = BandedViterbiTraceback::new(phmm_dim, band_width);

        for (j, layer) in layers.iter().enumerate() {
            let (start_col, end_col) = band_bounds(j, band_width, query_dim);
            if start_col >= end_col {
                return Err(ViterbiError::NoAlignmentFound);
            }

            // An insertion at the first column in this row would have to come
            // from outside the band. This also discards the insertion score
            // belonging to the preceding pHMM layer.
            v_i[start_col] = T::INFINITY;

            // The match score at the current column is replaced with the score
            // for the next row, so preserve it as the next diagonal value.
            let mut cur_m = v_m[start_col];
            let next_start_col = (j + 1).saturating_sub(band_width);

            for i in start_col..end_col {
                let match_val = cur_m;
                let delete_val = v_d[i];
                let insert_val = v_i[i];

                if i < seq.len() {
                    let x_idx = self.mapping().to_index(seq[i]);
                    let (state_m, match_score) = {
                        let (state, best) = best_state(
                            match_val + layer.transition[(Match, Match)],
                            delete_val + layer.transition[(Delete, Match)],
                            insert_val + layer.transition[(Insert, Match)],
                        );
                        (state, best + layer.emission_match[x_idx])
                    };
                    traceback.get_mut(i + 1, j + 1).set_match(state_m);
                    cur_m = std::mem::replace(&mut v_m[i + 1], match_score);

                    // Insertion remains on the current row, so its destination
                    // must also be inside this row's band.
                    if i + 1 < end_col {
                        let (state_i, insert_score) = update_insert(layer, x_idx, match_val, delete_val, insert_val);
                        traceback.get_mut(i + 1, j).set_insert(state_i);
                        v_i[i + 1] = insert_score;
                    }
                }

                // Delete advances to the next row without consuming a query
                // residue. At the lower edge it can leave the band.
                if i >= next_start_col {
                    let (state_d, delete_score) = update_delete(layer, match_val, delete_val, insert_val);
                    traceback.get_mut(i, j + 1).set_delete(state_d);
                    v_d[i] = delete_score;
                }
            }

            // While the band still touches the first query column, no match
            // state there can be reached after BEGIN.
            if next_start_col == 0 {
                v_m[0] = T::INFINITY;
            }
        }

        let j = layers.len();
        let (start_col, end_col) = band_bounds(j, band_width, query_dim);

        // Global alignment must finish at the bottom-right DP coordinate.
        if !(start_col..end_col).contains(&seq.len()) {
            return Err(ViterbiError::NoAlignmentFound);
        }

        // Calculate insertions in the final pHMM layer. As above, the first
        // insertion cell cannot be reached horizontally from outside the band.
        v_i[start_col] = T::INFINITY;
        for i in start_col..seq.len() {
            let x_idx = self.mapping().to_index(seq[i]);
            let (state_i, insert_score) = update_insert(end, x_idx, v_m[i], v_d[i], v_i[i]);
            traceback.get_mut(i + 1, j).set_insert(state_i);
            v_i[i + 1] = insert_score;
        }

        let i = seq.len();
        let mut best_score = GlobalBestScore::default();
        best_score.update_seq_end_last_layer(end, v_m[i], v_d[i], v_i[i]);

        if best_score.score == T::INFINITY {
            return Err(ViterbiError::NoAlignmentFound);
        }

        let GlobalBestScore { state, score } = best_score;
        let mut state = PhmmTracebackState::from(state);
        let (mut i, mut j) = (seq.len(), layers.len());
        let mut states = AlignmentStates::new();

        while i > 0 || j > 0 {
            let next_state = traceback.get(i, j).get_prev_state(state);

            if state.is_match() {
                debug_assert!(i > 0 && j > 0);
                states.add_state(b'M');
                i -= 1;
                j -= 1;
            } else if state.is_delete() {
                debug_assert!(j > 0);
                states.add_state(b'D');
                j -= 1;
            } else {
                debug_assert!(i > 0);
                states.add_state(b'I');
                i -= 1;
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
