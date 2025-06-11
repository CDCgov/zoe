use crate::{
    alignment::{
        Alignment, AlignmentStates,
        phmm::{
            CorePhmm, DpIndex, EnumArray, GlobalPhmm, LayerParams, PhmmError, PhmmIndex, PhmmState, PhmmStateArray,
            PhmmStateEnum, QueryIndex,
        },
    },
    data::ByteIndexMap,
    math::Float,
};

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the best score for the next match state and the state from which
/// it was reached.
///
/// `vals` may also contain any additional transitions, such as the "enter"
/// score for the case of local alignment.
#[inline]
#[must_use]
fn update_match<T: Float, const S: usize, E: PhmmStateEnum, const N: usize>(
    layer: &LayerParams<T, S>, x_idx: usize, mut vals: EnumArray<T, E, N>,
) -> (E, T) {
    use crate::alignment::phmm::state::PhmmState::*;

    vals[Delete.into()] += layer.transition[(Delete, Match)];
    vals[Match.into()] += layer.transition[(Match, Match)];
    vals[Insert.into()] += layer.transition[(Insert, Match)];

    let (ptr, best) = vals.locate_min();
    (ptr, best + layer.emission_match[x_idx])
}

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the best score for the next insert state and the state from which
/// it was reached.
#[inline]
#[must_use]
fn update_insert<T: Float, const S: usize>(
    layer: &LayerParams<T, S>, x_idx: usize, mut vals: PhmmStateArray<T>,
) -> (PhmmState, T) {
    use crate::alignment::phmm::PhmmState::*;

    vals[Delete] += layer.transition[(Delete, Insert)];
    vals[Match] += layer.transition[(Match, Insert)];
    vals[Insert] += layer.transition[(Insert, Insert)];

    let (ptr, best) = vals.locate_min();
    (ptr, best + layer.emission_insert[x_idx])
}

/// Given the current `vals` for delete/match/insert and the next `layer`,
/// calculate the best score for the next delete state and the state from which
/// it was reached.
#[inline]
#[must_use]
fn update_delete<T: Float, const S: usize>(layer: &LayerParams<T, S>, mut vals: PhmmStateArray<T>) -> (PhmmState, T) {
    use crate::alignment::phmm::PhmmState::*;

    vals[Delete] += layer.transition[(Delete, Delete)];
    vals[Match] += layer.transition[(Match, Delete)];
    vals[Insert] += layer.transition[(Insert, Delete)];

    vals.locate_min()
}

/// Specifications for how to properly run the Viterbi algorithm. By specifying
/// various components of the algorithm via the required trait methods, this
/// trait provides an implementation of the Viterbi algorithm.
pub(crate) trait ViterbiSpecs<'a, T: Float + 'a, const S: usize>: Sized {
    /// The pointer type used in the traceback matrix, either [`PhmmState`] or
    /// [`PhmmStateOrEnter`].
    ///
    /// [`PhmmStateOrEnter`]: super::PhmmStateOrEnter
    type Ptr: From<PhmmState>;

    /// The type used for tracking the best score and end of the alignment.
    type BestScore: Default + BestScore<T, S, Specs<'a> = Self>;

    /// Retrieves the underlying core pHMM.
    fn core(&self) -> &CorePhmm<T, S>;

    /// Retrieves the byte mapping used to represent the alphabet.
    fn mapping(&self) -> &ByteIndexMap<S>;

    /// Retrieves the query.
    fn query(&self) -> &[u8];

    /// Performs any initial checks on the pHMM
    fn initial_check(&self) -> Result<(), PhmmError>;

    /// Obtains the initialization for the first row of `v_m` (j=0).
    fn initialize_vm(&self) -> Vec<T>;

    /// Obtains the initialization of the traceback matrix.
    fn initialize_ptrs(&self) -> Vec<Vec<PhmmStateArray<Self::Ptr>>>;

    /// Updates the score for entering the match state.
    fn update_match_score(
        &self, layer: &LayerParams<T, S>, x_idx: usize, cur_vals: PhmmStateArray<T>, i: impl QueryIndex, j: impl PhmmIndex,
    ) -> (Self::Ptr, T);

    /// Performs the traceback.
    fn traceback(self, best_score: Self::BestScore, ptrs: Vec<Vec<PhmmStateArray<Self::Ptr>>>) -> Alignment<T>;

    /// Obtains the best scoring alignment along with its score via the Viterbi
    /// algorithm.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found, an error is given. An
    /// error is also returned if the model has no layers.
    fn viterbi(self) -> Result<Alignment<T>, PhmmError> {
        use crate::alignment::phmm::PhmmState::*;

        self.initial_check()?;

        let mut v_m = self.initialize_vm();
        let mut v_i = vec![T::INFINITY; self.query().len() + 1];
        let mut v_d = vec![T::INFINITY; self.query().len() + 1];

        let mut ptrs = self.initialize_ptrs();

        let mut i;
        let mut j = 0;

        let [layers @ .., end] = self.core().0.as_slice() else {
            return Err(PhmmError::EmptyModel);
        };

        let mut best_score = Self::BestScore::default();

        for layer in layers {
            let mut cur_m = v_m[0];

            i = 0;
            for x_idx in self.query().iter().map(|x| self.mapping().to_index(*x)) {
                let cur_vals = PhmmStateArray::new([v_d[i], cur_m, v_i[i]]);
                best_score.update(&self, cur_vals, DpIndex(i), DpIndex(j));

                let (ptr_m, match_score) = self.update_match_score(layer, x_idx, cur_vals, DpIndex(i), DpIndex(j));
                ptrs[j + 1][i + 1][Match] = ptr_m;
                cur_m = std::mem::replace(&mut v_m[i + 1], match_score);

                let (ptr_i, insert_score) = update_insert(layer, x_idx, cur_vals);
                ptrs[j][i + 1][Insert] = ptr_i.into();
                v_i[i + 1] = insert_score;

                let (ptr_d, delete_score) = update_delete(layer, cur_vals);
                ptrs[j + 1][i][Delete] = ptr_d.into();
                v_d[i] = delete_score;

                i += 1;
            }

            let cur_vals = PhmmStateArray::new([v_d[i], cur_m, v_i[i]]);
            best_score.update_seq_end(&self, cur_vals, DpIndex(j));

            let (ptr_d, delete_score) = update_delete(layer, cur_vals);
            ptrs[j + 1][i][Delete] = ptr_d.into();
            v_d[i] = delete_score;

            j += 1;
            // No version of alignment can enter a match state after BEGIN
            // without consuming a base
            v_m[0] = T::INFINITY;
        }

        // Calculate last row of insertion table and END state ptrs
        i = 0;
        for x_idx in self.query().iter().map(|x| self.mapping().to_index(*x)) {
            let cur_vals = PhmmStateArray::new([v_d[i], v_m[i], v_i[i]]);
            best_score.update_last_layer(&self, end, cur_vals, DpIndex(i));

            let (ptr_i, insert_score) = update_insert(end, x_idx, cur_vals);
            ptrs[j][i + 1][Insert] = ptr_i.into();
            v_i[i + 1] = insert_score;
            i += 1;
        }

        let cur_vals = PhmmStateArray::new([v_d[i], v_m[i], v_i[i]]);
        best_score.update_seq_end_last_layer(&self, end, cur_vals);

        // This is a necessary check, otherwise the traceback may panic
        if best_score.score() == T::INFINITY {
            return Err(PhmmError::NoAlignmentFound);
        }

        Ok(self.traceback(best_score, ptrs))
    }
}

/// Specifications for running a global Viterbi alignment.
pub struct GlobalViterbiSpecs<'a, T, const S: usize> {
    phmm:  &'a GlobalPhmm<T, S>,
    query: &'a [u8],
}

impl<'a, T: Float, const S: usize> GlobalViterbiSpecs<'a, T, S> {
    /// Create the Viterbi specifications for a global alignment from the pHMM
    /// and the sequence.
    #[inline]
    #[must_use]
    fn new(phmm: &'a GlobalPhmm<T, S>, query: &'a [u8]) -> Self {
        Self { phmm, query }
    }
}

impl<'a, T: Float, const S: usize> ViterbiSpecs<'a, T, S> for GlobalViterbiSpecs<'a, T, S> {
    type Ptr = PhmmState;
    type BestScore = GlobalBestScore<T>;

    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        &self.phmm.core
    }

    #[inline]
    fn mapping(&self) -> &ByteIndexMap<S> {
        self.phmm.mapping
    }

    #[inline]
    fn query(&self) -> &[u8] {
        self.query
    }

    #[inline]
    fn initial_check(&self) -> Result<(), PhmmError> {
        Ok(())
    }

    #[inline]
    fn initialize_vm(&self) -> Vec<T> {
        let mut out = vec![T::INFINITY; self.query.len() + 1];
        out[0] = T::ZERO;
        out
    }

    #[inline]
    fn initialize_ptrs(&self) -> Vec<Vec<PhmmStateArray<PhmmState>>> {
        vec![vec![PhmmStateArray::new([PhmmState::Match; 3]); self.query.len() + 1]; self.phmm.core.0.len()]
    }

    #[inline]
    fn update_match_score(
        &self, layer: &LayerParams<T, S>, x_idx: usize, cur_vals: PhmmStateArray<T>, _i: impl QueryIndex, _j: impl PhmmIndex,
    ) -> (PhmmState, T) {
        update_match(layer, x_idx, cur_vals)
    }

    fn traceback(self, best_score: GlobalBestScore<T>, ptrs: Vec<Vec<PhmmStateArray<PhmmState>>>) -> Alignment<T> {
        let mut i = self.query.len();
        let mut j = self.phmm.core.0.len() - 1;
        let GlobalBestScore { mut ptr, score } = best_score;

        let mut states = AlignmentStates::new();

        while i > 0 || j > 0 {
            states.add_state(ptr.to_op());
            let next_ptr = ptrs[j][i][ptr];
            match ptr {
                PhmmState::Match => {
                    i -= 1;
                    j -= 1;
                }
                PhmmState::Insert => {
                    i -= 1;
                }
                PhmmState::Delete => {
                    j -= 1;
                }
            }
            ptr = next_ptr;
        }

        states.make_reverse();

        Alignment {
            score,
            ref_range: 0..self.phmm.core.ref_length(),
            query_range: 0..self.query.len(),
            states,
            ref_len: self.phmm.core.ref_length(),
            query_len: self.query.len(),
        }
    }
}

/// Specifications for how to track the best score for the Viterbi algorithm.
///
/// This is designed to encapsulate behavior between the various types of
/// alignment. Four hooks are provided to allow the best score to update at
/// various points in the Viterbi algorithm.
pub(crate) trait BestScore<T, const S: usize> {
    type Specs<'a>
    where
        T: 'a;

    /// Retrieve the current best score
    fn score(&self) -> T;

    /// Update the best score before the end of the model and sequence.
    ///
    /// Provided to the function is the current number of bases consumed `i`,
    /// the current model layer `j`, and the current values for the delete,
    /// match, and insert states `vals`.
    fn update(&mut self, _specs: &Self::Specs<'_>, _vals: PhmmStateArray<T>, _i: impl QueryIndex, _j: impl PhmmIndex) {}

    /// Update the best score before the end of the model but at the end of the
    /// sequence.
    ///
    /// Provided to the function is the current model layer `j` and the current
    /// values for the delete, match, and insert states `vals`.
    fn update_seq_end(&mut self, _specs: &Self::Specs<'_>, _vals: PhmmStateArray<T>, _j: impl PhmmIndex) {}

    /// Update the best score before the end of the sequence but at the end of
    /// the model.
    ///
    /// ## Arguments
    /// * `layer`: the last layer, including transition probabilities into the
    ///   END state
    /// * `vals`: the values in the delete, match, and insert states of the
    ///   layer before the END state. This function must handle early exits out
    ///   of that state as well as exits through the END state.
    /// * `i`: The current number of bases consumed
    fn update_last_layer(
        &mut self, _specs: &Self::Specs<'_>, _layer: &LayerParams<T, S>, _vals: PhmmStateArray<T>, _i: impl QueryIndex,
    ) {
    }

    /// Update the best score at the end of the model and sequence.
    ///
    /// Provided to the function is the last layer `layer` and the current
    /// values for the delete, match, and insert states `vals`.
    fn update_seq_end_last_layer(&mut self, _specs: &Self::Specs<'_>, _layer: &LayerParams<T, S>, _vals: PhmmStateArray<T>) {
    }
}

/// A tracker for the best score for the global Viterbi algorithm.
pub(crate) struct GlobalBestScore<T> {
    /// A pointer to the state from which the END state was reached
    ptr:   PhmmState,
    /// The score in the END state
    score: T,
}

impl<T: Float, const S: usize> BestScore<T, S> for GlobalBestScore<T> {
    type Specs<'a>
        = GlobalViterbiSpecs<'a, T, S>
    where
        T: 'a;

    #[inline]
    fn score(&self) -> T {
        self.score
    }

    #[inline]
    fn update_seq_end_last_layer(
        &mut self, _specs: &Self::Specs<'_>, layer: &LayerParams<T, S>, mut vals: PhmmStateArray<T>,
    ) {
        use crate::alignment::phmm::PhmmState::*;

        vals[Delete] += layer.transition[(Delete, Match)];
        vals[Match] += layer.transition[(Match, Match)];
        vals[Insert] += layer.transition[(Insert, Match)];

        (self.ptr, self.score) = vals.locate_min();
    }
}

impl<T: Float> Default for GlobalBestScore<T> {
    #[inline]
    fn default() -> Self {
        Self {
            ptr:   PhmmState::Match,
            score: T::INFINITY,
        }
    }
}

impl<T: Float, const S: usize> GlobalPhmm<T, S> {
    /// Obtain the best scoring global alignment along with its score via the
    /// Viterbi algorithm.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found, an error is given. An
    /// error is also returned if the model has no layers.
    #[inline]
    pub fn viterbi<Q: AsRef<[u8]>>(&self, seq: Q) -> Result<Alignment<T>, PhmmError> {
        let seq = seq.as_ref();
        let specs = GlobalViterbiSpecs::new(self, seq);
        specs.viterbi()
    }
}
