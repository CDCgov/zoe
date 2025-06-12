use crate::{
    alignment::{
        Alignment, AlignmentStates,
        phmm::{
            BestScore, CorePhmm, GlobalPhmm, LayerParams, PhmmError, PhmmNumber, PhmmState, PhmmStateArray, QueryIndex,
            ViterbiStrategy, ViterbiTraceback, indexing::PhmmIndex, viterbi::update_match,
        },
    },
    data::ByteIndexMap,
};

/// Parameters for running a global Viterbi alignment, including the pHMM
/// information and the query.
pub struct GlobalViterbiParams<'a, T, const S: usize> {
    phmm:  &'a GlobalPhmm<T, S>,
    query: &'a [u8],
}

impl<'a, T: PhmmNumber, const S: usize> GlobalViterbiParams<'a, T, S> {
    /// Groups the parameters for a global Viterbi alignment in
    /// [`GlobalViterbiParams`].
    #[inline]
    #[must_use]
    fn new(phmm: &'a GlobalPhmm<T, S>, query: &'a [u8]) -> Self {
        Self { phmm, query }
    }
}

impl<'a, T: PhmmNumber, const S: usize> ViterbiStrategy<'a, T, S> for GlobalViterbiParams<'a, T, S> {
    type TracebackState = PhmmState;
    type BestScore = GlobalBestScore<T>;

    const TRACEBACK_DEFAULT: PhmmStateArray<PhmmState> = PhmmStateArray::new([PhmmState::Match; 3]);

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
    fn fill_vm(&self, v_m: &mut [T]) {
        v_m.fill(T::INFINITY);
        v_m[0] = T::ZERO;
    }

    #[inline]
    fn update_match_score(
        &self, layer: &LayerParams<T, S>, x_idx: usize, cur_vals: PhmmStateArray<T>, _i: impl QueryIndex, _j: impl PhmmIndex,
    ) -> (PhmmState, T) {
        update_match(layer, x_idx, cur_vals)
    }

    fn perform_traceback(
        self, best_score: GlobalBestScore<T>, traceback: ViterbiTraceback<PhmmStateArray<PhmmState>>,
    ) -> Alignment<T> {
        // i = len(query); j = len(layers);
        let mut cursor = traceback.data.len() - 1;
        let GlobalBestScore { mut state, score } = best_score;

        let mut states = AlignmentStates::new();

        while cursor > 0 {
            states.add_state(state.to_op());
            let next_state = traceback.data[cursor][state];
            match state {
                PhmmState::Match => {
                    // i -= 1; j -= 1;
                    cursor -= traceback.cols + 1;
                }
                PhmmState::Insert => {
                    // i -= 1;
                    cursor -= 1;
                }
                PhmmState::Delete => {
                    // j -= 1;
                    cursor -= traceback.cols;
                }
            }
            state = next_state;
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

/// A tracker for the best score for the global Viterbi algorithm.
pub(crate) struct GlobalBestScore<T> {
    /// The state from which the END state was reached
    state: PhmmState,
    /// The score in the END state
    score: T,
}

impl<T: PhmmNumber, const S: usize> BestScore<T, S> for GlobalBestScore<T> {
    type Specs<'a>
        = GlobalViterbiParams<'a, T, S>
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

        (self.state, self.score) = vals.locate_min();
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
        let specs = GlobalViterbiParams::new(self, seq);
        specs.viterbi()
    }
}
