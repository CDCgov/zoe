use crate::{
    alignment::{
        Alignment, AlignmentStates,
        phmm::{
            BestScore, CorePhmm, DomainPhmm, LayerParams, PhmmBacktrackFlags, PhmmError, PhmmNumber, PhmmState,
            PhmmTracebackState, PrecomputedDomainModule, ViterbiStrategy, ViterbiTraceback, best_state,
            indexing::{DpIndex, LastBase, LastMatch, PhmmIndex, PhmmIndexable, QueryIndex, QueryIndexable},
        },
    },
    data::ByteIndexMap,
};
use std::ops::Bound::{Excluded, Included};

/// Parameters for running a domain Viterbi alignment, including the pHMM
/// information and the query.
pub struct DomainViterbiParams<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a CorePhmm<T, S>,
    begin:   PrecomputedDomainModule<T, S>,
    end:     PrecomputedDomainModule<T, S>,
    query:   &'a [u8],
}

impl<'a, T: PhmmNumber, const S: usize> DomainViterbiParams<'a, T, S> {
    /// Create the Viterbi specifications for a local alignment from the pHMM
    /// and the sequence.
    #[inline]
    #[must_use]
    fn new(phmm: &'a DomainPhmm<T, S>, query: &'a [u8]) -> Self {
        Self {
            mapping: phmm.mapping,
            core: &phmm.core,
            begin: phmm.begin.precompute_begin_mod(query, phmm.mapping),
            end: phmm.end.precompute_end_mod(query, phmm.mapping),
            query,
        }
    }
}

impl<'a, T: PhmmNumber, const S: usize> ViterbiStrategy<'a, T, S> for DomainViterbiParams<'a, T, S> {
    type TracebackState = PhmmState;
    type BestScore = DomainBestScore<T>;

    #[inline]
    fn core(&self) -> &CorePhmm<T, S> {
        self.core
    }

    #[inline]
    fn mapping(&self) -> &ByteIndexMap<S> {
        self.mapping
    }

    #[inline]
    fn query(&self) -> &[u8] {
        self.query
    }

    #[inline]
    fn initial_check(&self) -> Result<(), PhmmError> {
        if self.begin.seq_len() != self.query.seq_len() || self.end.seq_len() != self.query.seq_len() {
            return Err(PhmmError::IncompatibleModule);
        }
        Ok(())
    }

    #[inline]
    fn fill_vm(&self, v_m: &mut [T]) {
        for (i, value) in v_m.iter_mut().enumerate() {
            *value = self.begin.get_score(DpIndex(i));
        }
    }

    fn update_match_score(
        &self, layer: &LayerParams<T, S>, x_idx: usize, mut match_val: T, mut delete_val: T, mut insert_val: T,
        _i: impl QueryIndex, _j: impl PhmmIndex,
    ) -> (PhmmState, T) {
        use crate::alignment::phmm::state::PhmmState::*;

        match_val += layer.transition[(Match, Match)];
        delete_val += layer.transition[(Delete, Match)];
        insert_val += layer.transition[(Insert, Match)];

        let (state, best) = best_state(match_val, delete_val, insert_val);
        (state, best + layer.emission_match[x_idx])
    }

    fn perform_traceback(
        self, best_score: DomainBestScore<T>, traceback: ViterbiTraceback<PhmmBacktrackFlags>,
    ) -> Alignment<T> {
        let DomainBestScore {
            score,
            i: DpIndex(end_i),
            state,
        } = best_score;
        let end_j = self.core.get_dp_index(LastMatch);

        let mut state = PhmmTracebackState::from(state);
        let mut i = end_i;
        let mut j = end_j;

        let mut states = AlignmentStates::new();
        states.soft_clip(self.query.len() - i);

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

        Alignment {
            score,
            ref_range: self.core.get_seq_range(start_j, end_j),
            query_range: self.query.get_seq_range(start_i, end_i),
            states,
            ref_len: self.core.ref_length(),
            query_len: self.query.len(),
        }
    }
}

/// A tracker for the best score for the domain Viterbi algorithm.
pub struct DomainBestScore<T> {
    /// The best score found at the end of the model
    score: T,
    /// The last query index consumed by the [`CorePhmm`]
    i:     DpIndex,
    /// The state from which the END state was reached
    state: PhmmState,
}

impl<T: PhmmNumber, const S: usize> BestScore<T, S> for DomainBestScore<T> {
    type Strategy<'a>
        = DomainViterbiParams<'a, T, S>
    where
        T: 'a;

    #[inline]
    fn score(&self) -> T {
        self.score
    }

    #[inline]
    fn update_last_layer(
        &mut self, strategy: &Self::Strategy<'_>, layer: &LayerParams<T, S>, mut match_val: T, mut delete_val: T,
        mut insert_val: T, i: impl QueryIndex,
    ) {
        use crate::alignment::phmm::PhmmState::*;

        match_val += layer.transition[(Match, Match)];
        delete_val += layer.transition[(Delete, Match)];
        insert_val += layer.transition[(Insert, Match)];

        let (state, mut score) = best_state(match_val, delete_val, insert_val);
        score += strategy.end.get_score(i);

        if score < self.score {
            self.score = score;
            self.i = strategy.query.to_dp_index(i);
            self.state = state;
        }
    }

    #[inline]
    fn update_seq_end_last_layer(
        &mut self, strategy: &Self::Strategy<'_>, layer: &LayerParams<T, S>, match_val: T, delete_val: T, insert_val: T,
    ) {
        self.update_last_layer(strategy, layer, match_val, delete_val, insert_val, LastBase);
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
    /// Obtain the best scoring local alignment along with its score via the
    /// Viterbi algorithm.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found, an error is given. An
    /// error is also returned if the model has no layers.
    pub fn viterbi<Q: AsRef<[u8]>>(&self, seq: Q) -> Result<Alignment<T>, PhmmError> {
        let seq = seq.as_ref();
        let specs = DomainViterbiParams::new(self, seq);
        specs.viterbi()
    }
}
