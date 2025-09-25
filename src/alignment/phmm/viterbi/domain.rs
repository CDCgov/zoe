use crate::{
    alignment::{
        Alignment, AlignmentStates,
        phmm::{
            BestScore, CorePhmm, DomainPhmm, DpIndex, LastBase, LayerParams, PhmmError, PhmmIndexable, PhmmNumber,
            PhmmState, PhmmStateArray, PrecomputedDomainModule, QueryIndex, QueryIndexable, ViterbiStrategy,
            ViterbiTraceback,
            indexing::{LastMatch, PhmmIndex},
            viterbi::update_match,
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

    const TRACEBACK_DEFAULT: PhmmStateArray<PhmmState> = PhmmStateArray::new([PhmmState::Match; 3]);

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
        &self, layer: &LayerParams<T, S>, x_idx: usize, cur_vals: PhmmStateArray<T>, _i: impl QueryIndex, _j: impl PhmmIndex,
    ) -> (PhmmState, T) {
        update_match(layer, x_idx, cur_vals)
    }

    fn perform_traceback(
        self, best_score: DomainBestScore<T>, traceback: ViterbiTraceback<PhmmStateArray<PhmmState>>,
    ) -> Alignment<T> {
        use PhmmState::*;

        let DomainBestScore {
            score,
            i: DpIndex(end_i),
            ptr: mut state,
        } = best_score;
        let end_j = self.core.get_dp_index(LastMatch);

        let mut i = end_i;
        let mut j = end_j;

        let mut states = AlignmentStates::new();
        states.soft_clip(self.query.len() - i);

        // Continue traceback until BEGIN state is reached
        while j > 0 || state != Match {
            states.add_state(state.to_op());
            let next_state = traceback.get(i, j)[state];

            match state {
                Match => {
                    i -= 1;
                    j -= 1;
                }
                Insert => {
                    i -= 1;
                }
                Delete => {
                    j -= 1;
                }
            }

            state = match next_state {
                PhmmState::Delete => Delete,
                PhmmState::Match => Match,
                PhmmState::Insert => Insert,
            };
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
    ptr:   PhmmState,
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
        &mut self, strategy: &Self::Strategy<'_>, layer: &LayerParams<T, S>, mut vals: PhmmStateArray<T>, i: impl QueryIndex,
    ) {
        use crate::alignment::phmm::PhmmState::*;

        vals[Delete] += layer.transition[(Delete, Match)];
        vals[Match] += layer.transition[(Match, Match)];
        vals[Insert] += layer.transition[(Insert, Match)];

        let (ptr, mut score) = vals.locate_min();
        score += strategy.end.get_score(i);

        if score < self.score {
            self.score = score;
            self.i = strategy.query.to_dp_index(i);
            self.ptr = ptr;
        }
    }

    #[inline]
    fn update_seq_end_last_layer(
        &mut self, strategy: &Self::Strategy<'_>, layer: &LayerParams<T, S>, vals: PhmmStateArray<T>,
    ) {
        self.update_last_layer(strategy, layer, vals, LastBase);
    }
}

impl<T: PhmmNumber> Default for DomainBestScore<T> {
    #[inline]
    fn default() -> Self {
        Self {
            score: T::INFINITY,
            i:     DpIndex(0),
            ptr:   PhmmState::Match,
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
