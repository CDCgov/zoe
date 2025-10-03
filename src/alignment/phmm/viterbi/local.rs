use super::{BestScore, ViterbiStrategy, ViterbiTraceback};
use crate::{
    alignment::{
        Alignment, AlignmentStates,
        phmm::{
            CorePhmm, GetLayer, GetModule, LayerParams, LocalPhmm, PhmmBacktrackFlags, PhmmError, PhmmNumber, PhmmState,
            PhmmStateOrEnter, PhmmTracebackState, best_state_or_enter,
            indexing::{Begin, DpIndex, End, LastBase, LastMatch, PhmmIndex, PhmmIndexable, QueryIndex, QueryIndexable},
            modules::PrecomputedLocalModule,
            viterbi::ExitLocation,
        },
    },
    data::ByteIndexMap,
};
use std::ops::Bound::{Excluded, Included};

/// Parameters for running a local Viterbi alignment, including the pHMM
/// information and the query.
struct LocalViterbiParams<'a, T, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    core:    &'a CorePhmm<T, S>,
    begin:   PrecomputedLocalModule<'a, T, S>,
    end:     PrecomputedLocalModule<'a, T, S>,
    query:   &'a [u8],
}

impl<'a, T: PhmmNumber, const S: usize> LocalViterbiParams<'a, T, S> {
    /// Groups the parameters for a local Viterbi alignment in
    /// [`LocalViterbiParams`].
    #[inline]
    #[must_use]
    fn new(phmm: &'a LocalPhmm<T, S>, query: &'a [u8]) -> Self {
        Self {
            mapping: phmm.mapping(),
            core: phmm.core(),
            begin: phmm.begin().precompute_begin_mod(query, phmm.mapping()),
            end: phmm.end().precompute_end_mod(query, phmm.mapping()),
            query,
        }
    }
}

impl<'a, T: PhmmNumber, const S: usize> ViterbiStrategy<'a, T, S> for LocalViterbiParams<'a, T, S> {
    type TracebackState = PhmmStateOrEnter;
    type BestScore = LocalBestScore<T>;

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
        if self.begin.num_pseudomatch() != self.core.num_pseudomatch()
            || self.end.num_pseudomatch() != self.core.num_pseudomatch()
            || self.begin.internal_params.seq_len() != self.query.seq_len()
            || self.end.internal_params.seq_len() != self.query.seq_len()
        {
            return Err(PhmmError::IncompatibleModule);
        }
        Ok(())
    }

    #[inline]
    fn fill_vm(&self, v_m: &mut [T]) {
        for (i, value) in v_m.iter_mut().enumerate() {
            *value = self.begin.get_score(DpIndex(i), Begin);
        }
    }

    fn update_match_score(
        &self, layer: &LayerParams<T, S>, x_idx: usize, mut match_val: T, mut delete_val: T, mut insert_val: T,
        i: impl QueryIndex, j: impl PhmmIndex,
    ) -> (Self::TracebackState, T) {
        use crate::alignment::phmm::state::PhmmState::*;

        // We are updating the match score of the layer after j after consuming
        // one more base after i. To reach this cell via entering, one must
        // consume all bases up to i in the begin module, then the (i+1)st is
        // consumed in this match state. The emission parameter is added within
        // `update_match`.
        let enter_val = self.begin.get_score(i, j.next_index());

        match_val += layer.transition[(Match, Match)];
        delete_val += layer.transition[(Delete, Match)];
        insert_val += layer.transition[(Insert, Match)];

        let (state, best) = best_state_or_enter(match_val, delete_val, insert_val, enter_val);

        (state, best + layer.emission_match[x_idx])
    }

    fn perform_traceback(
        self, best_score: LocalBestScore<T>, traceback: ViterbiTraceback<PhmmBacktrackFlags>,
    ) -> Alignment<T> {
        use PhmmState::*;

        let LocalBestScore {
            score,
            i: DpIndex(end_i),
            loc,
        } = best_score;
        let (state, end_j) = match loc {
            ExitLocation::Begin => {
                let mut states = AlignmentStates::new();
                states.soft_clip(self.query.len());
                return Alignment {
                    score,
                    ref_range: 0..0,
                    query_range: 0..0,
                    states,
                    ref_len: self.core.seq_len(),
                    query_len: self.query.len(),
                };
            }
            ExitLocation::Match(j) => (Match, self.core.get_dp_index(j)),
            ExitLocation::End(state) => {
                let Some(state) = PhmmState::get_from(state) else {
                    let mut states = AlignmentStates::new();
                    states.soft_clip(self.query.len());
                    return Alignment {
                        score,
                        ref_range: self.core.get_seq_range(End..End),
                        query_range: 0..0,
                        states,
                        ref_len: self.core.seq_len(),
                        query_len: self.query.len(),
                    };
                };
                (state, self.core.get_dp_index(LastMatch))
            }
        };

        let mut state = PhmmTracebackState::from(state);
        let mut i = end_i;
        let mut j = end_j;

        let mut states = AlignmentStates::new();
        states.soft_clip(self.query.len() - i);

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

        Alignment {
            score,
            ref_range: self.core.get_seq_range((start_j, end_j)),
            query_range: self.query.get_seq_range(start_i, end_i),
            states,
            ref_len: self.core.seq_len(),
            query_len: self.query.len(),
        }
    }
}

/// A tracker for the best score for the local Viterbi algorithm.
struct LocalBestScore<T> {
    /// The best score found at the end of the model
    score: T,
    /// The last query index consumed by the [`CorePhmm`]
    i:     DpIndex,
    /// The location from which the alignment exits the [`CorePhmm`]
    loc:   ExitLocation,
}

impl<T: PhmmNumber, const S: usize> BestScore<T, S> for LocalBestScore<T> {
    type Strategy<'a>
        = LocalViterbiParams<'a, T, S>
    where
        T: 'a;

    #[inline]
    fn score(&self) -> T {
        self.score
    }

    #[inline]
    fn update(
        &mut self, strategy: &Self::Strategy<'_>, match_val: T, _delete_val: T, _insert_val: T, i: impl QueryIndex,
        j: impl PhmmIndex,
    ) {
        let score = match_val + strategy.end.get_score(i, j);
        if score < self.score {
            self.score = score;
            self.i = strategy.query.to_dp_index(i);
            self.loc = match strategy.core.to_seq_index(j) {
                Some(loc) => ExitLocation::Match(loc),
                None => ExitLocation::Begin,
            }
        }
    }

    #[inline]
    fn update_seq_end(
        &mut self, strategy: &Self::Strategy<'_>, match_val: T, delete_val: T, insert_val: T, j: impl PhmmIndex,
    ) {
        self.update(strategy, match_val, delete_val, insert_val, LastBase, j);
    }

    fn update_last_layer(
        &mut self, strategy: &Self::Strategy<'_>, layer: &LayerParams<T, S>, mut match_val: T, mut delete_val: T,
        mut insert_val: T, i: impl QueryIndex,
    ) {
        use crate::alignment::phmm::PhmmState::*;

        // Option 1: Early exit from this layer
        self.update(strategy, match_val, delete_val, insert_val, i, LastMatch);

        // Option 2: Go through END state
        match_val += layer.transition[(Match, Match)];
        delete_val += layer.transition[(Delete, Match)];
        insert_val += layer.transition[(Insert, Match)];
        let enter_val = strategy.begin.get_score(i, End);

        let (state, mut score) = best_state_or_enter(match_val, delete_val, insert_val, enter_val);
        score += strategy.end.get_score(i, End);

        if score < self.score {
            self.score = score;
            self.i = strategy.query.to_dp_index(i);
            self.loc = ExitLocation::End(state);
        }
    }

    #[inline]
    fn update_seq_end_last_layer(
        &mut self, strategy: &Self::Strategy<'_>, layer: &LayerParams<T, S>, match_val: T, delete_val: T, insert_val: T,
    ) {
        self.update_last_layer(strategy, layer, match_val, delete_val, insert_val, LastBase);
    }
}

impl<T: PhmmNumber> Default for LocalBestScore<T> {
    #[inline]
    fn default() -> Self {
        Self {
            score: T::INFINITY,
            i:     DpIndex(0),
            loc:   ExitLocation::End(PhmmStateOrEnter::Match),
        }
    }
}

impl<T: PhmmNumber, const S: usize> LocalPhmm<T, S> {
    /// Obtain the best scoring local alignment along with its score via the
    /// Viterbi algorithm.
    ///
    /// ## Errors
    ///
    /// If no alignment with nonzero probability is found, an error is given. An
    /// error is also returned if the model has no layers.
    pub fn viterbi<Q: AsRef<[u8]>>(&self, seq: Q) -> Result<Alignment<T>, PhmmError> {
        let seq = seq.as_ref();
        let specs = LocalViterbiParams::new(self, seq);
        specs.viterbi()
    }
}
