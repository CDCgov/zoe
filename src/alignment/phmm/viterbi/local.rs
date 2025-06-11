use crate::{
    alignment::{
        Alignment, AlignmentStates,
        phmm::{
            Begin, BestScore, CorePhmm, DpIndex, End, LastBase, LayerParams, LocalPhmm, PhmmError, PhmmNumber, PhmmState,
            PhmmStateArray, PhmmStateOrEnter, PrecomputedLocalModule, QueryIndex, QueryIndexable, ViterbiStrategy,
            ViterbiTraceback,
            indexing::{LastMatch, PhmmIndex, PhmmIndexable},
            viterbi::{ExitLocation, update_match},
        },
    },
    data::ByteIndexMap,
};
use std::ops::Bound::{Excluded, Included};

/// Parameters for running a local Viterbi alignment, including the pHMM
/// information and the query.
pub struct LocalViterbiParams<'a, T, const S: usize> {
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
            mapping: phmm.mapping,
            core: &phmm.core,
            begin: phmm.begin.precompute_begin_mod(query, phmm.mapping),
            end: phmm.end.precompute_end_mod(query, phmm.mapping),
            query,
        }
    }
}

impl<'a, T: PhmmNumber, const S: usize> ViterbiStrategy<'a, T, S> for LocalViterbiParams<'a, T, S> {
    type TracebackState = PhmmStateOrEnter;
    type BestScore = LocalBestScore<T>;

    const TRACEBACK_DEFAULT: PhmmStateArray<PhmmStateOrEnter> = PhmmStateArray::new([PhmmStateOrEnter::Enter; 3]);

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
        &self, layer: &LayerParams<T, S>, x_idx: usize, cur_vals: PhmmStateArray<T>, i: impl QueryIndex, j: impl PhmmIndex,
    ) -> (PhmmStateOrEnter, T) {
        // We are updating the match score of the layer after j after consuming
        // one more base after i. To reach this cell via entering, one must
        // consume all bases up to i in the begin module, then the (i+1)st is
        // consumed in this match state. The emission parameter is added within
        // `update_match`.
        let enter_score = self.begin.get_score(i, j.next_index());
        update_match(layer, x_idx, cur_vals.with_enter(enter_score))
    }

    fn perform_traceback(
        self, best_score: LocalBestScore<T>, traceback: ViterbiTraceback<PhmmStateArray<PhmmStateOrEnter>>,
    ) -> Alignment<T> {
        use PhmmState::*;

        let LocalBestScore {
            score,
            i: DpIndex(end_i),
            loc,
        } = best_score;
        let (mut state, end_j) = match loc {
            ExitLocation::Begin => {
                let mut states = AlignmentStates::new();
                states.soft_clip(self.query.len());
                return Alignment {
                    score,
                    ref_range: 0..0,
                    query_range: 0..0,
                    states,
                    ref_len: self.core.ref_length(),
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
                        ref_range: self.core.get_seq_range(Included(End), Excluded(End)),
                        query_range: 0..0,
                        states,
                        ref_len: self.core.ref_length(),
                        query_len: self.query.len(),
                    };
                };
                (state, self.core.get_dp_index(LastMatch))
            }
        };
        let mut i = end_i;
        let mut j = end_j;

        let mut states = AlignmentStates::new();
        states.soft_clip(self.query.len() - i);

        // Continue traceback until BEGIN state is reached, or until Enter is
        // encountered (see break statement)
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
                PhmmStateOrEnter::Delete => Delete,
                PhmmStateOrEnter::Match => Match,
                PhmmStateOrEnter::Insert => Insert,
                PhmmStateOrEnter::Enter => {
                    break;
                }
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

/// A tracker for the best score for the local Viterbi algorithm.
pub struct LocalBestScore<T> {
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
    fn update(&mut self, strategy: &Self::Strategy<'_>, vals: PhmmStateArray<T>, i: impl QueryIndex, j: impl PhmmIndex) {
        let score = vals[PhmmState::Match] + strategy.end.get_score(i, j);
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
    fn update_seq_end(&mut self, strategy: &Self::Strategy<'_>, vals: PhmmStateArray<T>, j: impl PhmmIndex) {
        self.update(strategy, vals, LastBase, j);
    }

    #[inline]
    fn update_last_layer(
        &mut self, strategy: &Self::Strategy<'_>, layer: &LayerParams<T, S>, mut vals: PhmmStateArray<T>, i: impl QueryIndex,
    ) {
        use crate::alignment::phmm::PhmmState::*;

        // Option 1: Early exit from this layer
        self.update(strategy, vals, i, LastMatch);

        // Option 2: Go through END state
        vals[Delete] += layer.transition[(Delete, Match)];
        vals[Match] += layer.transition[(Match, Match)];
        vals[Insert] += layer.transition[(Insert, Match)];

        let (state, mut score) = vals.with_enter(strategy.begin.get_score(i, End)).locate_min();
        score += strategy.end.get_score(i, End);

        if score < self.score {
            self.score = score;
            self.i = strategy.query.to_dp_index(i);
            self.loc = ExitLocation::End(state);
        }
    }

    #[inline]
    fn update_seq_end_last_layer(
        &mut self, strategy: &Self::Strategy<'_>, layer: &LayerParams<T, S>, vals: PhmmStateArray<T>,
    ) {
        self.update_last_layer(strategy, layer, vals, LastBase);
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
