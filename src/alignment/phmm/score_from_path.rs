use crate::{
    alignment::{
        AlignmentStates, StatesSequence,
        phmm::{
            Begin, CorePhmm, DomainPhmm, DpIndex, End, FirstMatch, GlobalPhmm, IndexOffset, LocalPhmm, PhmmError, PhmmState,
            SemiLocalPhmm, SeqIndex, indexing::LastMatch,
        },
    },
    data::{ByteIndexMap, cigar::Ciglet},
    math::Float,
};
use std::ops::Range;

/// Given a PHMM alignment specified by `ciglets`, get the score for the portion
/// corresponding to the [`CorePhmm`].
///
/// `ref_range` specifies the range of the PHMM's underlying "reference" that
/// should be consumed while scoring. For global alignment, this is the full
/// length of the PHMM (as can be found with [`ref_length`]). For other forms of
/// alignment, this may be a smaller range.
///
/// `seq_in_alignment` should contain exactly the portion of the sequence
/// matched by the model in the given alignment.
///
/// `score` is the starting score to tally into (provided since order of
/// operations may matter with floating point arithmetic).
///
/// This function does not include the scores caused by:
/// * Transitions from the BEGIN state
/// * Transitions into the END state
/// * Transitions from a starting or ending module into a match state
///
/// [`ref_length`]: CorePhmm::ref_length
fn score_from_path_core<T: Float, const S: usize>(
    core: &CorePhmm<T, S>, mapping: &ByteIndexMap<S>, seq_in_alignment: &[u8], ref_range: Range<usize>,
    ciglets: impl IntoIterator<Item = Ciglet>, score: T,
) -> Result<(T, PhmmState), PhmmError> {
    use PhmmState::*;

    // TODO: Improve ergonomics later...
    let mut op_iter = ciglets
        .into_iter()
        .flat_map(|ciglet| std::iter::repeat_n(ciglet.op, ciglet.inc));

    let Some(op) = op_iter.next() else {
        if ref_range.is_empty() {
            // If the CIGAR string is empty, then either the path entered the
            // BEGIN state and then exited, or entered the END state and then
            // exited. In either case, a match state was passed through.
            return Ok((T::ZERO, Match));
        }
        return Err(PhmmError::FullModelNotUsed);
    };

    // Index into the query
    let mut i = 0;

    // Initialize current score, starting state, and starting layer index
    let (mut score, mut state, mut j) = match PhmmState::from_op(op)? {
        Delete => {
            // Cannot enter into a delete state except from BEGIN.
            if ref_range.start != 0 {
                return Err(PhmmError::InvalidPath);
            }
            // After entering the delete state, we are now in layer 1
            (score, Delete, 1)
        }
        Match => {
            let x_idx = mapping.to_index(seq_in_alignment[i]);
            // Need to look at the previous layer to get transitions into this
            // layer
            let layer_idx = SeqIndex(ref_range.start).prev_index();
            let score = score + core.get_layer(layer_idx).emission_match[x_idx];
            i += 1;
            (score, Match, ref_range.start + 1)
        }
        Insert => {
            // Cannot enter into a insert state except from BEGIN.
            if ref_range.start != 0 {
                return Err(PhmmError::InvalidPath);
            }
            // Insertion before reference position n has emission parameters and
            // transition parameters stored in layer n
            let x_idx = mapping.to_index(seq_in_alignment[i]);
            let score = score + core.get_layer(Begin).emission_insert[x_idx];
            i += 1;
            // After entering the delete state, we are still in begin layer
            (score, Insert, 0)
        }
    };

    // Main loop: first score transition out of current layer index j, then
    // score emission in new state if applicable
    for op in op_iter {
        let layer = core.get_layer(DpIndex(j));

        match PhmmState::from_op(op)? {
            Match => {
                let x_idx = mapping.to_index(seq_in_alignment[i]);
                // Must perform the two additions separately, since floating
                // point addition is not associative
                score += layer.transition[(state, Match)];
                score += layer.emission_match[x_idx];
                state = Match;
                i += 1;
                j += 1;
            }
            Insert => {
                let x_idx = mapping.to_index(seq_in_alignment[i]);
                // Must perform the two additions separately, since floating
                // point addition is not associative
                score += layer.transition[(state, Insert)];
                score += layer.emission_insert[x_idx];
                state = Insert;
                i += 1;
            }
            Delete => {
                score += layer.transition[(state, Delete)];
                state = Delete;
                j += 1;
            }
        }
    }

    // i and j were incremented 1 past their last value
    if i != seq_in_alignment.len() {
        return Err(PhmmError::FullSeqNotUsed);
    }
    if j != ref_range.end {
        return Err(PhmmError::FullModelNotUsed);
    }

    Ok((score, state))
}

impl<T: Float, const S: usize> GlobalPhmm<T, S> {
    /// Get the score for a particular alignment.
    ///
    /// This is designed to give the exact same score as [`viterbi`] when the
    /// best `alignment` is passed, performing all arithmetic operations in the
    /// same order so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are M, =, X, I, and D.
    ///
    /// [`viterbi`]: GlobalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(&self, seq: Q, alignment: &AlignmentStates) -> Result<T, PhmmError> {
        use PhmmState::*;

        let seq = seq.as_ref();

        let first_op = alignment.iter().next().ok_or(PhmmError::FullModelNotUsed)?.op;
        let first_state = PhmmState::from_op(first_op)?;

        // Get transition from BEGIN state the actual first state
        let score = self.core.get_layer(Begin).transition[(Match, first_state)];

        // Add score from first state after START until final state before END
        let (mut score, final_state) =
            score_from_path_core(&self.core, self.mapping, seq, 0..self.core.ref_length(), alignment, score)?;

        // Add transition into END state
        score += self.core.get_layer(End).transition[(final_state, Match)];

        Ok(score)
    }
}

impl<T: Float, const S: usize> LocalPhmm<T, S> {
    /// Get the best score for an empty alignment (all soft clipping).
    fn score_empty_alignment<Q: AsRef<[u8]>>(&self, seq: Q) -> T {
        let seq = seq.as_ref();
        let mut best_score = T::INFINITY;

        for i in 0..=seq.len() {
            let (inserted_begin, inserted_end) = seq.split_at(i);
            let score_through_begin = {
                let begin_score = self.begin.internal_params.get_begin_score(inserted_begin, self.mapping)
                    + self.get_begin_external_score(Begin);
                let end_score = self.get_end_internal_score(inserted_end, self.mapping) + self.get_end_external_score(Begin);
                begin_score + end_score
            };
            let score_through_end = {
                let begin_score = self.begin.internal_params.get_begin_score(inserted_begin, self.mapping)
                    + self.get_begin_external_score(End);
                let end_score = self.get_end_internal_score(inserted_end, self.mapping) + self.get_end_external_score(End);
                begin_score + end_score
            };
            let score = score_through_begin.min(score_through_end);
            if score < best_score {
                best_score = score;
            }
        }

        best_score
    }

    /// Get the best score for a particular alignment.
    ///
    /// The score may be infinite if no paths have a nonzero probability. This
    /// is designed to give the exact same score as [`viterbi`] when the best
    /// `path` is passed, performing all arithmetic operations in the same order
    /// so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are M, =, X, I, and D. If `path` does not
    /// correspond to a valid path through a [`LocalPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    ///
    /// [`viterbi`]: LocalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(
        &self, seq: Q, alignment: &AlignmentStates, ref_range: Range<usize>,
    ) -> Result<T, PhmmError> {
        use PhmmState::*;

        let seq = seq.as_ref();
        let mut ciglets = alignment.as_slice();
        let mut score = T::ZERO;

        // Split query between begin module, core pHMM, and end module
        let (begin_seq, seq, end_seq) = {
            let begin_inserted = ciglets.next_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
            let end_inserted = ciglets.next_back_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
            if ciglets.is_empty() {
                return Ok(self.score_empty_alignment(seq));
            }

            let (seq, end_seq) = seq.split_at(seq.len() - end_inserted);
            let (begin_seq, seq) = seq.split_at(begin_inserted);
            (begin_seq, seq, end_seq)
        };

        // Add contribution from internal parameters of begin module
        score += self.get_begin_internal_score(begin_seq, self.mapping);

        // Add the contribution of the external parameters into the core pHMM
        // and the transitions out of the BEGIN state
        if ref_range.start == 0 {
            let first_op = ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?;
            let first_state = PhmmState::from_op(first_op)?;
            let score_through_begin =
                score + self.get_begin_external_score(Begin) + self.core.get_layer(Begin).transition[(Match, first_state)];

            match first_state {
                Delete | Insert => score = score_through_begin,
                Match => {
                    // Handle ambiguity: it is unclear whether we pass through
                    // the BEGIN state of the core PHMM or not
                    let score_skipping_begin = score + self.get_begin_external_score(FirstMatch);
                    score = score_through_begin.min(score_skipping_begin);
                }
            }
        } else {
            // We must enter into a match state if not going through BEGIN
            if PhmmState::from_op(ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?)? != Match {
                return Err(PhmmError::InvalidPath);
            }

            score += self.get_begin_external_score(SeqIndex(ref_range.start));
        }

        // Add the contribution from the core pHMM excluding transitions into
        // the END state, and obtain the final state before the END state
        let (mut score, final_state) = score_from_path_core(
            &self.core,
            self.mapping,
            seq,
            ref_range.clone(),
            ciglets.iter().copied(),
            score,
        )?;

        // Compute the contribution from the internal parameters of the end
        // module, which involves processing any soft clipping at the end of the
        // alignment
        let end_score_internal = self.get_end_internal_score(end_seq, self.mapping);

        score = if ref_range.end == self.core.ref_length() {
            let score_through_end = score
                + self.core.get_layer(End).transition[(final_state, Match)]
                + (end_score_internal + self.get_end_external_score(End));

            match final_state {
                Delete | Insert => score_through_end,
                Match => {
                    // Handle ambiguity: it is unclear whether we pass through
                    // the END state of the core PHMM or not
                    let score_skipping_end = score + (end_score_internal + self.get_end_external_score(LastMatch));
                    score_through_end.min(score_skipping_end)
                }
            }
        } else {
            // We must exit from a match state if not going through END
            if final_state != Match {
                return Err(PhmmError::InvalidPath);
            }

            // Subtract 1 since range is end-exclusive
            score + (end_score_internal + self.get_end_external_score(SeqIndex(ref_range.end - 1)))
        };

        Ok(score)
    }
}

impl<T: Float, const S: usize> DomainPhmm<T, S> {
    // TODO: Doc link
    /// Get the best score for a particular alignment.
    ///
    /// The score may be infinite if no paths have a nonzero probability. This
    /// is designed to give the exact same score as `viterbi` when the best
    /// `path` is passed, performing all arithmetic operations in the same order
    /// so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are M, =, X, I, and D. If `path` does not
    /// correspond to a valid path through a [`DomainPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    pub fn score_from_path<Q: AsRef<[u8]>>(&self, seq: Q, alignment: &AlignmentStates) -> Result<T, PhmmError> {
        use PhmmState::*;

        let seq = seq.as_ref();
        let mut ciglets = alignment.as_slice();
        let mut score = T::ZERO;

        // Split query between begin module, core pHMM, and end module
        let (begin_seq, seq, end_seq) = {
            let begin_inserted = ciglets.next_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
            let end_inserted = ciglets.next_back_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);

            let (seq, end_seq) = seq.split_at(seq.len() - end_inserted);
            let (begin_seq, seq) = seq.split_at(begin_inserted);
            (begin_seq, seq, end_seq)
        };

        // Add contribution from parameters of begin module
        score += self.get_begin_score(begin_seq, self.mapping);

        // Add the transitions out of the BEGIN state
        let first_op = ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?;
        let first_state = PhmmState::from_op(first_op)?;
        score += self.core.get_layer(Begin).transition[(Match, first_state)];

        // Add the contribution from the core pHMM excluding transitions into
        // the END state, and obtain the final state before the END state
        let (mut score, final_state) = score_from_path_core(
            &self.core,
            self.mapping,
            seq,
            0..self.core.ref_length(),
            ciglets.iter().copied(),
            score,
        )?;

        // Add the transition into the END state
        score += self.core.get_layer(End).transition[(final_state, Match)];

        // Add the contribution from parameters of the end module
        score += self.get_end_score(end_seq, self.mapping);

        Ok(score)
    }
}

impl<T: Float, const S: usize> SemiLocalPhmm<T, S> {
    /// Get the best score for an empty alignment (which means the sequence is
    /// also empty)
    fn score_empty_alignment(&self) -> T {
        let score_through_begin = self.get_begin_score(Begin) + self.get_end_score(Begin);
        let score_through_end = self.get_begin_score(End) + self.get_end_score(End);

        score_through_begin.min(score_through_end)
    }

    /// Get the best score for a particular path. The score may be infinite if
    /// no paths have a nonzero probability.
    ///
    /// This is designed to give the exact same score as [`viterbi`] when the
    /// best `path` is passed, performing all arithmetic operations in the same
    /// order so as not to change the floating point error.
    ///
    /// ## Errors
    ///
    /// The CIGAR string must consume the entire model and sequence, and the
    /// only supported operations are M, =, X, I, and D. If `path` does not
    /// correspond to a valid path through a [`LocalPhmm`], then
    /// [`PhmmError::InvalidPath`] is returned.
    ///
    /// [`viterbi`]: LocalPhmm::viterbi
    pub fn score_from_path<Q: AsRef<[u8]>>(
        &self, seq: Q, alignment: &AlignmentStates, ref_range: Range<usize>,
    ) -> Result<T, PhmmError> {
        use PhmmState::*;

        let seq = seq.as_ref();
        let ciglets = alignment.as_slice();
        let mut score = T::ZERO;

        if ciglets.peek_op().is_none() {
            return Ok(self.score_empty_alignment());
        }

        // Add the contribution of the external parameters into the core pHMM
        // and the transitions out of the BEGIN state
        if ref_range.start == 0 {
            let first_op = ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?;
            let first_state = PhmmState::from_op(first_op)?;
            let score_through_begin =
                score + self.get_begin_score(Begin) + self.core.get_layer(Begin).transition[(Match, first_state)];

            match first_state {
                Delete | Insert => score = score_through_begin,
                Match => {
                    // Handle ambiguity: it is unclear whether we pass through
                    // the BEGIN state of the core PHMM or not
                    let score_skipping_begin = score + self.get_begin_score(FirstMatch);
                    score = score_through_begin.min(score_skipping_begin);
                }
            }
        } else {
            // We must enter into a match state if not going through BEGIN
            if PhmmState::from_op(ciglets.peek_op().ok_or(PhmmError::FullModelNotUsed)?)? != Match {
                return Err(PhmmError::InvalidPath);
            }

            score += self.get_begin_score(SeqIndex(ref_range.start));
        }

        // Add the contribution from the core pHMM excluding transitions into
        // the END state, and obtain the final state before the END state
        let (mut score, final_state) = score_from_path_core(
            &self.core,
            self.mapping,
            seq,
            ref_range.clone(),
            ciglets.iter().copied(),
            score,
        )?;

        score = if ref_range.end == self.core.ref_length() {
            let score_through_end =
                score + self.core.get_layer(End).transition[(final_state, Match)] + self.get_end_score(End);

            match final_state {
                Delete | Insert => score_through_end,
                Match => {
                    // Handle ambiguity: it is unclear whether we pass through
                    // the END state of the core PHMM or not
                    let score_skipping_end = score + self.get_end_score(LastMatch);
                    score_through_end.min(score_skipping_end)
                }
            }
        } else {
            // We must exit from a match state if not going through END
            if final_state != Match {
                return Err(PhmmError::InvalidPath);
            }

            // Subtract 1 since range is end-exclusive
            score + self.get_end_score(SeqIndex(ref_range.end - 1))
        };

        Ok(score)
    }
}
