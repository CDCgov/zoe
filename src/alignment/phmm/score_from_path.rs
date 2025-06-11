use crate::{
    alignment::{
        AlignmentStates,
        phmm::{Begin, CorePhmm, DpIndex, End, GlobalPhmm, IndexOffset, PhmmError, PhmmState, SeqIndex},
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
    /// Get the score for a particular path.
    ///
    /// This is designed to give the exact same score as [`viterbi`] when the
    /// best `path` is passed, performing all arithmetic operations in the same
    /// order so as not to change the floating point error.
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
