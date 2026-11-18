use crate::{
    alignment::{AlignOp, Alignment, ScoreParams},
    data::cigar::Ciglet,
};
use std::{
    cmp::Ordering::{Equal, Greater, Less},
    error::Error,
    fmt::{Debug, Display},
};

impl<T> Alignment<T> {
    /// Re-computes the score for a local alignment. The score field of the
    /// alignment is not utilized (and could be any type, such as `()`).
    ///
    /// The full reference and query corresponding to the alignment must be
    /// passed.
    ///
    /// ## Errors
    ///
    /// See [`RescoreSwError`].
    ///
    /// ## Panics
    ///
    /// See [here](Alignment#validity-and-panics).
    ///
    /// ## Validity
    ///
    /// No repeating adjacent states should be present in the
    /// [`AlignmentStates`] in the [`Alignment`].
    ///
    /// [`AlignmentStates`]: crate::alignment::AlignmentStates
    pub fn rescore_sw<const S: usize>(
        &self, reference: &[u8], query: &[u8], params: &ScoreParams<i8, S>,
    ) -> Result<u32, RescoreSwError> {
        if query.len() != self.query_len {
            return Err(RescoreSwError::QueryLen {
                expected: self.query_len,
                actual:   query.len(),
            });
        }
        if reference.len() != self.ref_len {
            return Err(RescoreSwError::RefLen {
                expected: self.ref_len,
                actual:   reference.len(),
            });
        }

        let ScoreParams {
            matrix,
            gap_open,
            gap_extend,
        } = params;
        let [gap_open, gap_extend] = [*gap_open, *gap_extend].map(i32::from);
        let matrix = matrix.as_ref();

        let mut score = 0;
        let mut r = self.ref_range.start;
        let mut q = self.query_range.start;

        for Ciglet { inc, op } in &self.states {
            match op {
                AlignOp::M | AlignOp::Eq | AlignOp::X => {
                    for _ in 0..inc {
                        score += i32::from(matrix.get_weight(reference[r], query[q]));
                        q += 1;
                        r += 1;
                    }
                }
                AlignOp::I => {
                    score += gap_open + gap_extend * (inc - 1) as i32;
                    q += inc;
                }
                AlignOp::D => {
                    score += gap_open + gap_extend * (inc - 1) as i32;
                    r += inc;
                }
                AlignOp::N => r += inc,
                AlignOp::S | AlignOp::H | AlignOp::P => {}
            }
        }

        // score is i32, so this cast solely could fail due to it being negative
        match u32::try_from(score) {
            Ok(score) => Ok(score),
            Err(_) => Err(RescoreSwError::NegativeScore { score }),
        }
    }
}

/// An error type for [`Alignment::rescore_sw`].
#[derive(Eq, PartialEq)]
pub enum RescoreSwError {
    /// The length of the passed query was incorrect.
    QueryLen { expected: usize, actual: usize },
    /// The length of the passed reference was incorrect.
    RefLen { expected: usize, actual: usize },
    /// A negative score was obtained during rescoring.
    NegativeScore { score: i32 },
}

impl Display for RescoreSwError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match *self {
            RescoreSwError::QueryLen { expected, actual } => write!(
                f,
                "The alignment was for a query of length {expected}, but a query of length {actual} was passed"
            ),
            RescoreSwError::RefLen { expected, actual } => write!(
                f,
                "The alignment was for a reference of length {expected}, but a reference of length {actual} was passed"
            ),
            RescoreSwError::NegativeScore { score } => write!(f, "The alignment produced a negative score: {score}"),
        }
    }
}

impl Debug for RescoreSwError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for RescoreSwError {}

/// Computes the score for a local alignment.
///
/// When scoring an alignment, this function expects the full `query` to be
/// passed, as well as the slice of the reference which was aligned to (e.g.,
/// the indices in [`Alignment::ref_range`]).
///
/// ## Validity
///
/// - The iterator of ciglets should never contain two ciglets with the same
///   operation adjacent to each other.
///
/// ## Errors
///
/// - The `ciglets` must contain valid operations in `MIDNSHP=X`.
/// - All of `query`, `ref_in_alignment`, and `cigar` must be fully consumed.
/// - The final score should be nonnegative.
///
/// ## Panics
///
/// - All increments must be nonzero.
///
/// [`Alignment::ref_range`]: crate::alignment::types::Alignment::ref_range
#[deprecated(since = "0.0.37", note = "Use Alignment::rescore_sw instead")]
#[allow(deprecated)]
pub fn sw_score_from_path<const S: usize>(
    ciglets: impl IntoIterator<Item = Ciglet>, ref_in_alignment: &[u8], query: &[u8], params: &ScoreParams<i8, S>,
) -> Result<u32, SwScoringError> {
    let ScoreParams {
        matrix,
        gap_open,
        gap_extend,
    } = params;
    let [gap_open, gap_extend] = [*gap_open, *gap_extend].map(i32::from);
    let matrix = matrix.as_ref();

    let mut score = 0;
    let mut r = 0;
    let mut q = 0;

    for Ciglet { inc, op } in ciglets {
        match op {
            AlignOp::M | AlignOp::Eq | AlignOp::X => {
                for _ in 0..inc {
                    let Some(reference_base) = ref_in_alignment.get(r).copied() else {
                        return Err(SwScoringError::ReferenceEnded);
                    };
                    let Some(query_base) = query.get(q).copied() else {
                        return Err(SwScoringError::QueryEnded);
                    };
                    score += i32::from(matrix.get_weight(reference_base, query_base));
                    q += 1;
                    r += 1;
                }
            }
            AlignOp::I => {
                score += gap_open + gap_extend * (inc - 1) as i32;
                q += inc;
            }
            AlignOp::D => {
                score += gap_open + gap_extend * (inc - 1) as i32;
                r += inc;
            }
            AlignOp::S => q += inc,
            AlignOp::N => r += inc,
            AlignOp::H | AlignOp::P => {}
        }
    }

    match q.cmp(&query.len()) {
        Less => return Err(SwScoringError::FullQueryNotUsed),
        Greater => return Err(SwScoringError::QueryEnded),
        Equal => {}
    }

    match r.cmp(&ref_in_alignment.len()) {
        Less => return Err(SwScoringError::FullReferenceNotUsed),
        Greater => return Err(SwScoringError::ReferenceEnded),
        Equal => {}
    }

    // score is i32, so this cast solely could fail due to it being negative
    match u32::try_from(score) {
        Ok(score) => Ok(score),
        Err(_) => Err(SwScoringError::NegativeScore(score)),
    }
}

/// An enum representing errors that can happen when calculating an alignment
/// score for a particular CIGAR string.
#[derive(PartialEq)]
#[non_exhaustive]
#[deprecated(since = "0.0.37", note = "The new Alignment::rescore_sw returns RescoreSwError instead")]
pub enum SwScoringError {
    /// The CIGAR string produced a negative score
    NegativeScore(i32),
    /// Query ended before the entire CIGAR string was consumed
    QueryEnded,
    /// Reference ended before the entire CIGAR string was consumed
    ReferenceEnded,
    /// Failed to consume the full, provided query, which was expected to
    /// contain no more than what was represented by the CIGAR string
    FullQueryNotUsed,
    /// Failed to consume the full, provided reference, which was expected to
    /// contain only the aligned region of the original reference
    FullReferenceNotUsed,
}

#[allow(deprecated)]
impl Display for SwScoringError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SwScoringError::NegativeScore(score) => write!(f, "The alignment produced a negative score: {score}"),
            SwScoringError::QueryEnded => write!(f, "The query ended before the entire CIGAR string was consumed!"),
            SwScoringError::ReferenceEnded => write!(f, "The reference ended before the entire CIGAR string was consumed!"),
            SwScoringError::FullQueryNotUsed => write!(
                f,
                "Failed to consume the full, provided query, which was expected to contain no more than what was represented by the CIGAR string"
            ),
            SwScoringError::FullReferenceNotUsed => {
                write!(
                    f,
                    "Failed to consume the full, provided reference, which was expected to contain only the aligned region of the original reference"
                )
            }
        }
    }
}

#[allow(deprecated)]
impl Debug for SwScoringError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

#[allow(deprecated)]
impl Error for SwScoringError {}
