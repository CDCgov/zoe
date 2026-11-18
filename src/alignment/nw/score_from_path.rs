use crate::{
    alignment::{AlignOp, Alignment, ScoreParams},
    data::{cigar::Ciglet, err::impl_err_conv},
};
use std::{
    cmp::Ordering::{Equal, Greater, Less},
    error::Error,
    fmt::{Debug, Display},
};

impl<T> Alignment<T> {
    /// Re-computes the score for a global alignment. The score field of the
    /// alignment is not utilized (and could be any type, such as `()`).
    ///
    /// The full reference and query corresponding to the alignment must be
    /// passed.
    ///
    /// ## Errors
    ///
    /// See [`RescoreNwError`].
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
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    pub fn rescore_nw<const S: usize>(
        &self, reference: &[u8], query: &[u8], params: &ScoreParams<i8, S>,
    ) -> Result<i32, RescoreNwError> {
        let ScoreParams {
            matrix,
            gap_open,
            gap_extend,
        } = params;
        let [gap_open, gap_extend] = [*gap_open, *gap_extend].map(i32::from);
        let matrix = matrix.as_ref();

        if query.len() != self.query_len {
            return Err(RescoreNwError::QueryLen {
                expected: self.query_len,
                actual:   query.len(),
            });
        }
        if self.query_range.start != 0 {
            return Err(RescoreNwError::QueryRangeStart {
                start: self.query_range.start,
            });
        }
        if self.query_range.end != self.query_len {
            return Err(RescoreNwError::QueryRangeEnd {
                end: self.query_range.end,
                len: self.query_len,
            });
        }
        if reference.len() != self.ref_len {
            return Err(RescoreNwError::RefLen {
                expected: self.ref_len,
                actual:   reference.len(),
            });
        }
        if self.ref_range.start != 0 {
            return Err(RescoreNwError::RefRangeStart {
                start: self.ref_range.start,
            });
        }
        if self.ref_range.end != self.ref_len {
            return Err(RescoreNwError::RefRangeEnd {
                end: self.ref_range.end,
                len: self.ref_len,
            });
        }

        let mut score = 0;
        let mut r = 0;
        let mut q = 0;

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
                AlignOp::P => {}
                AlignOp::H | AlignOp::S => return Err(RescoreNwError::Clipping),
            }
        }

        Ok(score)
    }
}

/// An error type for [`Alignment::rescore_nw`].
#[derive(Eq, PartialEq)]
pub enum RescoreNwError {
    /// The length of the passed query was incorrect.
    QueryLen { expected: usize, actual: usize },
    /// The query range did not start at 0.
    QueryRangeStart { start: usize },
    /// The query range did not reach the end of the query.
    QueryRangeEnd { end: usize, len: usize },
    /// The length of the passed reference was incorrect.
    RefLen { expected: usize, actual: usize },
    /// The reference range did not start at 0.
    RefRangeStart { start: usize },
    /// The reference range did not reach the end of the reference.
    RefRangeEnd { end: usize, len: usize },
    /// Clipping was found inside the alignment states.
    Clipping,
}

impl Display for RescoreNwError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match *self {
            RescoreNwError::QueryLen { expected, actual } => write!(
                f,
                "The alignment was for a query of length {expected}, but a query of length {actual} was passed"
            ),
            RescoreNwError::QueryRangeStart { start } => write!(
                f,
                "Expected the query range to start at 0 for a global alignment, but found {start}"
            ),
            RescoreNwError::QueryRangeEnd { end, len } => write!(
                f,
                "Expected the query range to end at the query length ({len}) for a global alignment, but found {end}"
            ),
            RescoreNwError::RefLen { expected, actual } => write!(
                f,
                "The alignment was for a reference of length {expected}, but a reference of length {actual} was passed"
            ),
            RescoreNwError::RefRangeStart { start } => write!(
                f,
                "Expected the reference range to start at 0 for a global alignment, but found {start}"
            ),
            RescoreNwError::RefRangeEnd { end, len } => write!(
                f,
                "Expected the reference range to end at the reference length ({len}) for a global alignment, but found {end}"
            ),
            RescoreNwError::Clipping => write!(f, "Clipping was found inside the alignment states"),
        }
    }
}

impl Debug for RescoreNwError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for RescoreNwError {}

impl_err_conv!(RescoreNwError);

/// Computes the score for a global alignment.
///
/// ## Validity
///
/// The iterator of ciglets should never contain two ciglets with the same
/// operation adjacent to each other.
///
/// ## Errors
///
/// All of `query`, `ref_in_alignment`, and `cigar` must be fully consumed.
///
/// ## Panics
///
/// All increments must be nonzero.
#[deprecated(since = "0.0.37", note = "Use Alignment::rescore_nw instead")]
#[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation, deprecated)]
pub fn nw_score_from_path<const S: usize>(
    ciglets: impl IntoIterator<Item = Ciglet>, reference: &[u8], query: &[u8], params: &ScoreParams<i8, S>,
) -> Result<i32, NwScoringError> {
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
                    let Some(reference_base) = reference.get(r).copied() else {
                        return Err(NwScoringError::ReferenceEnded);
                    };
                    let Some(query_base) = query.get(q).copied() else {
                        return Err(NwScoringError::QueryEnded);
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
            AlignOp::N => r += inc,
            AlignOp::P => {}
            AlignOp::H | AlignOp::S => return Err(NwScoringError::Clipping),
        }
    }

    match q.cmp(&query.len()) {
        Less => return Err(NwScoringError::FullQueryNotUsed),
        Greater => return Err(NwScoringError::QueryEnded),
        Equal => {}
    }

    match r.cmp(&reference.len()) {
        Less => return Err(NwScoringError::FullReferenceNotUsed),
        Greater => return Err(NwScoringError::ReferenceEnded),
        Equal => {}
    }

    Ok(score)
}

/// An enum representing errors that can happen when calculating an alignment
/// score for a particular CIGAR string.
#[derive(PartialEq)]
#[non_exhaustive]
#[deprecated(since = "0.0.37", note = "The new Alignment::rescore_nw returns RescoreNwError instead")]
pub enum NwScoringError {
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
    /// Clipping was found inside the global alignment
    Clipping,
}

#[allow(deprecated)]
impl Display for NwScoringError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            NwScoringError::QueryEnded => write!(f, "The query ended before the entire CIGAR string was consumed!"),
            NwScoringError::ReferenceEnded => write!(f, "The reference ended before the entire CIGAR string was consumed!"),
            NwScoringError::FullQueryNotUsed => write!(
                f,
                "Failed to consume the full, provided query, which was expected to contain no more than what was represented by the CIGAR string"
            ),
            NwScoringError::FullReferenceNotUsed => {
                write!(
                    f,
                    "Failed to consume the full, provided reference, which was expected to contain only the aligned region of the original reference"
                )
            }
            NwScoringError::Clipping => write!(f, "Clipping was found inside the CIGAR string"),
        }
    }
}

#[allow(deprecated)]
impl Debug for NwScoringError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

#[allow(deprecated)]
impl Error for NwScoringError {}
