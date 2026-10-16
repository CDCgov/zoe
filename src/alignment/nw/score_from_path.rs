use crate::{
    alignment::{AlignOp, ScalarProfile},
    data::cigar::Ciglet,
};
use std::{
    cmp::Ordering::{Equal, Greater, Less},
    error::Error,
    fmt::{Debug, Display},
};

/// Computes the score for a global alignment.
///
/// When scoring an alignment, this function expects the `query` to be passed as
/// a [`ScalarProfile`].
///
/// ## Validity
///
/// - The iterator of ciglets should never contain two ciglets with the same
///   operation adjacent to each other.
///
/// ## Errors
///
/// All of `query`, `ref_in_alignment`, and `cigar` must be fully consumed.
///
/// ## Panics
///
/// - All increments must be nonzero.
#[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
pub fn nw_score_from_path<const S: usize>(
    ciglets: impl IntoIterator<Item = Ciglet>, reference: &[u8], query: &ScalarProfile<S>,
) -> Result<i32, NwScoringError> {
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
                    let Some(query_base) = query.seq.get(q).copied() else {
                        return Err(NwScoringError::QueryEnded);
                    };
                    score += i32::from(query.matrix.get_weight(reference_base, query_base));
                    q += 1;
                    r += 1;
                }
            }
            AlignOp::I => {
                score += query.gap_open + query.gap_extend * (inc - 1) as i32;
                q += inc;
            }
            AlignOp::D => {
                score += query.gap_open + query.gap_extend * (inc - 1) as i32;
                r += inc;
            }
            AlignOp::N => r += inc,
            AlignOp::P => {}
            AlignOp::H | AlignOp::S => return Err(NwScoringError::Clipping),
        }
    }

    match q.cmp(&query.seq.len()) {
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

impl Debug for NwScoringError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for NwScoringError {}
