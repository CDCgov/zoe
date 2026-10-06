use crate::{alignment::ScalarProfile, data::cigar::Ciglet};
use std::{
    cmp::Ordering::{Equal, Greater, Less},
    error::Error,
    fmt::{Debug, Display},
};

/// Computes the score for a local alignment.
///
/// When scoring an alignment, this function expects the full `query` to be
/// passed as a [`ScalarProfile`], as well as the slice of the reference which
/// was aligned to (e.g., the indices in [`Alignment::ref_range`]).
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
pub fn sw_score_from_path<const S: usize>(
    ciglets: impl IntoIterator<Item = Ciglet>, ref_in_alignment: &[u8], query: &ScalarProfile<S>,
) -> Result<u32, SwScoringError> {
    let mut score = 0;
    let mut r = 0;
    let mut q = 0;

    for Ciglet { inc, op } in ciglets {
        match op {
            b'M' | b'=' | b'X' => {
                for _ in 0..inc {
                    let Some(reference_base) = ref_in_alignment.get(r).copied() else {
                        return Err(SwScoringError::ReferenceEnded);
                    };
                    let Some(query_base) = query.seq.get(q).copied() else {
                        return Err(SwScoringError::QueryEnded);
                    };
                    score += i32::from(query.matrix.get_weight(reference_base, query_base));
                    q += 1;
                    r += 1;
                }
            }
            b'I' => {
                score += query.gap_open + query.gap_extend * (inc - 1) as i32;
                q += inc;
            }
            b'D' => {
                score += query.gap_open + query.gap_extend * (inc - 1) as i32;
                r += inc;
            }
            b'S' => q += inc,
            b'N' => r += inc,
            b'H' | b'P' => {}
            op => return Err(SwScoringError::InvalidCigarOp(op)),
        }
    }

    match q.cmp(&query.seq.len()) {
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
    /// Unsupported CIGAR opcode used in argument
    InvalidCigarOp(u8),
}

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
            SwScoringError::InvalidCigarOp(op) => write!(f, "An unsupported CIGAR opcode was encountered: {op}"),
        }
    }
}

impl Debug for SwScoringError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl Error for SwScoringError {}
