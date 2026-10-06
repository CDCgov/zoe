use crate::{
    alignment::{ScalarProfile, ScoringError},
    data::cigar::Ciglet,
};
use std::cmp::Ordering::{Equal, Greater, Less};

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
) -> Result<u32, ScoringError> {
    let mut score = 0;
    let mut r = 0;
    let mut q = 0;

    for Ciglet { inc, op } in ciglets {
        match op {
            b'M' | b'=' | b'X' => {
                for _ in 0..inc {
                    let Some(reference_base) = ref_in_alignment.get(r).copied() else {
                        return Err(ScoringError::ReferenceEnded);
                    };
                    let Some(query_base) = query.seq.get(q).copied() else {
                        return Err(ScoringError::QueryEnded);
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
            op => return Err(ScoringError::InvalidCigarOp(op)),
        }
    }

    match q.cmp(&query.seq.len()) {
        Less => return Err(ScoringError::FullQueryNotUsed),
        Greater => return Err(ScoringError::QueryEnded),
        Equal => {}
    }

    match r.cmp(&ref_in_alignment.len()) {
        Less => return Err(ScoringError::FullReferenceNotUsed),
        Greater => return Err(ScoringError::ReferenceEnded),
        Equal => {}
    }

    // score is i32, so this cast solely could fail due to it being negative
    match u32::try_from(score) {
        Ok(score) => Ok(score),
        Err(_) => Err(ScoringError::NegativeScore(score)),
    }
}
