//! ## Needleman–Wunsch Alignment
//!
//! For generating the optimal, global score use [`nw_scalar_score`]. For
//! computing the global alignment between, use [`nw_scalar_align`]. To then
//! obtain the aligned sequences (with gap characters inserted), use
//! [`Alignment::get_aligned_seqs`].
//!
//! ### Affine Gap Penalties
//!
//! We use the affine gap formula, $W(k) = u(k-1) + v$, where $k$ is the gap
//! length, $u$ is the gap extend penalty, $v$ is the gap open penalty, and
//! $W(k)$ is the total penalty for the gap. In order to use $W(k) = uk + v$,
//! set the gap open score as your desired gap open score plus the gap extension
//! score.
//!
//! ### Usage Note
//!
//! The following steps may be needed to perform an alignment:
//!
//! 1. Choose an alphabet. *Zoe* provides an easy interface to work with DNA,
//!    assuming the bases `ACGT` are case-insensitive and `N` is used as a
//!    catch-all. Otherwise, you can define your own alphabet using
//!    [`ByteIndexMap::new`].
//!
//! 2. Specify the [`WeightMatrix`] used for scoring matches and mismatches. For
//!    DNA, [`new_dna_matrix`] is a convenient constructor.
//!
//! 3. Build the query profile with [`ScalarProfile::new`]. This step combines
//!    the query, matrix of weights, and gap open and gap extend penalties. This
//!    step also performs some basic checks, and if any of the inputs are
//!    invalid, a [`ProfileError`] is returned.
//!
//! 4. Use the query profile to align against any number of different references
//!    using [`nw_scalar_score`] or [`nw_scalar_align`].
//!
//! Below is an example using DNA. We use a match score of 4 and a mismatch
//! score of -2, as defined in `WEIGHTS`.
//!
//! ```
//! # use zoe::{
//! #     alignment::{Alignment, AlignmentStates, ScalarProfile, nw::nw_scalar_align},
//! #     data::matrices::WeightMatrix
//! # };
//! let reference: &[u8] = b"GGCCACAGGATTGAG";
//! let query: &[u8] = b"TACCACAGTATTAG";
//!
//! const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
//! const GAP_OPEN: i8 = -3;
//! const GAP_EXTEND: i8 = -1;
//!
//! let profile = ScalarProfile::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let alignment = nw_scalar_align(&reference, &profile);
//!
//! let Alignment {
//!     score,
//!     ref_range,
//!     query_range,
//!     states,
//!     ..
//! } = alignment;
//! assert_eq!(score, 35);
//! assert_eq!(states, AlignmentStates::try_from(b"12M1D2M").unwrap());
//! ```
//!
//! [`ByteIndexMap::new`]: crate::data::ByteIndexMap::new
//! [`WeightMatrix`]: crate::data::matrices::WeightMatrix
//! [`WeightMatrix::new`]: crate::data::matrices::WeightMatrix::new
//! [`new_dna_matrix`]: crate::data::matrices::WeightMatrix::new_dna_matrix
//! [`ProfileError`]: crate::alignment::ProfileError
//! [`Alignment`]: super::Alignment
//! [`Alignment::get_aligned_seqs`]: super::Alignment::get_aligned_seqs

use crate::alignment::{Alignment, BackTrackable, BacktrackMatrix, ScalarProfile};
use std::ops::Add;

/// Needleman–Wunsch algorithm (non-vectorized), yielding the optimal score.
///
/// Provides the globally optimal sequence alignment score using affine gap
/// penalties. Our implementation is derived as an adaptation of
/// [`sw_scalar_score`].
///
/// ## Complexity
///
/// For query length $m$ and reference length $n$:
///
/// - Time: $O(mn)$
/// - Space: $O(m)$
///
/// ## Limitations
///
/// - Not suitable for production due to high runtimes.
/// - Useful for testing correctness.
///
/// ## Example
///
/// ```
/// # use zoe::{
/// #     alignment::{Alignment, AlignmentStates, ScalarProfile, nw::nw_scalar_score},
/// #     data::matrices::WeightMatrix
/// # };
/// let reference: &[u8] = b"GGCCACAGGATTGAG";
/// let query: &[u8] = b"TACCACAGTATTAG";
///
/// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
/// const GAP_OPEN: i8 = -3;
/// const GAP_EXTEND: i8 = -1;
///
/// let profile = ScalarProfile::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
/// let score = nw_scalar_score(&reference, &profile);
///
/// assert_eq!(score, 35);
/// ```
///
/// [`sw_scalar_score`]: super::sw::sw_scalar_score
#[must_use]
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
pub fn nw_scalar_score<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> i32 {
    // Definitions:
    // * Q[c] = c-th base in the query, zero-indexed
    // * R[r] = r-th base in the reference, zero-indexed
    // * W[r,c] = contribution to score of aligning Q[c] and R[r]
    // * H[r,c] = maximum score for aligning Q[0..c] to R[0..r]
    // * E[r,c] = maximum score for aligning Q[0..c] to R[0..r] such that the
    //   alignment ends with a gap consuming a base in the reference
    // * F[r,c] = maximum score for aligning Q[0..c] to R[0..r] such that the
    //   alignment ends with a gap consuming a base in the query

    // A row H[r, ..] of H. Initialized for row -1 (no bases in the reference
    // consumed). H[-1, 0] is gap_open (consuming one base of query). Each
    // subsequent value adds gap_extend to the previous.
    let mut h_row = vec![0; query.seq.len()];
    let mut h = query.gap_open;
    for h_val in &mut h_row {
        *h_val = h;
        h += query.gap_extend;
    }

    // A row E[r, ..] of E. Initialized for row 0 (one base in the reference
    // consumed, via a gap). E[0,0] is only reached by consuming a base in the
    // query via a gap_open, then conuming a base in the reference via a
    // gap_open. Each subsequent value adds gap_extend to the previous.
    let mut e_row = vec![0; query.seq.len()];
    let mut e = 2 * query.gap_open;
    for e_val in &mut e_row {
        *e_val = e;
        e += query.gap_extend;
    }

    // Iterate over R (index=r)
    for (r, reference_base) in reference.iter().copied().enumerate() {
        // A value in F, initialized to F[r,0]. F[r,0] is reached by consuming
        // r+1 bases in the reference via a gap, then consuming a base in the
        // query via a gap.
        let mut f = 2 * query.gap_open + r as i32 * query.gap_extend;
        // A value of H. Initialized to H[r-1,-1] (r bases consumed in R and no
        // bases consumed in Q)
        let mut h = if r == 0 {
            // H[-1,-1]: The alignment has just begun (no bases consumed in R or
            // Q)
            0
        } else {
            // H[r-1,-1]: Reached by a gap in the reference consuming r bases
            query.gap_open + query.gap_extend * (r - 1) as i32
        };

        // Iterate over Q (index=c). Assumes that h contains H[r-1,c-1].
        for c in 0..query.seq.len() {
            // W[r,c]
            let match_score = i32::from(query.matrix.get_weight(reference_base, query.seq[c]));
            // Temporarily update to H[r-1,c-1] + W[r,c]
            h += match_score;

            // Extract E[r,c]
            let mut e = e_row[c];
            // Update H[r-1,c-1] to H[r,c] by computing the maximum of
            // H[r-1,c-1] + W[r,c], E[r,c], and F[r,c]
            h = h.max(e).max(f);

            // Update E[r,c] to E[r,c+1]
            e = e.add(query.gap_extend).max(h + query.gap_open);
            // Update F[r,c] to F[r,c+1]
            f = f.add(query.gap_extend).max(h + query.gap_open);

            // Update H[r,c] in h_row, and set h to H[r-1,c] (becomes H[r-1,c-1]
            // next iteration)
            (h, h_row[c]) = (h_row[c], h);
            // Update E[r,c] in e_row
            e_row[c] = e;
        }
    }

    h_row[h_row.len() - 1]
}

/// Needleman–Wunsch algorithm (non-vectorized), yielding the optimal alignment.
///
/// Provides the globally optimal sequence alignment using affine gap penalties.
/// Our implementation is derived as an adaptation of [`sw_scalar_align`].
///
/// ## Complexity
///
/// For query length $m$ and reference length $n$:
///
/// - Time: $O(mn)$
/// - Space: $O(mn)$
///
/// ## Limitations
///
/// - Not suitable for production due to high runtimes.
/// - Useful for testing correctness.
///
/// ## Example
///
/// ```
/// # use zoe::{
/// #     alignment::{Alignment, AlignmentStates, ScalarProfile, nw::nw_scalar_align},
/// #     data::matrices::WeightMatrix
/// # };
/// let reference: &[u8] = b"GGCCACAGGATTGAG";
/// let query: &[u8] = b"TACCACAGTATTAG";
///
/// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
/// const GAP_OPEN: i8 = -3;
/// const GAP_EXTEND: i8 = -1;
///
/// let profile = ScalarProfile::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
/// let alignment = nw_scalar_align(&reference, &profile);
///
/// let Alignment {
///     score,
///     ref_range,
///     query_range,
///     states,
///     ..
/// } = alignment;
/// assert_eq!(score, 35);
/// assert_eq!(states, AlignmentStates::try_from(b"12M1D2M").unwrap());
/// ```
///
/// [`sw_scalar_align`]: super::sw::sw_scalar_align
#[must_use]
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
pub fn nw_scalar_align<const S: usize>(reference: &[u8], query: &ScalarProfile<S>) -> Alignment<i32> {
    // See dev comments in nw_scalar_score for more details

    let mut h_row = vec![0; query.seq.len()];
    let mut h = query.gap_open;
    for h_val in &mut h_row {
        *h_val = h;
        h += query.gap_extend;
    }

    let mut e_row = vec![0; query.seq.len()];
    let mut e = 2 * query.gap_open;
    for e_val in &mut e_row {
        *e_val = e;
        e += query.gap_extend;
    }

    let mut backtrack = BacktrackMatrix::new(reference.len(), query.seq.len());

    for (r, reference_base) in reference.iter().copied().enumerate() {
        let mut f = 2 * query.gap_open + r as i32 * query.gap_extend;
        let mut h = if r == 0 {
            0
        } else {
            query.gap_open + query.gap_extend * (r - 1) as i32
        };

        for c in 0..query.seq.len() {
            backtrack.move_to(r, c);

            let match_score = i32::from(query.matrix.get_weight(reference_base, query.seq[c]));
            h += match_score;

            let mut e = e_row[c];
            h = h.max(e).max(f);

            if e == h {
                backtrack.up();
            }

            if f == h {
                backtrack.left();
            }

            // H[r-1,c]
            let next_diag = h_row[c];
            // Update H[r,c] in h_row
            h_row[c] = h;

            // Temporarily store H[r,c]+go in h
            h += query.gap_open;
            // Update E[r,c] to E[r,c+1] by computing the maximum of E[r,c]+ge
            // and H[r,c]+go
            e = e.add(query.gap_extend).max(h);
            // Update F[r,c] to F[r,c+1] by computing the maximum of F[r,c]+ge
            // and F[r,c]+go
            f = f.add(query.gap_extend).max(h);

            // Mark the backtrack as up_extending only if it cannot be reached
            // via opening a gap. This check ensures E[r,c+1] > H[r,c]+go, so
            // that E[r,c+1] = E[r,c]+ge.
            if e > h {
                backtrack.up_extending();
            }

            // Mark the backtrack as left_extending only if it cannot be reached
            // via opening a gap. This check ensures F[r,c+1] > H[r,c]+go, so
            // that F[r,c+1] = F[r,c]+ge.
            if f > h {
                backtrack.left_extending();
            }

            // Store H[r-1,c] in h (becomes H[r-1,c-1] next iteration)
            h = next_diag;
            e_row[c] = e;
        }
    }

    backtrack.to_alignment_global(h_row[h_row.len() - 1], reference.len(), query.seq.len())
}

/// Computes the score for a global alignment.
///
/// When scoring an alignment, this function expects the `query` to be
/// passed as a [`ScalarProfile`].
///
/// ## Validity
///
/// - The iterator of ciglets should never contain two ciglets with the same
///   operation adjacent to each other.
///
/// ## Errors
///
/// - The `ciglets` must contain valid operations in `MIDNP=X`,
/// - All of `query`, `ref_in_alignment`, and `cigar` must be fully consumed.
///
/// ## Panics
///
/// - All increments must be nonzero.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// access this function.
///
/// </div>
#[cfg(feature = "alignment-diagnostics")]
#[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
pub fn nw_score_from_path<const S: usize>(
    ciglets: impl IntoIterator<Item = crate::data::cigar::Ciglet>, reference: &[u8], query: &ScalarProfile<S>,
) -> Result<i32, super::ScoringError> {
    use std::cmp::Ordering::{Equal, Greater, Less};

    let mut score = 0;
    let mut r = 0;
    let mut q = 0;

    for crate::data::cigar::Ciglet { inc, op } in ciglets {
        match op {
            b'M' | b'=' | b'X' => {
                for _ in 0..inc {
                    let Some(reference_base) = reference.get(r).copied() else {
                        return Err(super::ScoringError::ReferenceEnded);
                    };
                    let Some(query_base) = query.seq.get(q).copied() else {
                        return Err(super::ScoringError::QueryEnded);
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
            b'N' => r += inc,
            b'P' => {}
            op => return Err(super::ScoringError::InvalidCigarOp(op)),
        }
    }
    match q.cmp(&query.seq.len()) {
        Less => return Err(super::ScoringError::FullQueryNotUsed),
        Greater => return Err(super::ScoringError::QueryEnded),
        Equal => {}
    }

    match r.cmp(&reference.len()) {
        Less => return Err(super::ScoringError::FullReferenceNotUsed),
        Greater => return Err(super::ScoringError::ReferenceEnded),
        Equal => {}
    }

    Ok(score)
}
