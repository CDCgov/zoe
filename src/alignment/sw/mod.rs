//! ## Smith-Waterman Alignment
//!
//! For generating the optimal, local score use [`sw_simd_score`] or a wrapper
//! thereof. Likewise, for generating an optimal, local alignment use
//! [`sw_simd_alignment`].
//!
//! ### Affine Gap Penalties
//!
//! We use the affine gap formula, $W(k) = u(k-1) + v$, where $k$ is the gap
//! length, $u$ is the gap extend penalty, $v$ is the gap open penalty, and
//! $W(k)$ is the total penalty for the gap. In order to use $W(k) = uk + v$,
//! simply pass gap open as the gap open plus gap extension.
//!
//! ### Usage Note
//!
//! The following steps may be needed:
//!
//! 1. Choose an alphabet. *Zoe* provides an easy interface to work with DNA,
//!    assuming the bases `ACGT` are case-insensitive and `N` is used as a
//!    catch-all. Otherwise, you can define your own alphabet using
//!    [`ByteIndexMap::new`].
//!
//! 2. Specify the [`WeightMatrix`] used for scoring matches and mismatches. For
//!    DNA, [`new_dna_matrix`] is a convenient constructor.
//!
//! 3. Build the query profile with [`StripedProfile::new`]. This step combines
//!    the query, matrix of weights, and gap open and gap extend penalties. This
//!    step also performs some basic checks, and if any of the inputs are
//!    invalid, a [`QueryProfileError`] is returned.
//!
//! 4. Use the query profile to align against any number of different
//!    references. For the score only, use
//!    [`StripedProfile::smith_waterman_score`]. For the full alignment, use
//!    [`StripedProfile::smith_waterman_alignment`].
//!
//! Note: Scalar versions ([`sw_scalar_score`], [`sw_scalar_alignment`], and
//! [`ScalarProfile`]) are also available for fuzzing and testing purposes.
//!
//! Below is an example using DNA. We use a match score of 4 and a mismatch
//! score of -2, as defined in `WEIGHTS`. We also choose to use a
//! `StripedProfile` (so that the SIMD algorithm is used) with integer type `i8`
//! and 32 lanes.
//!
//! ```
//! # use zoe::{
//! #     alignment::{Alignment, AlignmentStates, StripedProfile},
//! #     data::matrices::WeightMatrix
//! # };
//! let reference: &[u8] = b"GGCCACAGGATTGAG";
//! let query: &[u8] = b"CTCAGATTG";
//!
//! const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
//! const GAP_OPEN: i8 = -3;
//! const GAP_EXTEND: i8 = -1;
//!
//! let profile = StripedProfile::<i8, 32, 5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let alignment = profile.smith_waterman_alignment(reference).unwrap();
//!
//! let Alignment {
//!     score,
//!     ref_range,
//!     query_range,
//!     states,
//!     ..
//! } = alignment;
//! assert_eq!(score, 27);
//! assert_eq!(states, AlignmentStates::try_from(b"5M1D4M").unwrap());
//! ```
//!
//! Below is an example using a different alphabet. Matches are given a score of
//! 1 and mismatches are given a score of -1.
//! ```
//! # use zoe::{
//! #     alignment::{Alignment, AlignmentStates, StripedProfile},
//! #     data::{matrices::WeightMatrix, ByteIndexMap},
//! # };
//! let reference: &[u8] = b"BDAACAABDDDB";
//! let query: &[u8] = b"AABDDAB";
//!
//! const MAPPING: ByteIndexMap<4> = ByteIndexMap::new(*b"ABCD", b'A');
//! const WEIGHTS: WeightMatrix<i8, 4> = WeightMatrix::new(&MAPPING, 1, -1, None);
//! const GAP_OPEN: i8 = -4;
//! const GAP_EXTEND: i8 = -2;
//!
//! let profile = StripedProfile::<i8, 32, 4>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let alignment = profile.smith_waterman_alignment(reference).unwrap();
//!
//! let Alignment {
//!     score,
//!     ref_range,
//!     query_range,
//!     states,
//!     ..
//! } = alignment;
//! assert_eq!(score, 5);
//! assert_eq!(states, AlignmentStates::try_from(b"5M2S").unwrap())
//! ```
//!
//! When using the SIMD algorithm, you must specify the number of lanes `N` and
//! integer type `T` (typically i8, i16, or i32). If the integer type has too
//! small of a range, it is possible the alignment score will overflow (in which
//! case `None` is returned). A higher-level abstraction to avoid this uses
//! profile sets: [`LocalProfiles`] or [`SharedProfiles`].
//!
//! The former is designed for use within a single thread, while the latter
//! allows multiple threads to access it. Both store a set of lazily-evaluated
//! query profiles for `i8`, `i16`, and `i32`. To create one of these, you can
//! call [`LocalProfiles::new_with_w256`], [`SharedProfiles::new_with_w256`], or
//! one of the other constructors. Then, call
//! [`LocalProfiles::smith_waterman_alignment_from_i8`],
//! [`SharedProfiles::smith_waterman_alignment_from_i8`], or one of the other
//! methods.
//!
//! When using DNA, you can also create a profile by using
//! [`Nucleotides::into_local_profile`] or [`Nucleotides::into_shared_profile`].
//!
//! Below is the previous example using DNA, but with this higher-level API:
//! ```
//! # use zoe::{
//! #     alignment::{Alignment, AlignmentStates, ProfileSets},
//! #     data::matrices::WeightMatrix,
//! #     prelude::Nucleotides
//! # };
//! let reference: &[u8] = b"GGCCACAGGATTGAG";
//! let query: Nucleotides = b"CTCAGATTG".into();
//!
//! const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
//! const GAP_OPEN: i8 = -3;
//! const GAP_EXTEND: i8 = -1;
//!
//! let profile = query.into_local_profile(&WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let alignment = profile.smith_waterman_alignment_from_i8(reference).unwrap();
//!
//! let Alignment {
//!     score,
//!     ref_range,
//!     query_range,
//!     states,
//!     ..
//! } = alignment;
//! assert_eq!(score, 27);
//! assert_eq!(states, AlignmentStates::try_from(b"5M1D4M").unwrap());
//! ```
//!
//! And similarly with the different alphabet:
//! ```
//! # use zoe::{
//! #     alignment::{Alignment, AlignmentStates, LocalProfiles, ProfileSets},
//! #     data::{matrices::WeightMatrix, ByteIndexMap},
//! # };
//! let reference: &[u8] = b"BDAACAABDDDB";
//! let query: &[u8] = b"AABDDAB";
//!
//! const MAPPING: ByteIndexMap<4> = ByteIndexMap::new(*b"ABCD", b'A');
//! const WEIGHTS: WeightMatrix<i8, 4> = WeightMatrix::new(&MAPPING, 1, -1, None);
//! const GAP_OPEN: i8 = -4;
//! const GAP_EXTEND: i8 = -2;
//!
//! let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let alignment = profile.smith_waterman_alignment_from_i8(reference).unwrap();
//!
//! let Alignment {
//!     score,
//!     ref_range,
//!     query_range,
//!     states,
//!     ..
//! } = alignment;
//! assert_eq!(score, 5);
//! assert_eq!(states, AlignmentStates::try_from(b"5M2S").unwrap());
//! ```
//!
//! ## Module Citations
//!
//! 1. Smith, Temple F. & Waterman, Michael S. (1981). "Identification of Common
//!    Molecular Subsequences" (PDF). Journal of Molecular Biology. 147 (1):
//!    195–197.
//!
//! 2. Osamu Gotoh (1982). "An improved algorithm for matching biological
//!    sequences". Journal of Molecular Biology. 162 (3): 705–708.
//!
//! 3. Stephen F. Altschul & Bruce W. Erickson (1986). "Optimal sequence
//!    alignment using affine gap costs". Bulletin of Mathematical Biology. 48
//!    (5–6): 603–616.
//!
//! 4. Tomáš Flouri, Kassian Kobert, Torbjørn Rognes, Alexandros
//!    Stamatakis(2015). "Are all global alignment algorithms and
//!    implementations correct?" bioRxiv 031500. doi:
//!    <https://doi.org/10.1101/031500>
//!
//! 5. Farrar, Michael (2006). "Striped Smith-Waterman speeds database searches
//!    six times over other SIMD implementations". Bioinformatics, 23(2),
//!    156-161. doi: <https://doi.org/10.1093/bioinformatics/btl582>
//!
//! 6. Szalkowski, Adam, Ledergerber, Christian, Krähenbühl, Philipp & Dessimoz,
//!    Christophe (2008). "SWPS3 - fast multi-threaded vectorized Smith-Waterman
//!    for IBM Cell/B.E. and x86/SSE2". BMC Research Notes. 1:
//!    107. doi: <https://doi.org/10.1186/1756-0500-1-107>
//!
//! 7. Zhao, Mengyao, Lee, Wan-Ping, Garrison, Erik P. & Marth, Gabor T. (2013).
//!    "SSW library: an SIMD Smith-Waterman C/C++ library for use in genomic
//!    applications". `PLoS One`. 8(12): e82138. doi:
//!    <https://doi.org/10.1371/journal.pone.0082138>
//!
//! 8. Daily, Jeff (2016). "Parasail: SIMD C library for global, semi-global,
//!    and local pairwise sequence alignments". BMC Bioinformatics. 17: 81. doi:
//!    <https://doi.org/10.1186/s12859-016-0930-z>
//!
//!
//! [`sw_scalar_score`]: sw::sw_scalar_score
//! [`sw_simd_score`]: sw::sw_simd_score
//! [`sw_scalar_alignment`]: sw::sw_scalar_alignment
//! [`sw_simd_alignment`]: sw::sw_simd_alignment
//! [`ByteIndexMap::new`]: crate::data::ByteIndexMap::new
//! [`WeightMatrix`]: crate::data::matrices::WeightMatrix
//! [`WeightMatrix::new`]: crate::data::matrices::WeightMatrix::new
//! [`new_biased_dna_matrix`]: crate::data::matrices::WeightMatrix::new_biased_dna_matrix
//! [`new_dna_matrix`]: crate::data::matrices::WeightMatrix::new_dna_matrix
//! [`QueryProfileError`]: crate::data::err::QueryProfileError
//! [`Nucleotides::into_local_profile`]:
//!     crate::data::types::nucleotides::Nucleotides::into_local_profile
//! [`Nucleotides::into_shared_profile`]:
//!     crate::data::types::nucleotides::Nucleotides::into_shared_profile
//! [`ScalarProfile::new`]: ScalarProfile::new
//! [`StripedProfile::new`]: StripedProfile::new
//! [`ScalarProfile::smith_waterman_score`]: ScalarProfile::smith_waterman_score
//! [`StripedProfile::smith_waterman_score`]:
//!     StripedProfile::smith_waterman_score
//! [`LocalProfiles`]: LocalProfiles
//! [`SharedProfiles`]: SharedProfiles
//! [`LocalProfiles::smith_waterman_score_from_i8`]:
//!     LocalProfiles::smith_waterman_score_from_i8
//! [`SharedProfiles::smith_waterman_score_from_i8`]:
//!     SharedProfiles::smith_waterman_score_from_i8

#![allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
use super::*;
use crate::data::cigar::Ciglet;
use std::{
    cmp::Ordering::{Equal, Greater, Less},
    simd::prelude::*,
};

/// Compute the Smith-Waterman score for a given alignment.
///
/// When scoring an alignment, this function expects the full `query` to be
/// passed as a [`ScalarProfile`], as well as the slice of the reference which
/// was aligned to (e.g., the indices in [`Alignment::ref_range`]).
///
/// ## Errors
///
/// - The `ciglets` must contain valid operations
/// - All of `query`, `ref_in_alignment`, and `cigar` must be fully consumed
/// - The final score should be nonnegative
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

#[cfg(test)]
pub(crate) mod test_data {
    use super::ScalarProfile;
    use crate::data::matrices::WeightMatrix;
    use std::sync::LazyLock;

    pub(crate) static REFERENCE: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/NC_007362.1.txt")); // H5 HA
    pub(crate) static QUERY: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/NC_026433.1.txt")); // H1 HA1

    pub(crate) static WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(2, -5, Some(b'N'));
    pub(crate) static BIASED_WEIGHTS: WeightMatrix<u8, 5> = WEIGHTS.to_biased_matrix();
    pub(crate) static SCALAR_PROFILE: LazyLock<ScalarProfile<5>> =
        LazyLock::new(|| ScalarProfile::new(QUERY, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap());

    pub(crate) const GAP_OPEN: i8 = -10;
    pub(crate) const GAP_EXTEND: i8 = -1;
}

#[cfg(test)]
mod test;

#[cfg(test)]
mod bench;

mod scalar;
mod striped;

pub use scalar::*;
pub use striped::*;
