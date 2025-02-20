//! ## Functions for aligning sequence data.
//!
//! *Zoe* supports efficient local alignment for DNA, protein, or any other
//! sequence data via the Smith-Waterman algorithm. For generating the score
//! (useful for database searches and determining whether sequences are
//! related), this includes both a traditional implementation (see
//! [`sw_scalar_score`]) and a striped SIMD implementation (see
//! [`sw_simd_score`]). The function [`sw_scalar_alignment`] can be used to
//! generate the alignment itself.
//!
//! To calculate an alignment score or generate an alignment, follow these
//! steps. The first two are designed to be done at compile time, and the
//! subsequent two to be done at runtime.
//!
//! 1. Choose an alphabet. *Zoe* provides an easy interface to work with DNA,
//!    assuming the bases `ACGT` are case-insensitive and `N` is used as a
//!    catch-all. Otherwise, you can define your own alphabet using
//!    [`ByteIndexMap::new`].
//!
//! 2. Specify the [`WeightMatrix`] used for scoring matches and mismatches. For
//!    DNA, [`new_dna_matrix`] is a convenient constructor.
//!
//! 3. Build the query profile, with either [`ScalarProfile::new`] or
//!    [`StripedProfile::new`] (for SIMD). This step combines the query, matrix
//!    of weights, and gap open and gap extend penalties. This step also
//!    performs some basic checks, and if any of the inputs are invalid, a
//!    [`QueryProfileError`] is returned.
//!
//! 4. Use the query profile to align against any number of different
//!    references. Simply call the [`ScalarProfile::smith_waterman_score`] or
//!    [`StripedProfile::smith_waterman_score`] methods. The profile can be
//!    reused in multiple different calls.
//!
//! Below is an example using DNA. We use a match score of 4 and a mismatch
//! score of -2, as defined in `WEIGHTS`. We also choose to use a
//! `StripedProfile` (so that the SIMD algorithm is used) with integer type `i8`
//! and 32 lanes.
//! ```
//! # use zoe::{alignment::StripedProfile, data::WeightMatrix};
//! let reference: &[u8] = b"GGCCACAGGATTGAG";
//! let query: &[u8] = b"CTCAGATTG";
//!
//! const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
//! const GAP_OPEN: i8 = -3;
//! const GAP_EXTEND: i8 = -1;
//!
//! let profile = StripedProfile::<i8, 32, 5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let score = profile.smith_waterman_score(reference).unwrap();
//! assert_eq!(score, 27);
//! ```
//!
//! Below is an example using a different alphabet. Matches are given a score of
//! 1 and mismatches are given a score of -1.
//! ```
//! # use zoe::{
//! #     alignment::StripedProfile,
//! #     data::{WeightMatrix, ByteIndexMap},
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
//! let score = profile.smith_waterman_score(reference).unwrap();
//! assert_eq!(score, 5);
//! ```
//!
//! When using the SIMD algorithm, you must specify the number of lanes `N` and
//! integer type `T` (typically i8, i16, i32, or i64). If the integer type has
//! too small of a range, it is possible the alignment score will overflow (in
//! which case `None` is returned). A higher-level abstraction to avoid this is
//! [`LocalProfiles`] and [`SharedProfiles`].
//!
//! The former is designed for use within a single thread, while the latter
//! allows multiple threads to access it. Both store a set of lazily-evaluated
//! query profiles for `i8`, `i16`, `i32`, and `i64`. To create one of these,
//! you can call [`LocalProfiles::new_with_i8`],
//! [`SharedProfiles::new_with_i8`], or one of the other constructors. Then,
//! call [`LocalProfiles::smith_waterman_score_from_i8`],
//! [`SharedProfiles::smith_waterman_score_from_i8`], or one of the other
//! methods.
//!
//! When using DNA, you can also create a profile by using
//! [`Nucleotides::into_local_profile`] or [`Nucleotides::into_shared_profile`].
//!
//! Below is the previous example using DNA, but with this higher-level API:
//! ```
//! # use zoe::{data::WeightMatrix, prelude::Nucleotides};
//! let reference: &[u8] = b"GGCCACAGGATTGAG";
//! let query: Nucleotides = b"CTCAGATTG".into();
//!
//! const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
//! const GAP_OPEN: i8 = -3;
//! const GAP_EXTEND: i8 = -1;
//!
//! let profile = query.into_local_profile::<32, 5>(&WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let score = profile.smith_waterman_score_from_i8(reference).unwrap();
//! assert_eq!(score, 27);
//! ```
//!
//! And similarly with the different alphabet:
//! ```
//! # use zoe::{
//! #     alignment::LocalProfiles,
//! #     data::{WeightMatrix, ByteIndexMap},
//! # };
//! let reference: &[u8] = b"BDAACAABDDDB";
//! let query: &[u8] = b"AABDDAB";
//!
//! const MAPPING: ByteIndexMap<4> = ByteIndexMap::new(*b"ABCD", b'A');
//! const WEIGHTS: WeightMatrix<i8, 4> = WeightMatrix::new(&MAPPING, 1, -1, None);
//! const GAP_OPEN: i8 = -4;
//! const GAP_EXTEND: i8 = -2;
//!
//! let profile = LocalProfiles::<32, 4>::new_with_i8(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
//! let score = profile.smith_waterman_score_from_i8(reference).unwrap();
//! assert_eq!(score, 5);
//! ```
//!
//! [`sw_scalar_score`]: sw::sw_scalar_score
//! [`sw_simd_score`]: sw::sw_simd_score
//! [`sw_scalar_alignment`]: sw::sw_scalar_alignment
//! [`ByteIndexMap::new`]: crate::data::ByteIndexMap::new
//! [`WeightMatrix`]: crate::data::WeightMatrix
//! [`WeightMatrix::new`]: crate::data::WeightMatrix::new
//! [`into_biased_matrix`]: crate::data::WeightMatrix::into_biased_matrix
//! [`new_biased_dna_matrix`]: crate::data::WeightMatrix::new_biased_dna_matrix
//! [`new_dna_matrix`]: crate::data::WeightMatrix::new_dna_matrix
//! [`QueryProfileError`]: crate::data::err::QueryProfileError
//! [`Nucleotides::into_local_profile`]:
//!     crate::data::types::nucleotides::Nucleotides::into_local_profile
//! [`Nucleotides::into_shared_profile`]:
//!     crate::data::types::nucleotides::Nucleotides::into_shared_profile
//! [`ScalarProfile::new`]: ScalarProfile::new
//! [`StripedProfile::new`]: StripedProfile::new
//! [`ScalarProfile::smith_waterman_score`]: ScalarProfile::smith_waterman_score
//! [`StripedProfile::smith_waterman_score`]:
//!    StripedProfile::smith_waterman_score
//! [`LocalProfiles`]: LocalProfiles
//! [`SharedProfiles`]: SharedProfiles
//! [`LocalProfiles::new_with_i8`]: LocalProfiles::new_with_i8
//! [`SharedProfiles::new_with_i8`]: SharedProfiles::new_with_i8
//! [`LocalProfiles::smith_waterman_score_from_i8`]:
//!     LocalProfiles::smith_waterman_score_from_i8
//! [`SharedProfiles::smith_waterman_score_from_i8`]:
//!     SharedProfiles::smith_waterman_score_from_i8

pub mod sw;

mod profile;
mod profile_set;
mod state;
mod std_traits;

pub use profile::*;
pub use profile_set::*;
pub use state::*;
