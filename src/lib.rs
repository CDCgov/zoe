#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]
#![warn(clippy::all, clippy::pedantic)]
#![allow(
    clippy::module_name_repetitions,
    clippy::similar_names,
    clippy::wildcard_imports,
    clippy::enum_glob_use,
    clippy::comparison_chain,
    stable_features
)]
#![cfg_attr(test, feature(test))]
#![feature(portable_simd, const_fn_floating_point_arithmetic, let_chains, lazy_cell)]
#![forbid(incomplete_features)]

/// ## Functions for aligning sequence data.
///
/// *Zoe* supports efficient local alignment for DNA, protein, or any other
/// sequence data via the Smith-Waterman algorithm. For generating the score
/// (useful for database searches and determining whether sequences are
/// related), this includes both a traditional implementation (see
/// [`sw_scalar_score`]) and the usual striped SIMD implementation (see
/// [`sw_simd_score`]). The function [`sw_scalar_alignment`] can be used to
/// generate the alignment itself.
///
/// To calculate an alignment score or generate an alignment, follow these
/// steps. The first two are designed to be possible at compile time, and the
/// subsequent two to execute at runtime.
///
/// 1. Choose an alphabet. *Zoe* provides an easy interface to work with DNA,
///    assuming the bases `ACGT` are case-insensitive and `N` is used as a
///    catch-all. Otherwise, you can define your own alphabet using
///    [`ByteIndexMap::new`].
///
/// 2. Specify the matrix of weights used for scoring matches and mismatches.
///    SIMD algorithms require a [`BiasedWeightMatrix`], while scalar algorithms
///    use a [`SimpleWeightMatrix`]. For DNA, either can be quickly constructed
///    with [`new_biased_dna_matrix`] and [`new_dna_matrix`] respectively. For a
///    custom alphabet, use [`SimpleWeightMatrix::new`] and then convert it into
///    a [`BiasedWeightMatrix`] using [`into_biased_matrix`] if using SIMD.
///
/// 3. Build the query profile, with either [`ScalarProfile::new`] or
///    [`StripedProfile::new`] (for SIMD). This step combines the query, matrix
///    of weights, and gap open and gap extend penalties. This step also
///    performs some basic checks, and if any of the inputs are invalid, a
///    [`QueryProfileError`] is returned.
///
/// 4. Use the query profile to align against any number of different
///    references. Simply call the [`ScalarProfile::smith_waterman_score`] or
///    [`StripedProfile::smith_waterman_score`] methods. The profile can be
///    reused in multiple different calls.
///
/// Using DNA:
/// ```
/// # use zoe::{alignment::StripedProfile, data::BiasedWeightMatrix};
/// let reference: &[u8] = b"GGCCACAGGATTGAG";
/// let query: &[u8] = b"CTCAGATTG";
///
/// const WEIGHTS: BiasedWeightMatrix<5> = BiasedWeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
///
/// let profile = StripedProfile::<u8, 32, 5>::new(query, &WEIGHTS, 3, 1).unwrap();
/// let score = profile.smith_waterman_score(reference).unwrap();
/// assert_eq!(score, 27);
/// ```
///
/// Using a different alphabet:
/// ```
/// # use zoe::{
/// #     alignment::StripedProfile,
/// #     data::{BiasedWeightMatrix, ByteIndexMap, SimpleWeightMatrix},
/// # };
/// let reference: &[u8] = b"BDAACAABDDDB";
/// let query: &[u8] = b"AABDDAB";
///
/// const MAPPING: ByteIndexMap<4> = ByteIndexMap::new(*b"ABCD", b'A');
/// const WEIGHTS: BiasedWeightMatrix<4> = SimpleWeightMatrix::new(&MAPPING, 1, -1, None).into_biased_matrix();
///
/// let profile = StripedProfile::<u8, 32, 4>::new(query, &WEIGHTS, 4, 2).unwrap();
/// let score = profile.smith_waterman_score(reference).unwrap();
/// assert_eq!(score, 5);
/// ```
///
/// When using the SIMD algorithm, you must specify the number of lanes `N` and
/// integer type `T`. If the integer type has too small of a range, it is
/// possible the alignment will overflow (in which case `None` is returned). A
/// higher-level abstraction to avoid this is [`LocalProfile`] and
/// [`SharedProfile`].
///
/// The former is designed for use within a single thread,
/// while the latter allows multiple threads to access it. Both store a set of
/// lazily-evaluated query profiles for `u8`, `u16`, `u32`, and `u64`. To create
/// one of these, you can call [`LocalProfile::new_with_u8`],
/// [`SharedProfile::new_with_u8`], or one of the other constructors. Then, call
/// [`LocalProfile::smith_waterman_score_from_u8`],
/// [`SharedProfile::smith_waterman_score_from_u8`], or one of the other
/// methods.
///
/// When using DNA, you can also create a profile by using
/// [`Nucleotides::into_local_profile`] or [`Nucleotides::into_shared_profile`].
///
/// Using DNA:
/// ```
/// # use zoe::{data::BiasedWeightMatrix, prelude::Nucleotides};
/// let reference: &[u8] = b"GGCCACAGGATTGAG";
/// let query: Nucleotides = b"CTCAGATTG".into();
///
/// const WEIGHTS: BiasedWeightMatrix<5> = BiasedWeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
///
/// let profile = query.into_local_profile::<32, 5>(&WEIGHTS, 3, 1).unwrap();
/// let score = profile.smith_waterman_score_from_u8(reference).unwrap();
/// assert_eq!(score, 27);
/// ```
///
/// Using a different alphabet:
/// ```
/// # use zoe::{
/// #     alignment::LocalProfile,
/// #     data::{BiasedWeightMatrix, ByteIndexMap, SimpleWeightMatrix},
/// # };
/// let reference: &[u8] = b"BDAACAABDDDB";
/// let query: &[u8] = b"AABDDAB";
///
/// const MAPPING: ByteIndexMap<4> = ByteIndexMap::new(*b"ABCD", b'A');
/// const WEIGHTS: BiasedWeightMatrix<4> = SimpleWeightMatrix::new(&MAPPING, 1, -1, None).into_biased_matrix();
///
/// let profile = LocalProfile::<32, 4>::new_with_u8(query, &WEIGHTS, 4, 2).unwrap();
/// let score = profile.smith_waterman_score_from_u8(reference).unwrap();
/// assert_eq!(score, 5);
/// ```
///
/// [`sw_scalar_score`]: alignment::sw::sw_scalar_score
/// [`sw_simd_score`]: alignment::sw::sw_simd_score
/// [`sw_scalar_alignment`]: alignment::sw::sw_scalar_alignment
/// [`ByteIndexMap::new`]: data::ByteIndexMap::new
/// [`BiasedWeightMatrix`]: data::BiasedWeightMatrix
/// [`SimpleWeightMatrix`]: data::SimpleWeightMatrix
/// [`SimpleWeightMatrix::new`]: data::SimpleWeightMatrix::new
/// [`into_biased_matrix`]: data::SimpleWeightMatrix::into_biased_matrix
/// [`new_biased_dna_matrix`]: data::BiasedWeightMatrix::new_biased_dna_matrix
/// [`new_dna_matrix`]: data::SimpleWeightMatrix::new_dna_matrix
/// [`QueryProfileError`]: data::err::QueryProfileError
/// [`Nucleotides::into_local_profile`]: data::types::nucleotides::Nucleotides::into_local_profile
/// [`Nucleotides::into_shared_profile`]: data::types::nucleotides::Nucleotides::into_shared_profile
/// [`ScalarProfile::new`]: alignment::ScalarProfile::new
/// [`StripedProfile::new`]: alignment::StripedProfile::new
/// [`ScalarProfile::smith_waterman_score`]: alignment::ScalarProfile::smith_waterman_score
/// [`StripedProfile::smith_waterman_score`]: alignment::StripedProfile::smith_waterman_score
/// [`LocalProfile`]: alignment::LocalProfile
/// [`SharedProfile`]: alignment::SharedProfile
/// [`LocalProfile::new_with_u8`]: alignment::LocalProfile::new_with_u8
/// [`SharedProfile::new_with_u8`]: alignment::SharedProfile::new_with_u8
/// [`LocalProfile::smith_waterman_score_from_u8`]: alignment::LocalProfile::smith_waterman_score_from_u8
/// [`SharedProfile::smith_waterman_score_from_u8`]: alignment::SharedProfile::smith_waterman_score_from_u8
pub mod alignment;
/// Composition and consensus functions.
pub mod composition;
/// Data import, export, and manipulation functions.
pub mod data;
/// Distance functions, especially for sequence data.
pub mod distance;
/// K-mer processing, searching, and counting.
pub mod kmer;
/// Sequence search and/or replacement.
pub mod search;

#[cfg(feature = "rand")]
/// Generate sequences and other data.
pub(crate) mod generate;
/// Iterator utilities.
pub(crate) mod iter_utils;
/// Mathematical utilities.
pub(crate) mod math;
/// SIMD traits to extend portable SIMD.
pub(crate) mod simd;

/// Common structures and traits re-exported
pub mod prelude {
    pub use crate::composition::{AlignmentComposition, CreateConsensus, NucleotideCounts};
    pub use crate::data::types::{amino_acids::AminoAcids, nucleotides::Nucleotides};
    pub use crate::data::{
        convert::ToDNA, err::OrFail, fasta::FastaReader, fastq::FastQReader, PairwiseSequence, StdForSequences,
    };
    #[cfg(feature = "rand")]
    pub use crate::generate::rand_sequence;
    pub use crate::kmer::{encoder::KmerEncoder, kmer_counter::KmerCounter, kmer_set::KmerSet};
    pub use crate::search::{ByteSubstring, ByteSubstringMut};
}
