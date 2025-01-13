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
/// Below is an example using DNA. We use a match score of 4 and a mismatch
/// score of -2, as defined in `WEIGHTS`. The last two arguments to
/// [`StripedProfile::new`] specify the gap open penalty as 3 and the gap extend
/// penalty as 1.
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
/// Below is an example using a different alphabet. Matches are given a score of
/// 1, mismatches are given a score of -1, the gap open penalty is 4, and the
/// gap extend penalty is 2.
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
/// higher-level abstraction to avoid this is [`LocalProfiles`] and
/// [`SharedProfiles`].
///
/// The former is designed for use within a single thread, while the latter
/// allows multiple threads to access it. Both store a set of lazily-evaluated
/// query profiles for `u8`, `u16`, `u32`, and `u64`. To create one of these,
/// you can call [`LocalProfiles::new_with_u8`], [`SharedProfiles::new_with_u8`],
/// or one of the other constructors. Then, call
/// [`LocalProfiles::smith_waterman_score_from_u8`],
/// [`SharedProfiles::smith_waterman_score_from_u8`], or one of the other
/// methods.
///
/// When using DNA, you can also create a profile by using
/// [`Nucleotides::into_local_profile`] or [`Nucleotides::into_shared_profile`].
///
/// Below is the previous example using DNA, but with this higher-level API:
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
/// And similarly with the different alphabet:
/// ```
/// # use zoe::{
/// #     alignment::LocalProfiles,
/// #     data::{BiasedWeightMatrix, ByteIndexMap, SimpleWeightMatrix},
/// # };
/// let reference: &[u8] = b"BDAACAABDDDB";
/// let query: &[u8] = b"AABDDAB";
///
/// const MAPPING: ByteIndexMap<4> = ByteIndexMap::new(*b"ABCD", b'A');
/// const WEIGHTS: BiasedWeightMatrix<4> = SimpleWeightMatrix::new(&MAPPING, 1, -1, None).into_biased_matrix();
///
/// let profile = LocalProfiles::<32, 4>::new_with_u8(query, &WEIGHTS, 4, 2).unwrap();
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
/// [`Nucleotides::into_local_profile`]:
///     data::types::nucleotides::Nucleotides::into_local_profile
/// [`Nucleotides::into_shared_profile`]:
///     data::types::nucleotides::Nucleotides::into_shared_profile
/// [`ScalarProfile::new`]: alignment::ScalarProfile::new
/// [`StripedProfile::new`]: alignment::StripedProfile::new
/// [`ScalarProfile::smith_waterman_score`]:
///     alignment::ScalarProfile::smith_waterman_score
/// [`StripedProfile::smith_waterman_score`]:
///     alignment::StripedProfile::smith_waterman_score
/// [`LocalProfiles`]: alignment::LocalProfiles
/// [`SharedProfiles`]: alignment::SharedProfiles
/// [`LocalProfiles::new_with_u8`]: alignment::LocalProfiles::new_with_u8
/// [`SharedProfiles::new_with_u8`]: alignment::SharedProfiles::new_with_u8
/// [`LocalProfiles::smith_waterman_score_from_u8`]:
///     alignment::LocalProfiles::smith_waterman_score_from_u8
/// [`SharedProfiles::smith_waterman_score_from_u8`]:
///     alignment::SharedProfiles::smith_waterman_score_from_u8
pub mod alignment;
/// Composition and consensus functions.
pub mod composition;
// TODO: Add examples, expand docs, fix links
/// ## Data import, export, and manipulation functions.
///
/// ## Views
///
/// Many of the data type provided by *Zoe* have versions holding owned data as
/// well as versions holding references. The latter are called *views*, and are
/// useful when performing operations on a subsequence of the original data. For
/// example, [`Nucleotides`] (which is a wrapper around [`Vec<u8>`]) has the
/// corresponding types [`NucleotidesView`] and [`NucleotidesViewMut`], which
/// are wrappers around immutable and mutable byte slices respectively.
///
/// A view can be constructed directly from an existing slice, such as
/// [`NucleotidesView::from_bytes_unchecked`], but more often a view is created
/// from an owned instance of the data by calling [`as_view`] or
/// [`as_view_mut`]. If only a range of the data is desired, then [`slice`] and
/// [`slice_mut`] are used.
///
/// A view can be directly displayed/printed, or it can be copied into owned
/// data again by using [`to_owned_data`].
///
/// Views can also be re-sliced in-place using the [`restrict`] method. This is
/// a useful way to avoid needing to have extra let bindings.
///
/// [`Nucleotides`]: data::types::nucleotides::Nucleotides
/// [`NucleotidesView`]: data::types::nucleotides::NucleotidesView
/// [`NucleotidesViewMut`]: data::types::nucleotides::NucleotidesViewMut
/// [`NucleotidesView::from_bytes_unchecked`]:
///     data::types::nucleotides::NucleotidesView::from_bytes_unchecked
/// [`as_view`]: prelude::DataOwned::as_view
/// [`as_view_mut`]: prelude::DataOwned::as_view_mut
/// [`slice`]: prelude::Slice::slice
/// [`slice_mut`]: prelude::SliceMut::slice_mut
/// [`to_owned_data`]: prelude::DataView::to_owned_data
/// [`restrict`]: prelude::DataView::restrict
pub mod data;
/// Distance functions, especially for sequence data.
pub mod distance;

/// ## K-mer processing, searching, and counting.
///
/// This module provides structs and methods for handling common k-mer
/// operations efficiently using integer encodings. A k-mer is a short
/// subsequence of nucleotides, which are represented in *Zoe* by [`Kmer`].
/// Although the length of the k-mer might not be known at compile time, a
/// maximum possible length (`MAX_LEN`) is required, which is used to determine
/// the appropriate integer type to store the encoded k-mer.
///
/// *Zoe* has two types of structs to store k-mers:
/// * A [`KmerSet`], which is a [`HashSet`] for k-mers
/// * A [`KmerCounter`] to store k-mers and their counts
///
/// Several convenience methods are provided, such as:
/// * [`insert_from_sequence`], which quickly inserts/counts all overlapping
///   k-mers from a sequence
/// * [`find_kmers`] and [`find_kmers_rev`], which search a sequence for any of
///   the stored k-mers
///
/// Any particular implementation of [`KmerSet`] or [`KmerCounter`] relies on an
/// encoder. *Zoe* currently provides a single encoder, [`ThreeBitKmerEncoder`],
/// which uses three bits to store each base. It allows for `A`, `C`, `G`, `T`,
/// and `N` to all be represented. It does not preserve case or the distinction
/// between `T` and `U`. `N` is used as a catch-all for bases that are not
/// `ACGTUNacgtun`.
///
/// The corresponding types with this encoder are [`ThreeBitKmerSet`] and
/// [`ThreeBitKmerCounter`]. A more specialized [`ThreeBitOneMismatchKmerSet`]
/// is also provided, which automatically stores all k-mers within a Hamming
/// distance of 1 upon insertion.
///
/// For guidance on picking the appropriate `MAX_LEN`, see
/// [`SupportedThreeBitKmerLen`].
///
/// When performing more specialized k-mer operations, you may need to directly
/// encode and decode k-mers, rather than relying on the methods in [`KmerSet`]
/// or [`KmerCounter`]. It is important to use the same [`KmerEncoder`] to both
/// encode and decode the k-mer. For [`ThreeBitEncodedKmer`], there also exists
/// a [`Display`] implementation that does not require access to the encoder,
/// but this is suboptimal; it is preferred to decode the k-mer first.
///
/// ## Examples
///
/// Count the 3-mers in a sequence:
/// ```
/// # use zoe::{kmer::ThreeBitKmerCounter, prelude::*};
/// let sequence = b"GGCCACCAAGGCCA";
/// let mut kmer_counter = ThreeBitKmerCounter::<3>::new(3).unwrap();
/// kmer_counter.insert_from_sequence(sequence);
/// for (kmer, count) in kmer_counter {
///     println!("{kmer}\t{count}");
/// }
/// ```
///
/// Search for the 17-mers of a primer within a sequence, with up to one
/// mismatch:
/// ```
/// # use zoe::{kmer::ThreeBitOneMismatchKmerSet, prelude::*};
/// let primer = b"TGATAGTTTTAGAGTTAGGTAG";
/// let sequence = b"TGCCCGTAACGTACAGTTTTACAGTTAGGTACCC";
/// let mut kmer_set = ThreeBitOneMismatchKmerSet::<17>::new(17).unwrap();
/// kmer_set.insert_from_sequence(primer);
/// let kmer_pos = kmer_set.find_kmers(sequence);
/// assert_eq!(kmer_pos, Some(14..31));
/// ```
///
/// [`HashSet`]: std::collections::HashSet
/// [`Display`]: std::fmt::Display
/// [`Kmer`]: kmer::Kmer
/// [`KmerEncoder`]: kmer::KmerEncoder
/// [`KmerSet`]: kmer::KmerSet
/// [`KmerCounter`]: kmer::KmerCounter
/// [`ThreeBitKmerEncoder`]: kmer::ThreeBitKmerEncoder
/// [`ThreeBitEncodedKmer`]: kmer::ThreeBitEncodedKmer
/// [`ThreeBitKmerSet`]: kmer::ThreeBitKmerSet
/// [`ThreeBitOneMismatchKmerSet`]: kmer::ThreeBitOneMismatchKmerSet
/// [`ThreeBitKmerCounter`]: kmer::ThreeBitKmerCounter
/// [`SupportedThreeBitKmerLen`]: kmer::three_bit::SupportedThreeBitKmerLen
/// [`insert_from_sequence`]: kmer::KmerSet::insert_from_sequence
/// [`find_kmers`]: kmer::KmerSet::find_kmers
/// [`find_kmers_rev`]: kmer::KmerSet::find_kmers_rev
pub mod kmer;
/// Mathematical utilities.
pub mod math;
/// Sequence search and/or replacement.
pub mod search;

#[cfg(feature = "rand")]
/// Generate sequences and other data.
pub(crate) mod generate;
/// Iterator utilities.
pub(crate) mod iter_utils;

/// SIMD traits to extend portable SIMD.
pub(crate) mod simd {
    pub(crate) use crate::data::extension::simd::*;
}

/// Common structures and traits re-exported
pub mod prelude {
    pub use crate::alignment::PairwiseSequence;
    pub use crate::composition::{AlignmentComposition, CreateConsensus, GcContent, NucleotideCounts, ToBaseCounts};
    pub use crate::data::{
        StdForSequences,
        err::OrFail,
        records::{
            fasta::FastaReader,
            fastq::{FastQ, FastQReader, FastQView, FastQViewMut},
        },
        types::{
            amino_acids::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
            nucleotides::{
                Nucleotides, NucleotidesView, NucleotidesViewMut, RecodeNucleotides, RetainNucleotides, Translate,
            },
            phred::{QualityScores, QualityScoresView, QualityScoresViewMut, QualityStats},
        },
        view_traits::{DataOwned, DataView, DataViewMut, Len, Slice, SliceMut},
    };
    #[cfg(feature = "rand")]
    pub use crate::generate::rand_sequence;
    pub use crate::kmer::{encoder::KmerEncoder, kmer_counter::KmerCounter, kmer_set::KmerSet};
    pub use crate::search::{ByteSubstring, ByteSubstringMut};
}
