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
/// For example:
/// ```
/// # use zoe::prelude::*;
/// let owned_sequence: Nucleotides = b"GGCCACCAAGGCCA".into();
///
/// let view_of_sequence = owned_sequence.as_view();
/// assert_eq!(view_of_sequence.as_bytes(), b"GGCCACCAAGGCCA");
/// assert_eq!(view_of_sequence.gc_content(), 10);
///
/// let slice_of_middle = owned_sequence.slice(4..9);
/// assert_eq!(slice_of_middle.as_bytes(), b"ACCAA");
/// assert_eq!(slice_of_middle.gc_content(), 2);
///
/// let slice_of_middle_again = slice_of_middle.slice(2..);
/// assert_eq!(slice_of_middle_again.as_bytes(), b"CAA");
///
/// let view_to_owned = slice_of_middle.to_owned();
/// assert_eq!(view_to_owned, b"ACCAA".into());
///
/// let mut owned_sequence: Nucleotides = b"GGCCACCAAGGCCA".into();
/// let mut mutable_slice = owned_sequence.slice_mut(1..4);
/// mutable_slice.as_mut_bytes().fill(b'T');
/// assert_eq!(owned_sequence.as_bytes(), b"GTTTACCAAGGCCA");
/// ```
///
/// Views can also be re-sliced in-place using the [`restrict`] method. This is
/// a useful way to avoid needing to have extra let bindings. For example:
/// ```
/// # use zoe::prelude::*;
/// let mut owned_sequence: Nucleotides = b"GGCCACCAAGGCCA".into();
/// let mut mutable_slice = owned_sequence.as_view_mut();
/// mutable_slice.restrict(5..);
/// mutable_slice.restrict(..2);
/// assert_eq!(mutable_slice.as_bytes(), b"CC");
/// ```
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
/// ## Distance functions for sequence data and byte strings.
///
/// This module provides distance functions for DNA in the [`dna`] module, amino
/// acids in the [`aa`] module, and byte substrings (such as [`hamming`] and
/// [`hamming_simd`]).
///
/// Within both [`dna`] and [`aa`], there are standalone functions accepting
/// byte slices, and there are trait methods [`NucleotidesDistance`] and
/// [`AminoAcidsDistance`].
///
/// [`dna`]: distance::dna
/// [`aa`]: distance::aa
/// [`hamming`]: distance::hamming
/// [`hamming_simd`]: distance::hamming_simd
/// [`NucleotidesDistance`]: distance::dna::NucleotidesDistance
/// [`AminoAcidsDistance`]: distance::aa::AminoAcidsDistance
pub mod distance;

/// ## K-mer processing, searching, and counting.
///
/// This module provides structs and methods for handling common k-mer
/// operations efficiently using integer encodings. A k-mer is a short
/// subsequence of nucleotides, which are represented in *Zoe* by [`Kmer`].
///
/// *Zoe* has two structs to store k-mers:
/// * A [`KmerSet`], which is a [`HashSet`] for k-mers
/// * A [`KmerCounter`] to store k-mers and their counts
///
/// Several convenience methods are provided, such as:
/// * [`insert_from_sequence`], which quickly inserts/counts all overlapping
///   k-mers from a sequence
/// * [`find_in_seq`] and [`find_in_seq_rev`], which search a sequence for any
///   of the stored k-mers. The [`FindKmers`] trait provides similar methods
///   [`find_kmers`] and [`find_kmers_rev`] directly on the sequence, and
///   accepting the kmer collection as an argument.
///
/// Many structs related to k-mers, including [`KmerSet`] and [`KmerCounter`],
/// are generic over a maximum possible k-mer length `MAX_LEN`. This is used to
/// determine (at compile-time) the appropriate integer type to store the
/// encoded k-mer. The actual k-mer length can be set to a different value at
/// runtime as long as it is less than `MAX_LEN`. For guidance on picking the
/// appropriate `MAX_LEN`, see [`SupportedKmerLen`].
///
/// [`KmerSet`] and [`KmerCounter`] are also generic over the type of
/// [`KmerEncoder`]. *Zoe* currently provides a single encoder,
/// [`ThreeBitKmerEncoder`], which uses three bits to store each base. It allows
/// for `A`, `C`, `G`, `T`, and `N` to all be represented. It does not preserve
/// case or the distinction between `T` and `U`. `N` is used as a catch-all for
/// bases that are not `ACGTUNacgtun`.
///
/// For convenience, the type aliases [`ThreeBitKmerSet`] and
/// [`ThreeBitKmerCounter`] are provided.
///
/// <div class="warning important">
///
/// **Important**
///
/// When performing more specialized k-mer operations, you may need to directly
/// encode and decode k-mers, rather than relying on the methods in [`KmerSet`]
/// or [`KmerCounter`]. It is important to use the same [`KmerEncoder`] to both
/// encode and decode the k-mer.
///
/// </div>
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
/// # use zoe::{kmer::ThreeBitKmerSet, prelude::*};
/// let primer = b"TGATAGTTTTAGAGTTAGGTAG";
/// let sequence = b"TGCCCGTAACGTACAGTTTTACAGTTAGGTACCC";
/// let mut kmer_set = ThreeBitKmerSet::<17>::new(17).unwrap();
/// kmer_set.insert_from_sequence_with_variants::<1>(primer);
/// let kmer_pos = kmer_set.find_in_seq(sequence);
/// assert_eq!(kmer_pos, Some(14..31));
/// ```
///
/// This can be equivalent performed using:
/// ```
/// # use zoe::{kmer::{ThreeBitKmerSet, FindKmers}, prelude::*};
/// let primer: Nucleotides = b"TGATAGTTTTAGAGTTAGGTAG".into();
/// let sequence: Nucleotides = b"TGCCCGTAACGTACAGTTTTACAGTTAGGTACCC".into();
/// let mut kmer_set = ThreeBitKmerSet::<17>::new(17).unwrap();
/// kmer_set.insert_from_sequence_with_variants::<1>(primer);
/// let kmer_pos = sequence.find_kmers(&kmer_set);
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
/// [`ThreeBitKmerCounter`]: kmer::ThreeBitKmerCounter
/// [`SupportedKmerLen`]: kmer::SupportedKmerLen
/// [`insert_from_sequence`]: kmer::KmerSet::insert_from_sequence
/// [`find_in_seq`]: kmer::KmerCollectionContains::find_in_seq
/// [`find_in_seq_rev`]: kmer::KmerCollectionContains::find_in_seq_rev
/// [`FindKmers`]: kmer::FindKmers
/// [`find_kmers`]: kmer::FindKmers::find_kmers
/// [`find_kmers_rev`]: kmer::FindKmers::find_kmers_rev
pub mod kmer;
/// Mathematical utilities.
pub mod math;
/// ## Sequence search and/or replacement.
///
/// This module provides functions for performing searches and/or replacements
/// within byte slices or vectors. In addition to many customized search
/// functions, two general-purpose search algorithms are [`substring_match`] and
/// [`substring_match_simd`]. The latter is more efficient for larger needles.
///
/// Similarly, for fuzzy matching (permitting some number of mismatches between
/// the desired and located needle), [`fuzzy_substring_match`] and
/// [`fuzzy_substring_match_simd`] are provided.
///
/// The [`ByteSubstring`] trait provides methods which wrap these underlying
/// algorithms, such as [`contains_substring`], [`find_substring`], and
/// [`find_fuzzy_substring`]. These methods return a range instead of just the
/// starting index.
///
/// ## Restricting the search range
///
/// Sometimes, you may want to search a particular region of a string. There are
/// two main scenarios:
/// 1. The string search should return an index/range with respect to the
///    subsequence. For example, searching `ACGT` for `G` in the range `1..3`
///    would return a starting index of 1.
/// 2. The string search should return an index/range with respect to the full,
///    original sequence. Using the same example, the starting index would be 2.
///
/// The first scenario is easily resolved by creating a byte slice or a
/// [view](crate::data#views), and then searching within that. For example:
/// ```
/// # use zoe::prelude::*;
/// let haystack = b"CACATAACGTACAGTTTTACACAGTTAGGT";
/// let needle = b"CACA";
/// let position = haystack[10..25].find_substring(needle);
/// assert_eq!(position, Some(9..13))
/// ```
///
/// The second scenario could also be solved the same way, but it would require
/// shifting the output by the starting index. An alternative is the
/// [`RangeSearch`] struct, which implements [`ByteSubstring`] and automatically
/// adjusts the range afterwards. An example usage is below. Notice that the
/// potential match at the beginning of the haystack is not returned, since it
/// is not in the provided range.
/// ```
/// # use zoe::{search::ToRangeSearch, prelude::*};
/// let haystack = b"CACATAACGTACAGTTTTACACAGTTAGGT";
/// let needle = b"CACA";
/// let position = haystack.search_in(10..25).find_substring(needle);
/// assert_eq!(position, Some(19..23))
/// ```
///
/// [`RangeSearch`] also provides additional convenience methods. One use-case
/// is when a needle must be found in at most the first `n` elements of the
/// haystack, or at most the last `n` elements of the haystack. Ordinarily, one
/// would have to ensure the range is in bounds before slicing or using
/// `search_in`. To simplify this, *Zoe* provides `search_in_first` and
/// `search_in_last`, as shown below. Once again, the second occurrence of the
/// needle is returned.
/// ```
/// # use zoe::{search::ToRangeSearch, prelude::*};
/// let haystack = b"CACATAACGTACAGTTTTACACAGTTAGGT";
/// let needle = b"CACA";
/// let position = haystack.search_in_last(15).find_substring(needle);
/// assert_eq!(position, Some(19..23))
/// ```
///
/// [`substring_match`]: crate::search::substring_match
/// [`substring_match_simd`]: crate::search::substring_match_simd
/// [`fuzzy_substring_match`]: crate::search::fuzzy_substring_match
/// [`fuzzy_substring_match_simd`]: crate::search::fuzzy_substring_match_simd
/// [`ByteSubstring`]: crate::search::ByteSubstring
/// [`contains_substring`]: crate::search::ByteSubstring::contains_substring
/// [`find_substring`]: crate::search::ByteSubstring::find_substring
/// [`find_fuzzy_substring`]: crate::search::ByteSubstring::find_fuzzy_substring
/// [`RangeSearch`]: crate::search::RangeSearch
pub mod search;

#[cfg(feature = "rand")]
/// Generate sequences and other data.
pub(crate) mod generate;
/// Iterator utilities.
pub(crate) mod iter_utils;

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
    pub use crate::kmer::{EncodedKmerCollection, FindKmers, KmerCollectionContains, KmerCounter, KmerEncoder, KmerSet};
    pub use crate::search::{ByteSubstring, ByteSubstringMut};
}

pub(crate) use crate::data::extension::simd;
