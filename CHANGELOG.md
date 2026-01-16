# Zoe Changelog

All notable changes to this project will be documented in this file. The format
is roughly based on [Keep a Changelog], and this project tries to adheres to
[Semantic Versioning].

## [0.0.24] - TBD

### Added

- Added `to_alignment` to convert `SamData` to `Alignment`
- Added ability to push tags to `SamTags`
- Added `prepend_state`, `prepend_ciglet`, and `prepend_inc_op` to `AlignmentStates`
- Added `StatesSequenceMut` trait and `AlignmentStates::as_mut_slice` to provide mutable access to ciglets while iterating/processing
- Added `smith_waterman_alignment_from_i16_3pass` and `smith_waterman_alignment_from_i32_3pass` (behind `dev-3pass` feature gate)
- Added `get_first_codon`, `get_last_codon`, and `get_tail_codon` to `GetCodons` trait, with mutable equivalents in `GetCodonsMut`
- Added `CodonExtension` trait for `[u8; 3]` with convenience codon characterization methods
- Added `ToCigletIterator` trait, implemented on structs which can be converted to an iterator of ciglets
- Added `LenInAlignment` trait, replacing the previous `match_len` and offering a more flexible API
- Added `VARIANTS` constants to pHMM state enums (behind `dev-phmm` feature gate)
- Added index comparison methods to `PhmmIndex` (behind `dev-phmm` feature gate)
- Condensed the `PhmmParamKind` variants in the `visit_params` API (behind `dev-phmm` feature gate)
- Added row and column accessors to `TransitionParams` (behind `dev-phmm` and `alignment-diagnostics` feature gates)
- Added `WithErrorContext` and `ResultWithErrorContext` traits using the `ErrorWithContext` struct for updated error-handling.
- Added methods and implementations to iterate over CIGAR views

### Changed

- Modified alignment methods on profiles to take a `SeqSrc` enum as an argument, allowing the returned alignment to properly reflect which sequence is the query and which is the reference
- `SharedProfiles` no longer holds the sequence as `Box<u8>`, instead using `&'a [u8]`
- The `arbitrary` module now uses specification structs to allow various restrictions to be enforced when generating *Zoe* data types. See the module level documentation for more details.
- `CheckedCigar` no longer has a blanket implementation, and is instead manually implemented on relevant types (behind `fuzzing` feature gate)
- Updated error handling logic and style for `OrFail` trait.
- Standalone function `zoe::data::err::open_nonempty_file` is now public
- Alignment methods returning indices or ranges now return specialized structs
- `smith_waterman_alignment_3pass` is now private on `StripedProfile` (behind `dev-3pass` feature gate)

### Fixes

- Fixed spelling of `filter_to_dna_uanligned` to `filter_to_dna_unaligned` in trait `ToDNA`
- Fixed bug for parsing -128 in a weight matrix using `MatParser`

### Removed

- Removed the `match_len` method from CIGAR strings in favor of `LenInAlignment`
- Removed `CheckedCigar` in favor of `LenInAlignment` (checked methods are behind the `fuzzing` feature gate)

## [0.0.23] - 2025-11-24

### Added

- Added `nw_scalar_score` and `nw_scalar_alignment` for performing global sequence alignment
- Added random downsampling for iterators of known and unknown lengths in `iter_utils::sampling`
- Added the method `process_results` to more robustly work with iterators of results
- Added `get_aligned_query`, similar to `get_aligned_seqs` but returning only the query sequence
- Added additional getters and convenience methods on `SamTags` and `SamTagValue`
- Added `to_eq_x` to `Alignment` and `AlignmentStates`, which converts `M` to either `=` or `X` in the alignment
- Added `SequenceReadable` getter trait to record types and sequence types
- Added `from_ciglets_unchecked` and `from_cigar_unchecked` to `AlignmentStates`
- Added `decode_iter` to `KmerEncoder`
- Added `nw_score_from_path` for rescoring a global alignment (behind `alignment-diagnostics` feature gate)
- Introduced the `AlignmentAndSeqs` and `AminoAcidsIupacX` arbitrary wrappers (behind `fuzzing` feature gate)
- Added `visit_params` methods to pHMMs as a diagnostic tool for inspecting the parameters that are visited along a given path (behind `dev-phmm` and `alignment-diagnostics` feature gates)
- Added the ability to remove layers from a pHMM (behind `dev-phmm` and `alignment-diagnostics` feature gates)

### Changed

- `TryFrom` for `AlignmentStates` now merges adjacent ciglets with the same operation
- An iterator of ciglets now collects into a `Result<AlignmentStates, CigarError>` instead of `AlignmentStates`, and performs merging as above
- The view traits have been refactored so that all associated types are defined in the `ViewAssocTypes` trait
- `reborrow_view` has been added to allow invariant lifetimes to be shortened in views
- `Vec<u8>`, `&[u8]`, and `&mut [u8]` now implement the various view traits
- Pairs of methods accepting encoded and decoded k-mers have been merged
- `KmerCounter` has had its methods renamed from `insert` to `tally`
- Iterators over encoded k-mers from `KmerSet` now return owned values
- `KmerCollectionContains` is renamed to `FindKmersInSeq`
- The `kmer` module has been restructured to simplify the docs
- Implemented `std::fmt::Binary` for `ThreeBitEncodedKmer`, and added it as a trait bound to `AnyInt`
- `sw_score_from_path` and `sw_scalar_alignment_override` are now behind the `alignment-diagnostics` feature gate
- Modified the behavior of `AlignmentArbitrary` and `AlignmentStatesArbitrary` (behind `fuzzing` feature gate)
- Refactored the organization/visibility of pHMM code; modified `SamHmmParser` constructors; applied compact bit representation for Viterbi traceback; switched to eager evaluation for `next_index` and `prev_index` (behind `dev-phmm` feature gate)
- The inner fields of `SAMReader` and `FastQReader` are now properly private, and `SliceRange` is properly public

### Fixed

- The `fuzzing` feature now successfully compiles even when `dev-phmm` is not enabled
- Fixed another NaN error that would occur when certain transitions in a pHMM were zero probability (behind `dev-phmm` feature gate)

### Removed

- Removed checked methods from `KmerEncoder` (they only provided partial checking and proved not useful)

## [0.0.22] - 2025-09-25

### Added

- Added banded Smith Waterman implementation as a more memory-efficient alternative for long sequence alignment
- Added `HeaderReadable` getter trait for record types
- Added `slice_to_ref_range` method for `Alignment`, in order to slice/clip an alignment to a desired portion of the reference
- Added `CheckedCigar` trait for summing the lengths for operations consuming the query/reference (behind `fuzzing` feature fate)
- Implemented `Arbitrary` for `AlignmentStates` and `Alignment` along with additional wrapper types (behind `fuzzing` feature gate)
- Added Viterbi algorithm for `DomainPhmm` (behind `dev-phmm` feature gate)
- Implemented `Arbitrary` for `SemiLocalPhmm` and `DomainPhmm` (behind `fuzzing` and `dev-phmm` feature gates)
- Added `sw_scalar_alignment_override` for more flexible testing of dynamic programming alignment algorithms (behind `fuzzing` feature gate)

### Changed

- The wrapper types for generating arbitrary CIGAR strings have been modified (behind `fuzzing` feature gate)
- pHMM types and algorithms now rely on the more general `PhmmNumber` trait instead of `Float` (behind `dev-phmm` feature gate)

### Fixes

- `ProfileSets` is re-exported in the prelude
- Fixed the initialization of the first match score for the semilocal Viterbi algorithm (behind `dev-phmm` feature gate)
- Fixed the calculation of the end score for semilocal Viterbi algorithm to allow passing through the END state (behind `dev-phmm` feature gate)
- Avoid error with NaN that may appear when using the `score_from_path` function on local or domain pHMMs (behind `dev-phmm` feature gate)

## [0.0.21] - 2025-09-05

### Added

- Added `sw_simd_score_ends` for getting the score _and end coordinates_ of the alignment
- Added byte index map for amino acids, BLOSUM scoring matrices, and PAM scoring matrices
- Added display implementation for `WeightMatrix`
- Added `get_subset` to `WeightMatrix` for subsetting a weight matrix with a different alphabet
- Added `parse_matrix` for parsing a `WeightMatrix` from a file at compile time
- Added accessor `byte_keys()` for `ByteIndexMap`
- Added unchecked constructors for `ScalarProfile` and `StripedProfile`
- Added `SamDataView` and `SamDataViewMut`
- Added `map` attribute syntax to `define_whichever` and `impl_traits`, so that a map can be applied when implementing iterator
- Added `dev-3pass` feature gate for additional alignment functions that return coordinates, as well as a three-pass alignment function with reduced memory usage
- Added `score_from_path` and `viterbi` for `LocalPhmm` (behind `dev-phmm` feature gate)
- Added `DomainPhmm`, along with `score_from_path` for it (behind `dev-phmm` feature gate)
- Added `SemiLocalPhmm`, along with `score_from_path` and `viterbi` for it (behind `dev-phmm` feature gate)

### Changed

- `MaybeAligned` now can hold other types besides `Alignment`, and score-only functions now return `MaybeAligned<u32>`
- `into_biased_matrix` has been renamed to `to_biased_matrix` and takes `&self`
- `WeightMatrix` is now imported under `data::matrices` in addition to `data`
- `WeightMatrix` now has an additional lifetime generic for the `ByteIndexMap`
- `WeightMatrix` no longer implements `Copy`
- The `@` symbol is no longer included at the beginning of a header in a `FastQ` record (the reader strips it, and the display implementation re-adds it)
- Added `tags` field to `SamData`, as well as the ability to opt out of storing the tags in `SamReader`
- Changed `Debug` impl of `CigarViewMut` and derived more traits
- Various Smith-Waterman profile set methods have been moved to trait methods
  for `ProfileSets`.
- Updated traceback in `viterbi` to not use nested vecs (behind `dev-phmm` feature gate)
- Removed `as_u32` from `Alignment`, since the integer type for scores has been standardized

### Fixes

- Fixes a bug in `AlignmentStates::add_ciglet()` for `inc` > 1

## [0.0.20] - 2025-08-18

### Added

- Added a portable, generic striped Smith-Waterman alignment function (full matrix version) and related convenience functions.
- Added `MaybeAligned` as the return type for alignments.
- Added `new_with_w*` constructors for `SharedProfiles` and `LocalProfiles` for balanced profile creation for a segment width.
- Added `ref_len` and `query_len` fields to the `Alignment` struct, as well as `to_reverse` and `make_reverse` methods.
- Added `IntoIterator` and `FromIterator` implementations for `AlignmentStates` as well as an `iter()` method.
- Added `StatesSequence` trait for working with slices/iterators of `Ciglet` elements.
- Added `PartialEq` between `Cigar` and `AlignmentStates` without allocating on the heap.
- Added structs `CigarView` and `CigarViewMut`.
- Added macros for easily implementing traits on enums and wrapper types (`define_whichever` and `impl_traits`).
- Added `score_from_path` and the Viterbi algorithm for `GlobalPhmm`.
- Added `CorePhmm` struct storing PHMM layers (without other end modules or model information).
- Added `LocalPhmm` and structs towards supporting other alignment modes in pHMMs.

### Changed

- `Alignment` now uses `AlignmentStates` instead of `Cigar`, and alignment algorithms return `MaybeAligned` with scores as `u32` instead of `u64`.
- Modified `ScalarProfile` for alignments to hold a reference to the weight matrix.
- Methods `into_*_profile()` for `Nucleotides` now returns balanced profiles for register widths of 256.
- Improved performance of profile creation.
- All pHMM functionality is behind the feature gate `dev-phmm`.
- Moved `clear` and `restrict` methods for views into a separate `Restrict` trait.
- Extend the lifetimes returned by `as_bytes` and other methods for non-mutable views.
- Changed the `Debug` format for `AlignmentStates`, and converted `AlignmentStates::reverse()` to `make_reverse` and `to_reverse`.
- Replace ad-hoc casting methods in `AnyInt` and `Float` with generic `cast_as` and `cast_from`.
- Improved error messages for `ScoringError`.
- Errors for FASTQ, FASTA, and SAM will save `PathBuf` information for display when opening a file.
- Implement `Ord` and `PartialOrd` for `Cigar` and `QualityScores`.
- Removes an extra allocation when calling `Cigar::is_valid`.

### Removed

- Removed `new_with_i8()` and `new_with_i16()` from `SharedProfiles` and `LocalProfiles`.
- Removed `i64` profiles from `LocalProfiles` and `SharedProfiles`.

### Fixed

- `GetCode` now correctly retrieves nested IO error codes even if the outer IO error doesn't have one.
- `sw_score_from_path` no longer panics if the query or reference are not long enough

## [0.0.19] - 2025-06-11

### Added

- Added accessor methods to `SamRow` enum variants; added `is_unmapped` to `SamData`.
- Added macros for facilitating floating point comparisons under the `fuzzing` feature
- Implemented `Arbitrary` for PHMM-related types under the `fuzzing feature
- Added `SamHmmWriter` for writing pHMMs to SAM model files
- Added `expanded_cigar_iter` for iterating over the operations in a CIGAR string
- Added `peek_op`, `next_if_op`, and other related helper methods to `CigletIterator`
- Added `IsValidDNA` strategies for non-ambiguous bases
- Exposed `AlignmentStates`, a struct for incrementally building up CIGAR strings
- Add `invert` method to `Alignment`
- Added `sw_score_from_path` to obtain the local smith waterman alignment score from a CIGAR string

### Changed

- MSRV bumped to 1.88 (nightly)
- Renamed `arbitrary` feature to `fuzzing`
- Add additional error checks and improve error messages for `FastQReader` and `FastaReader`.
- `SamData` is no longer boxed in `SamRow`; `SamAligned` and `SamInsertion` are now private.
- `merge_pair_using_reference` now includes clipping in the output.

## [0.0.18] - 2025-05-09

### Added

- An iterator of `Ciglet`s can now be collected into a `Result<Cigar, CigarError>`

### Fixed

- Specified dependencies.
- Corrected `SamParser` to handle case when all parameters are on one line

## [0.0.17] - 2025-05-06

### Added

- Introduce `ByteIndexCounts`, a generalization of `NucleotideCounts`
- Added an `Alignment` struct to improve the ergonomics of _Zoe_'s alignment algorithms

### Changed

- `from_filename` on the different readers now throw an error for empty files, and `from_readable` is introduced to enable the same check for any input type
- Implement `DoubleEndedIterator` for `CigletIterator`

## [0.0.16] - 2025-04-28

### Added

- Added new `From` and `TryFrom` implementations to sequence types
- Implemented `Arbitrary` for `KmerSet` when the `arbitrary` feature is enabled
- Added functions for retrieving the codons from `Nucleotides`
- Add `make_uppercase` to `RecodeNucleotides`
- Add `to_view` for downgrading mutable views to immutable views
- Add preliminary struct for representing profile hidden Markov models and parsing them from SAM model files.
- Added functions such as `find_next_aa` for locating aa translations within nucleotide sequences

### Changed

- `ToDNA` now has a blanket impl to support the trait being used on more types
- The `Arbitrary` implementation for `Kmer` can now generate non-graphic ASCII
- Many traits have been sealed to prevent SemVer breakage. They will be unsealed as the API stabilizes.
- Reduces restrictions for implementing `KmerEncoder`.

### Fixed

## [0.0.15] - 2025-04-11

### Added

- Added new DNA validation / recoding functions, including a SIMD-accelerated `is_acgtn_uc` function
- Added constant `DEFAULT_SIMD_LANES` to standardize SIMD lane count
- Added more comprehensive documentation and tests for DNA mapping and sanitzation functions
- Added impl for `TryFrom`, method `is_valid`, and enum `CigarError` for the CIGAR module.

### Changed

- Restructured DNA nucleotide validation and recoding using strategy enums:
  - `IsValidDNA` for validation strategies
  - `RefineDNAStrat` for filtering and/or recoding strategies
  - `RecodeDNAStrat` for recoding strategies
- Improved mapping functions for handling gaps and non-standard characters
- Improved translation code to handle user specified partial codons

### Fixed

- Fixed documentation inconsistencies in alphabet and mapping modules
- Apply U to T recoding where possible for DNA mappings

## [0.0.14] - 2025-03-18

### Added

- Implemented `Arbitrary` traits for amino acids, vecs of bounded length, and kmers
- Added `clear()` to various structs and traits.
- Added `starts_with_repeating` and `ends_with_repeating` to `ByteSubstring`
- Added consuming kmer sequence iterators
- Added `position` and `rposition` functions for range searching
- Added `align_with_cigar_iter` function
- Added `find_all_kmers` and `find_all_kmers_rev` with supporting methods

### Changed

- MSRV is now nightly 1.87; we switch to using `shift_elements_{left,right}` from upstream.
- Restructured `FrequencyTable` to use struct instead of trait
- Improved quality score documentation and `QualityStats` implementations
- Improved naming consistency and testing of alphabet mapping/checking functions.

### Fixed

- Fixes functions that map IUPAC to ACGTN like `recode_iupac_to_acgtn_uc()`.

## [0.0.13] - 2025-02-18

### Added

- Updated README notices for open-sourcing Zoe on Github!
- Added optional feature "arbitrary" to support fuzzing Zoe types
- Added supported for the _signed_ striped SW scoring function (generally faster)
- Added set operations for k-mer sets and anew trait for finding k-mers in sequences
- Added support for restricted range searching in `ByteSubstring` and `FindKmers`
- Added `make_reverse_complement_simd` and `translate_to_stop` functions
- Added `SimdAnyInt` + `FromSameSignedness` traits to improve integer generics
- Added note / warning documentation styles

### Changed

- Profile set structs now use signed integers
- `SimpleWeightMatrix` and `BiasedWeightMatrix` are now `WeightMatrix`
- Improved AA translation efficiency
- Updates integer traits and implementations for `Cigar`
- Doc improvements

## [0.0.12] - 2025-01-22

### Added

- Added a k-mer module for membership and counting queries. A 3-bit k-mer encoder is initially provided.
- View types have been added to fundamental types (e.g., `Nucleotides`) and some record types (e.g., FASTQ).
- Added new documentation and this changelog.
- Adds license info, notices, and disclaimers

### Changed

- Zoe now uses Rust 2024.
- Major refactor of data modules, including visibility and traits implemented.
- The standard genetic code now uses `static` instead of lazy statics.
- Upates Cargo.toml metadata and revises how contributors are accounted

## [0.0.11] - 2024-12-20

- **Changed**: Updated documentation.
- **Fixed**: Bug fix for `ByteIndexMap` by adding synonym functionality to T and U.

## [0.0.10] - 2024-12-12

### Added

- Adds an initial Striped Smith-Waterman score algorithm.
- Lazy profile building ergonomics for unsigned integer types are also included.
- Created `ByteIndexMap` type and related ergonomics for mapping bytes.
- Adds a Bibliography document and updates README with a FAQ.

## [0.0.9] - 2024-12-11

### Changed

- Optimizes quality score median calculations via a new trait.
- The `test` feature is now conditional and not required for building Zoe.
- Updates the documentation.

## [0.0.8] - 2024-10-30

- **Changed**: Reduces allocations in `fuzzy_substring_match_simd`.
- **Fixed**: Corrects a performance-related off-by-one error in the fuzzy search function.

## [0.0.7] - 2024-10-04

- **Added**: Added optional dependency `multiversion` for SIMD functions.
- **Changed**: Improves correct handling of FASTA and FASTQ formats. The `rand` feature is now optional.

## [0.0.6] - 2024-09-18

- **Changed**: Improves amino acid translation using a custom perfect hash.
- **Removed**: Removes `ahash` dependency.

## [0.0.5] - 2024-09-17

- **Added**: Introduces the `find_k_repeating` algorithm to the search module.
- **Changed**: Improves correctness of Cigar string processing.
- **Removed**: Removes `atoi` dependency.

## [0.0.4] - 2024-09-12

- **Added**: DNA distance functions such as JC69 and TN93.

## [0.0.3] - 2024-09-04

- **Changed**: Improved string search capability and some re-organization.

## [0.0.2] - 2024-04-23

- **Changed**: A major refactor of Zoe.
- **Added**: Some search, validation, and alignment features put in place.

## [0.0.1] - 2023-12-07

- **Added**: Initial internal release. Provides various readers and types for bioinformatics data manipulation.

<!-- Versions -->
[0.0.24]: https://github.com/CDCgov/zoe/compare/v0.0.23...v0.0.24
[0.0.23]: https://github.com/CDCgov/zoe/compare/v0.0.22...v0.0.23
[0.0.22]: https://github.com/CDCgov/zoe/compare/v0.0.21...v0.0.22
[0.0.21]: https://github.com/CDCgov/zoe/compare/v0.0.20...v0.0.21
[0.0.20]: https://github.com/CDCgov/zoe/compare/v0.0.19...v0.0.20
[0.0.19]: https://github.com/CDCgov/zoe/compare/v0.0.18...v0.0.19
[0.0.18]: https://github.com/CDCgov/zoe/compare/v0.0.17...v0.0.18
[0.0.17]: https://github.com/CDCgov/zoe/compare/v0.0.16...v0.0.17
[0.0.16]: https://github.com/CDCgov/zoe/compare/v0.0.15...v0.0.16
[0.0.15]: https://github.com/CDCgov/zoe/compare/v0.0.14...v0.0.15
[0.0.14]: https://github.com/CDCgov/zoe/compare/v0.0.13...v0.0.14
[0.0.13]: https://github.com/CDCgov/zoe/compare/v0.0.12...v0.0.13
[0.0.12]: https://github.com/CDCgov/zoe/compare/v0.0.11...v0.0.12
[0.0.11]: https://github.com/CDCgov/zoe/compare/v0.0.10...v0.0.11
[0.0.10]: https://github.com/CDCgov/zoe/compare/v0.0.9...v0.0.10
[0.0.9]: https://github.com/CDCgov/zoe/compare/v0.0.8...v0.0.9
[0.0.8]: https://github.com/CDCgov/zoe/compare/v0.0.7...v0.0.8
[0.0.7]: https://github.com/CDCgov/zoe/compare/v0.0.6...v0.0.7
[0.0.6]: https://github.com/CDCgov/zoe/compare/v0.0.5...v0.0.6
[0.0.5]: https://github.com/CDCgov/zoe/compare/v0.0.4...v0.0.5
[0.0.4]: https://github.com/CDCgov/zoe/compare/v0.0.3...v0.0.4
[0.0.3]: https://github.com/CDCgov/zoe/compare/v0.0.2...v0.0.3
[0.0.2]: https://github.com/CDCgov/zoe/compare/v0.0.1...v0.0.2
[0.0.1]: https://github.com/CDCgov/zoe/releases/tag/v0.0.1

<!-- Links -->
[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html
