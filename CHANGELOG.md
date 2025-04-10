# Zoe Changelog

All notable changes to this project will be documented in this file. The format is roughly based on [Keep a Changelog], and this project tries to adheres to [Semantic Versioning].

## [0.0.16] - TBD

### Added

- Added new `From` and `TryFrom` implementations to sequence types
- Implemented `Arbitrary` for `KmerSet` when the `arbitrary` feature is enabled
- Added functions for retrieving the codons from `Nucleotides`
- Add `make_uppercase` to `RecodeNucleotides`

### Changed

- `ToDNA` now has a blanket impl to support the trait being used on more types
- The `Arbitrary` implementation for `Kmer` can now generate non-graphic ASCII

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
- Added supported for the *signed* striped SW scoring function (generally faster)
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
