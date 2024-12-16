# Zoe

[Zoe](https://en.wikipedia.org/wiki/Zoe_(name)) provides both broad and highly specialized implementations for
bioinformatics. In particular, we focus on common data formats and methods
relevant for the sequencing of RNA viruses.

**This library runs on Rust *nightly*.**

## Distinctives

*Zoe* has design goals that may or may not be compatible with your project:

- Embracing Rust **Nightly** to enable early access to features and fixes, especially [portable SIMD](https://github.com/rust-lang/portable-simd)
- Providing broad support for the authors' activities over creating a multiplicity of smaller projects.
- Avoiding `unsafe` as well as non-portable implementations
- Avoiding large transitive dependency trees and binary sizes
- Providing reasonably performant solutions while not sacrificing the above
- Porting diverse, relevant methods and providing clear attribution to the original authors

## FAQ

### Yes, yes, but will *Zoe* ever move to stable Rust?

*Zoe* has used a number of nightly features over its lifetime. Some have been stabilized (e.g., `let_chains`,`lazy_cell`,`const_fn_floating_point_arithmetic`) while others have been removed (`iter_collect_into`), avoided (`array_chunks`, `iter_array_chunks`), or forbidden due to incompleteness (`generic_const_exprs`).

In order for *Zoe* to move to stable, the portable SIMD [nightly feature](https://github.com/rust-lang/rust/issues/86656) would first need to be stabilized since it is an essential component of the library. After that, *Zoe* could march towards stable Rust by gating, working around or removing any remaining nightly features. Notwithstanding, development of *Zoe* is likely to continue on nightly, if only for convenient access to nightly bencher, the latest LLVM code generation, and convenient use of rustfmt's unstable [features](https://rust-lang.github.io/rustfmt/?version=master&search=).

### What functionality uses Portable SIMD?

This list may not be up-to-date or exhaustive, but some example usage:

- Alignment functions: `sw_score_simd` (TODO)
- Composition functions: [`composition::gc_content_simd`], [`composition::AlignmentComposition`]
- Distance functions: [`distance::hamming_simd`], [`distance::p_distance_acgt`]
- Search
  - Substring: [`search::fuzzy_substring_match_simd`], [`search::substring_match_simd`], [`search::find_k_repeating`]
  - Byte: [`search::position_by_byte`] (and other variants) and [`search::replace_all_bytes_simd`]
- Transformation functions: [`data::types::nucleotides::reverse_complement_simd`]
- Validation functions: [`data::CheckSequence`]

### How should this library be cited?

While it is always welcome to cite *Zoe* as the implementation (see `CITATION.bib`), we encourage authors to also cite the **original works** for the method being used. To that end, the innermost function (not convenience methods) usually contains all relevant citations to it (raise an issue if not). In addition, the complete bibliography of methods used in this library is included in [BibTex](https://en.wikipedia.org/wiki/BibTeX) format for easy inclusion in a mansuscript (see `BIBLIOGRAPHY.bib`).

## Crate Features

Zoe is designed to minimize large dependency trees and binary sizes while providing
performant solutions. Several dependencies are optional and can be enabled as features
at the user's discretion.

### Default Features

- **rand**: Enabled by default. This optional feature provides support for random sequence generation for nucleotides and amino acids. Adds the dependency [rand_xoshiro](https://docs.rs/rand_xoshiro/latest/rand_xoshiro/).

### Optional Features

- **multiversion**: This optional feature adds multiversioning of select Zoe SIMD functions using the [multiversion](https://docs.rs/multiversion/latest/multiversion/) crate.
