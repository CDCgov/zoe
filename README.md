# The Zoe Anthology

[Zoe](https://en.wikipedia.org/wiki/Zoe_(name)) provides both broad and highly specialized implementations for
bioinformatics. In particular, we focus on common data formats and methods relevant for the sequencing of RNA viruses. Read the [documentation](https://cdcgov.github.io/zoe).

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

### How do I add this crate?

Provided *Zoe* has published the latest version to [crates.io](https://crates.io/), within your project you can simply use: `cargo add zoe`

If you would like to pin to a Github tag or commit, you could use the following within your `Cargo.toml` project file:

```toml
[dependencies]
# Fill in VERSION_TAG_HERE with tagged version:
zoe = { git = "https://github.com/CDCgov/zoe", tag = "VERSION_TAG_HERE"}

# OR, fill in MAIN_COMMIT_HERE:
zoe = { git = "https://github.com/CDCgov/zoe", commit = "MAIN_COMMIT_HERE"}
```

*CAUTION: because feature branches get routinely deleted, pinning to feature branches (via branch, commit or tag) should only be used for active development or forks!*

### Yes, yes, but will *Zoe* ever move to stable Rust?

*Zoe* has used a number of nightly features over its lifetime. Some have been stabilized (e.g., <!--`let_chains`, -->`lazy_cell`,`const_fn_floating_point_arithmetic`) while others have been removed, avoided, or forbidden due to incompleteness.

In order for *Zoe* to move to stable, the portable SIMD [nightly feature](https://github.com/rust-lang/rust/issues/86656) would first need to be stabilized since it is an essential component of the library. After that, *Zoe* could march towards stable Rust by gating, working around or removing any remaining nightly features. Notwithstanding, development of *Zoe* is likely to continue on nightly, if only for convenient access to nightly bencher, the latest LLVM code generation, and convenient use of rustfmt's unstable [features](https://rust-lang.github.io/rustfmt/?version=master&search=).

### What functionality uses Portable SIMD?

This list may not be up-to-date or exhaustive, but some example usage:

- Alignment functions: [`alignment::sw::sw_simd_score`]
- Composition functions: [`composition::gc_content_simd`], [`composition::AlignmentComposition`]
- Distance functions: [`distance::hamming_simd`], [`distance::dna::p_distance_acgt`]
- Search
  - Substring: [`search::fuzzy_substring_match_simd`], [`search::substring_match_simd`], [`search::find_k_repeating`]
  - Byte: [`search::position_by_byte`] (and other variants) and [`search::replace_all_bytes_simd`]
- Transformation functions: [`data::types::nucleotides::reverse_complement_simd`]
- Validation functions: [`data::CheckSequence`]

### How should this library be cited?

While it is always welcome to cite *Zoe* as the implementation (see [`CITATION.bib`](https://github.com/CDCgov/zoe/blob/main/BIBLIOGRAPHY.bib))for the method being used. To that end, the innermost function (not convenience methods) usually contains all relevant citations to it (raise an issue if not). In addition, the complete bibliography of methods used in this library is included in [BibTex](https://en.wikipedia.org/wiki/BibTeX) format for easy inclusion in a mansuscript (see [`BIBLIOGRAPHY.bib`](https://github.com/CDCgov/zoe/blob/main/BIBLIOGRAPHY.bib)).

## Crate Features

Zoe is designed to minimize large dependency trees and binary sizes while providing
performant solutions. Several dependencies are optional and can be enabled as features
at the user's discretion.

### Default Features

- **rand**: Enabled by default. This optional feature provides support for random sequence generation for nucleotides and amino acids. Adds the dependency [rand_xoshiro](https://docs.rs/rand_xoshiro/latest/rand_xoshiro/).

### Optional Features

- **multiversion**: This optional feature adds multiversioning of select Zoe SIMD functions using the [multiversion](https://docs.rs/multiversion/latest/multiversion/) crate.
- **fuzzing**: This optional feature adds several testing features to Zoe, including [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html) implementations for many of Zoe's data types to support fuzz testing. It also introduces many wrapper types with custom [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html) implementations providing data with stronger assumptions.

## Notices

### Contact Info

For direct correspondence on the project, feel free to contact: [Samuel S. Shepard](mailto:sshepard@cdc.gov), Centers for Disease Control and Prevention or reach out to other [contributors](CONTRIBUTORS.md).

### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).  All contributions to this repository will be released under the CC0 dedication.  By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

### License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later. This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version. This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details. You should have received a copy of the Apache Software License along with this program. If not, see: <http://www.apache.org/licenses/LICENSE-2.0.html>. The source code forked from other open source projects will inherit its license.

### Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the [Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md). For more information about CDC's privacy policy, please visit <http://www.cdc.gov/other/privacy.html>.

### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git).) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
