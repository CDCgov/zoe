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
