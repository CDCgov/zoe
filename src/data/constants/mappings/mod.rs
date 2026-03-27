//! Mappings for transforming sequence-like data.
//!
//! This includes byte-to-byte mappings with [`ByteMap`] and [`ByteIndexMap`],
//! as well as amino acid translation from nucleotide sequences in
//! [`StdGeneticCode`].

pub(crate) mod aa;
pub(crate) mod byte_index;
mod dna;
pub(crate) mod gc;

#[cfg(test)]
mod test;

pub use aa::*;
pub use byte_index::*;
pub use dna::*;
pub use gc::*;

/// Utility function for building *is alpha*-like maps
const fn make_is_alpha_mapping<const N: usize>(alpha: &[u8; N]) -> [bool; 256] {
    let mut mapping = [false; 256];
    let mut i = 0;

    while i < N {
        mapping[alpha[i] as usize] = true;
        i += 1;
    }
    mapping
}
