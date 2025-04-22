pub(crate) mod aa;
pub(crate) mod byte_index;
pub(crate) mod dna;
pub(crate) mod gc;

#[cfg(test)]
mod test;

pub use aa::*;
pub use byte_index::*;
pub use dna::*;
pub use gc::*;

/// A boolean mapping of possible CIGAR string characters (not counting `*`)
pub(crate) const IS_CIGAR: [bool; 256] = make_is_alpha_mapping(b"MDNX=ISHP0123456789");

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

/// Utility function for making a mapping with a default value.
const fn make_mapping_with_default<const N: usize>(from_byte: &[u8; N], dest_byte: &[u8; N], all_others: u8) -> [u8; 256] {
    let mut mapping = [all_others; 256];
    let mut i = 0;

    while i < N {
        mapping[from_byte[i] as usize] = dest_byte[i];
        i += 1;
    }
    mapping
}

/// Utility function for making a mapping but keeping the original value otherwise.
#[allow(clippy::cast_possible_truncation)]
const fn make_mapping_otherwise_self<const N: usize>(from_byte: &[u8; N], dest_byte: &[u8; N]) -> [u8; 256] {
    let mut mapping = [0u8; 256];

    let mut i = 0;
    while i < mapping.len() {
        mapping[i] = i as u8;
        i += 1;
    }

    i = 0;
    while i < N {
        mapping[from_byte[i] as usize] = dest_byte[i];
        i += 1;
    }
    mapping
}
