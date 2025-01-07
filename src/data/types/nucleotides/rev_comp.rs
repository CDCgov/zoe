use crate::{
    data::mappings::TO_REVERSE_COMPLEMENT,
    simd::{SimdByteFunctions, SimdMaskFunctions},
};
use std::simd::{LaneCount, SupportedLaneCount};

/// Performs the DNA reverse complement of the byte slice in place.
/// Assumes ASCII input.
#[inline]
pub fn make_reverse_complement(bases: &mut [u8]) {
    bases.reverse();
    for x in bases {
        *x = TO_REVERSE_COMPLEMENT[*x as usize];
    }
}

/// Performs the DNA reverse complement of the byte slice into a new vector.
/// Assumes ASCII input.
#[inline]
#[must_use]
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .copied()
        .map(|x| TO_REVERSE_COMPLEMENT[x as usize])
        .collect()
}

/// Reverse complement of a nucleotide sequences using explicit SIMD
/// instructions.
///
/// # Note
///
/// Works well for Haswell and above architectures. Recommend 32 lanes for x86.
/// Both [`swizzle_dyn`](std::simd::prelude::Simd::swizzle_dyn) and
/// [`gather_or`](std::simd::prelude::Simd::gather_or) are too slow to be
/// relevant versus the scalar implementation.
#[inline]
#[must_use]
#[allow(dead_code)]
pub fn reverse_complement_simd<const N: usize>(bases: &[u8]) -> Vec<u8>
where
    LaneCount<N>: SupportedLaneCount, {
    let (pre, mid, sfx) = bases.as_simd::<N>();
    let mut reverse_complement = Vec::with_capacity(bases.len());

    reverse_complement.extend(sfx.iter().rev().copied().map(|x| TO_REVERSE_COMPLEMENT[x as usize]));

    reverse_complement.extend(
        mid.iter()
            .map(|&v| {
                let mut rev = v.reverse();
                let lowercase = rev.is_ascii_lowercase();
                rev = lowercase.make_selected_ascii_uppercase(&rev);

                rev.exchange_byte_pairs(b'T', b'A');
                rev.exchange_byte_pairs(b'G', b'C');
                rev.exchange_byte_pairs(b'R', b'Y');
                rev.exchange_byte_pairs(b'K', b'M');
                rev.exchange_byte_pairs(b'B', b'V');
                rev.exchange_byte_pairs(b'H', b'D');
                rev.if_value_then_replace(b'U', b'A');

                rev = lowercase.make_selected_ascii_lowercase(&rev);
                rev.to_array()
            })
            .rev()
            .flatten(),
    );

    reverse_complement.extend(pre.iter().rev().copied().map(|x| TO_REVERSE_COMPLEMENT[x as usize]));

    reverse_complement
}
