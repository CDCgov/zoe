use crate::data::err::{DistanceError, DistanceError::*};
use std::simd::prelude::*;

#[must_use]
/// Calculates the number of differences at the byte (or base/residue) level.
///
/// # Example
/// ```
/// use zoe::distance::hamming;
///
/// let s1 = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
/// let s2 = b"ATGCATnGATCGATCGATCGAnCGATCGATnC";
///
/// assert!(3 == hamming(s1, s2));
/// ```
///
pub fn hamming(x: &[u8], y: &[u8]) -> usize {
    x.iter().zip(y).filter(|(a, b)| a != b).count()
}

/// Calculates the number of differences at the byte (or base/residue) level up to `u32`
/// differenes. This version is faster than the `usize` version.
///
/// # Example
/// ```
/// use zoe::distance::hamming_u32;
///
/// let s1 = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
/// let s2 = b"ATGCATnGATCGATCGATCGAnCGATCGATnC";
///
/// assert!(3 == hamming_u32(s1, s2));
/// ```
// Originally proposed by ScottMCM on Zulip
#[must_use]
pub fn hamming_u32(a: &[u8], b: &[u8]) -> u32 {
    std::iter::zip(a, b).map(|(a, b)| a != b).fold(0, |a, b| a + u32::from(b))
}

use std::cmp::min;
#[allow(clippy::wildcard_imports)]
use std::simd::*;

/// Calculates the number of differences at the byte (or base/residue) level. The
/// argument `N` is the number of SIMD lanes to use at a time. Unless compiling
/// for a specific platform, a choice of 16 is a good value.
///
/// # Example
/// ```
/// use zoe::distance::hamming_simd;
///
/// let s1 = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
/// let s2 = b"ATGCATnGATCGATCGATCGAnCGATCGATnC";
///
/// assert!(3 == hamming_simd::<16>(s1, s2));
/// ```
///
#[must_use]
pub fn hamming_simd<const N: usize>(x: &[u8], y: &[u8]) -> usize
where
    LaneCount<N>: SupportedLaneCount, {
    let mut matches: usize = 0;
    let limit = min(x.len(), y.len());

    let mut x = x.chunks_exact(N * 255);
    let mut y = y.chunks_exact(N * 255);

    for (c1, c2) in x.by_ref().zip(y.by_ref()) {
        let mut accum: Simd<u8, N> = Simd::splat(0);

        let mut c1 = c1.chunks_exact(N);
        let mut c2 = c2.chunks_exact(N);

        for (v1, v2) in c1.by_ref().zip(c2.by_ref()) {
            let v1: Simd<u8, N> = Simd::from_slice(v1);
            let v2: Simd<u8, N> = Simd::from_slice(v2);
            let m = v1.simd_eq(v2).to_int();
            // True => -1, so - -1 => +1
            accum -= m.cast();
        }

        let accum2: Simd<u16, N> = accum.cast();
        matches += accum2.reduce_sum() as usize;
    }

    let x = x.remainder();
    let y = y.remainder();
    let mut accum: Simd<u8, N> = Simd::splat(0);
    let mut c1 = x.chunks_exact(N);
    let mut c2 = y.chunks_exact(N);

    for (v1, v2) in c1.by_ref().zip(c2.by_ref()) {
        let v1: Simd<u8, N> = Simd::from_slice(v1);
        let v2: Simd<u8, N> = Simd::from_slice(v2);
        let m = v1.simd_eq(v2).to_int();
        // True => -1, so - -1 => +1
        accum -= m.cast();
    }
    let accum2: Simd<u16, N> = accum.cast();
    matches += accum2.reduce_sum() as usize;

    let r1 = c1.remainder();
    let r2 = c2.remainder();
    matches += r1.iter().zip(r2).filter(|(a, b)| a == b).count();

    limit - matches
}

/// Calculates a "physiochemical" distance measure using the euclidean distances
/// over physiochemical factors. Only valid amino acids permitted in the
/// denominator.
///
/// # Example
/// ```
/// use zoe::distance::physiochemical;
///
/// let s1 = b"MANATEE";
/// let s2 = b"MANGAEE";
///
/// assert!( (physiochemical(s1,s2).unwrap()  - 1.1464952).abs() < 0.01 )
/// ```
///
/// # Citation
/// For factor analysis used by the function:
///
/// > Atchley et al. 2008. "Solving the protein sequence metric problem." Proc
/// > Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Published 2005 Apr 25.
///
///
/// # Errors
/// If either argument is empty or the sequence characters are invalid, an error
/// is thrown. See [`DistanceError`].
//
// Could make a generic function that depends on distance matrix.
#[allow(clippy::cast_precision_loss)]
pub fn physiochemical(seq1: &[u8], seq2: &[u8]) -> Result<f32, DistanceError> {
    use crate::data::matrices::PHYSIOCHEMICAL_FACTORS as dm;

    if seq1.is_empty() || seq2.is_empty() {
        return Err(NoData);
    }

    if seq1 == seq2 {
        return Ok(0.0);
    }

    let (number_valid, mut pcd_distance) = seq1
        .iter()
        .zip(seq2)
        .filter_map(|(a1, a2)| dm[*a1 as usize][*a2 as usize])
        .map(|val| (1, val))
        .fold((0u32, 0.0), |(a_valid, a_value), (i_valid, i_value)| {
            (a_valid + i_valid, a_value + i_value)
        });

    if number_valid > 0 {
        pcd_distance /= number_valid as f32;
        Ok(pcd_distance)
    } else {
        Err(NotComparable)
    }
}
