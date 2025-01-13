use std::{
    cmp::min,
    simd::{LaneCount, SupportedLaneCount, prelude::*},
};

use crate::math::Uint;

/// Calculates the number of differences at the byte (or base/residue) level. An
/// unsigned integer may be supplied as a generic argument. Smaller sizes like
/// `u32` have been known to auto-vectorize more nicely.
///
///
/// # Example
/// ```
/// # use zoe::distance::hamming;
///
/// let s1 = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
/// let s2 = b"ATGCATnGATCGATCGATCGAnCGATCGATnC";
///
/// assert!(3 == hamming::<u32>(s1, s2));
/// ```
#[must_use]
pub fn hamming<T: Uint>(a: &[u8], b: &[u8]) -> usize {
    std::iter::zip(a, b)
        .map(|(a, b)| a != b)
        .fold(T::ZERO, |a, b| a + T::from_bool(b))
        .as_usize()
}

/// Calculates the number of differences at the byte (or base/residue) level.
///
/// The const argument `N` is the number of SIMD lanes to use at a time. Unless
/// compiling for a specific platform, a choice of 16 or 32 is a good value.
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
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn hamming_simd<const N: usize>(x: &[u8], y: &[u8]) -> usize
where
    LaneCount<N>: SupportedLaneCount, {
    let mut matches: usize = 0;
    let limit = min(x.len(), y.len());

    let (x, y) = if x.len() < y.len() {
        (x, &y[..x.len()])
    } else {
        (&x[..y.len()], y)
    };

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

#[cfg(test)]
mod test {
    use crate::distance::hamming;

    use super::hamming_simd;

    #[test]
    fn test_hamming_simd() {
        use crate::distance::dna::test::{LONG_READ_1 as x, LONG_READ_2 as y};

        for simd_distance in [hamming_simd::<16>(x, y), hamming_simd::<32>(x, y), hamming_simd::<4>(x, y)] {
            assert_eq!(simd_distance, hamming::<usize>(x, y));
            assert_eq!(simd_distance, hamming::<u32>(x, y));
        }

        assert_eq!(hamming_simd::<32>(&x[..20], y), hamming::<usize>(&x[..20], y));
    }
}
