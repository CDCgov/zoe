use std::{
    cmp::min,
    simd::{prelude::*, LaneCount, SupportedLaneCount},
};

use crate::simd::SimdByteFunctions;

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

#[allow(clippy::cast_precision_loss)]
#[must_use]
pub fn p_distance_acgt<const N: usize>(x: &[u8], y: &[u8]) -> Option<f64>
where
    LaneCount<N>: SupportedLaneCount, {
    let alpha_v = [Simd::splat(b'A'), Simd::splat(b'C'), Simd::splat(b'G'), Simd::splat(b'T')];

    let (x, y) = if x.len() < y.len() {
        (x, &y[..x.len()])
    } else {
        (&x[..y.len()], y)
    };

    let mut mismatches: usize = 0;
    let mut valid_length = 0;

    let mut x = x.chunks_exact(N * 255);
    let mut y = y.chunks_exact(N * 255);

    for (c1, c2) in x.by_ref().zip(y.by_ref()) {
        let mut accum: Simd<u8, N> = Simd::splat(0);

        let mut c1 = c1.chunks_exact(N);
        let mut c2 = c2.chunks_exact(N);

        for (v1, v2) in c1.by_ref().zip(c2.by_ref()) {
            let mut v1: Simd<u8, N> = Simd::from_slice(v1).to_ascii_uppercase();
            let mut v2: Simd<u8, N> = Simd::from_slice(v2).to_ascii_uppercase();
            v1.exchange_byte_pairs(b'U', b'T');
            v2.exchange_byte_pairs(b'U', b'T');

            let mut valid = Mask::from_array([false; N]);
            for a in alpha_v {
                valid |= a.simd_eq(v1);
                valid |= a.simd_eq(v2);
            }

            valid_length += valid.to_bitmask().count_ones() as usize;

            let m = (v1.simd_ne(v2) & valid).to_int();
            // True => -1, so - -1 => +1
            accum -= m.cast();
        }

        let accum2: Simd<u16, N> = accum.cast();
        mismatches += accum2.reduce_sum() as usize;
    }

    let x = x.remainder();
    let y = y.remainder();
    let mut accum: Simd<u8, N> = Simd::splat(0);
    let mut c1 = x.chunks_exact(N);
    let mut c2 = y.chunks_exact(N);

    for (v1, v2) in c1.by_ref().zip(c2.by_ref()) {
        let mut v1: Simd<u8, N> = Simd::from_slice(v1).to_ascii_uppercase();
        let mut v2: Simd<u8, N> = Simd::from_slice(v2).to_ascii_uppercase();
        v1.exchange_byte_pairs(b'U', b'T');
        v2.exchange_byte_pairs(b'U', b'T');
        let mut valid1 = Mask::from_array([false; N]);
        let mut valid2 = Mask::from_array([false; N]);
        for a in alpha_v {
            valid1 |= a.simd_eq(v1);
            valid2 |= a.simd_eq(v2);
        }
        let valid = valid1 & valid2;

        valid_length += valid.to_bitmask().count_ones() as usize;

        let m = (v1.simd_ne(v2) & valid).to_int();
        // True => -1, so - -1 => +1
        accum -= m.cast();
    }
    let accum2: Simd<u16, N> = accum.cast();
    mismatches += accum2.reduce_sum() as usize;

    let r1 = c1.remainder();
    let r2 = c2.remainder();

    mismatches += r1
        .iter()
        .zip(r2)
        .filter(|(a, b)| {
            let mut a = a.to_ascii_uppercase();
            let mut b = b.to_ascii_uppercase();
            if a == b'U' {
                a = b'T';
            }

            if b == b'U' {
                b = b'T';
            }
            let valid = matches!(a, b'A' | b'G' | b'T' | b'C') && matches!(b, b'A' | b'G' | b'T' | b'C');
            valid_length += usize::from(valid);
            a != b && valid
        })
        .count();

    if valid_length > 0 {
        Some(mismatches as f64 / valid_length as f64)
    } else {
        None
    }
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
    use crate::distance::{hamming, hamming_u32};

    use super::hamming_simd;

    #[test]
    fn test_hamming_simd() {
        use crate::distance::dna::test::{LONG_READ_1 as x, LONG_READ_2 as y};

        for simd_distance in [hamming_simd::<16>(x, y), hamming_simd::<32>(x, y), hamming_simd::<4>(x, y)] {
            assert_eq!(simd_distance, hamming(x, y));
            assert_eq!(simd_distance, hamming_u32(x, y) as usize);
        }

        assert_eq!(hamming_simd::<32>(&x[..20], y), hamming(&x[..20], y));
    }
}
