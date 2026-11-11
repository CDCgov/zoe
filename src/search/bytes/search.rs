use crate::simd::SimdMaskFunctions;
use std::simd::prelude::*;

/// Finds the index of the byte `b` in the `haystack`.
///
/// See [`position_simd`] for the SIMD byte search leveraged internally.
///
/// ## Parameters
///
/// `N`: The number of SIMD lanes to use.
#[inline]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn position_by_byte<const N: usize>(haystack: &[u8], b: u8) -> Option<usize> {
    position_by_byte_inner::<N>(haystack, b)
}

/// Finds the index of the byte `b` in the `haystack`. The `haystack` is lazily
/// mapped using `simd_transform` and `byte_transform`.
///
/// See [`position_simd`] for the SIMD byte search leveraged internally.
///
/// ## Parameters
///
/// - `N`: The number of SIMD lanes to use.
/// - `S`: The SIMD-vectorized mapping to lazily apply to the `haystack`.
/// - `B`: The non-vectorized mapping to lazily apply to the `haystack`.
#[inline]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn position_by_byte_mapped<const N: usize, S, B>(
    haystack: &[u8], b: u8, simd_transform: S, byte_transform: B,
) -> Option<usize>
where
    S: Fn(Simd<u8, N>) -> Simd<u8, N>,
    B: Fn(u8) -> u8, {
    position_by_byte_mapped_inner(haystack, b, simd_transform, byte_transform)
}

/// Similar to [`position_by_byte`], but without multiversioning (for use as a
/// helper function inside other multiversioned functions).
#[inline]
#[must_use]
pub(crate) fn position_by_byte_inner<const N: usize>(haystack: &[u8], b: u8) -> Option<usize> {
    position_by_byte_mapped_inner(haystack, b, |v: Simd<u8, N>| v, |b| b)
}

/// Similar to [`position_by_byte_mapped`], but without multiversioning (for use
/// as a helper function inside other multiversioned functions).
#[inline]
#[must_use]
pub(crate) fn position_by_byte_mapped_inner<const N: usize, S, B>(
    haystack: &[u8], b: u8, simd_transform: S, byte_transform: B,
) -> Option<usize>
where
    S: Fn(Simd<u8, N>) -> Simd<u8, N>,
    B: Fn(u8) -> u8, {
    let (pre, mid, post) = haystack.as_simd();

    if let Some(p) = pre.iter().copied().map(&byte_transform).position(|x| x == b) {
        Some(p)
    } else if let Some(p) = position_simd::<N, 4, _>(mid, |v| simd_transform(v).simd_eq(Simd::from_array([b; N]))) {
        Some(p + pre.len())
    } else {
        post.iter()
            .copied()
            .map(byte_transform)
            .position(|x| x == b)
            .map(|p| pre.len() + (mid.len() * N) + p)
    }
}

/// Finds the index of the bytes `b1` or `b2` in the `haystack`.
///
/// See [`position_simd`] for the SIMD byte search leveraged internally.
///
/// ## Parameters
///
/// `N`: The number of SIMD lanes to use.
#[inline]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn position_by_byte2<const N: usize>(haystack: &[u8], b1: u8, b2: u8) -> Option<usize> {
    position_by_byte2_inner::<N>(haystack, b1, b2)
}

/// Finds the index of the bytes `b1` or `b2` in the `haystack`. The `haystack`
/// is lazily mapped using `simd_transform` and `byte_transform`.
///
/// See [`position_simd`] for the SIMD byte search leveraged internally.
///
/// ## Parameters
///
/// - `N`: The number of SIMD lanes to use.
/// - `S`: The SIMD-vectorized mapping to lazily apply to the `haystack`.
/// - `B`: The non-vectorized mapping to lazily apply to the `haystack`.
#[inline]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn position_by_byte2_mapped<const N: usize, S, B>(
    haystack: &[u8], b1: u8, b2: u8, simd_transform: S, byte_transform: B,
) -> Option<usize>
where
    S: Fn(Simd<u8, N>) -> Simd<u8, N>,
    B: Fn(u8) -> u8, {
    position_by_byte2_mapped_inner(haystack, b1, b2, simd_transform, byte_transform)
}

/// Similar to [`position_by_byte2`], but without multiversioning (for use as a
/// helper function inside other multiversioned functions).
#[inline]
#[must_use]
pub(crate) fn position_by_byte2_inner<const N: usize>(haystack: &[u8], b1: u8, b2: u8) -> Option<usize> {
    position_by_byte2_mapped_inner(haystack, b1, b2, |v: Simd<u8, N>| v, |b| b)
}

/// Similar to [`position_by_byte2`], but without multiversioning (for use as a
/// helper function inside other multiversioned functions).
#[inline]
#[must_use]
pub(crate) fn position_by_byte2_mapped_inner<const N: usize, S, B>(
    haystack: &[u8], b1: u8, b2: u8, simd_transform: S, byte_transform: B,
) -> Option<usize>
where
    S: Fn(Simd<u8, N>) -> Simd<u8, N>,
    B: Fn(u8) -> u8, {
    let (pre, mid, post) = haystack.as_simd();

    if let Some(p) = pre.iter().copied().map(&byte_transform).position(|x| x == b1 || x == b2) {
        Some(p)
    } else if let Some(p) = position_simd::<N, 3, _>(mid, |v| {
        let mapped = simd_transform(v);
        mapped.simd_eq(Simd::from_array([b1; N])) | mapped.simd_eq(Simd::from_array([b2; N]))
    }) {
        Some(p + pre.len())
    } else {
        post.iter()
            .copied()
            .map(byte_transform)
            .position(|x| x == b1 || x == b2)
            .map(|p| pre.len() + mid.len() * N + p)
    }
}

/// Finds the index of the bytes `b1` or `b2` or `b3` in the `haystack`.
///
/// See [`position_simd`] for the SIMD byte search leveraged internally.
///
/// ## Parameters
///
/// `N`: The number of SIMD lanes to use.
#[inline]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn position_by_byte3<const N: usize>(haystack: &[u8], b1: u8, b2: u8, b3: u8) -> Option<usize> {
    let (pre, mid, post) = haystack.as_simd();

    if let Some(p) = pre.iter().position(|x| *x == b1 || *x == b2 || *x == b3) {
        Some(p)
    } else if let Some(p) = position_simd::<N, 1, _>(mid, |v| {
        v.simd_eq(Simd::from_array([b1; N])) | v.simd_eq(Simd::from_array([b2; N])) | v.simd_eq(Simd::from_array([b3; N]))
    }) {
        Some(p + pre.len())
    } else {
        post.iter()
            .position(|x| *x == b1 || *x == b2 || *x == b3)
            .map(|p| pre.len() + mid.len() * N + p)
    }
}

/// Searches the `haystack` using the provide SIMD byte `predicate` returning
/// the index found, or [`None`] otherwise.
///
/// ## Acknowledgements
///
/// The unrolling algorithm inspired from previous work in the excellent [memchr
/// crate](https://crates.io/crates/memchr).
///
/// ## Parameters
///
/// The following parameters can be used to adjust performance characteristics:
///
/// - `N`: The number of SIMD lanes to use.
/// - `UF`: The unroll factor, which must be non-zero.
/// - `P`: The SIMD-vectorized predicate for identifying matching bytes.
#[inline]
#[must_use]
#[allow(clippy::needless_range_loop)]
pub fn position_simd<const N: usize, const UF: usize, P>(haystack: &[Simd<u8, N>], predicate: P) -> Option<usize>
where
    P: Fn(Simd<u8, N>) -> Mask<i8, N>, {
    const { assert!(UF > 0, "`UF` must be non-zero") }

    // unroll factor

    let chunk_size = N * UF;

    let (chunks, rem) = haystack.as_chunks::<UF>();
    let mut mask_buffer = [Mask::from_array([false; N]); UF];

    for (i, c) in chunks.iter().enumerate() {
        mask_buffer[0] = predicate(c[0]);
        let mut mask = mask_buffer[0];
        for j in 1..UF {
            mask_buffer[j] = predicate(c[j]);
            mask |= mask_buffer[j];
        }

        if mask.any() {
            let offset = i * chunk_size;

            for j in 0..(UF - 1) {
                if mask_buffer[j].any() {
                    return Some(offset + j * N + mask_buffer[j].bitmask_offset());
                }
            }
            return Some(offset + (UF - 1) * N + mask_buffer[UF - 1].bitmask_offset());
        }
    }

    for (i, &v) in rem.iter().enumerate() {
        let mask = predicate(v);
        if mask.any() {
            return Some(N * (haystack.len() - rem.len()) + (i * N) + mask.bitmask_offset());
        }
    }

    None
}

#[cfg(test)]
mod test {
    use crate::search::bytes::position_by_byte;

    #[test]
    fn byte_position_test() {
        let haystack = b"0000000000000001".to_vec();
        assert_eq!(
            position_by_byte::<4>(&haystack, b'1'),
            haystack.iter().position(|x| *x == b'1'),
            "result v. expected using {h}",
            h = String::from_utf8_lossy(&haystack)
        );

        for i in (0..169).step_by(3) {
            let mut haystack = vec![b'0'; i];
            haystack.push(b'1');
            assert_eq!(
                position_by_byte::<4>(&haystack, b'1'),
                haystack.iter().position(|x| *x == b'1'),
                "result v. expected using {h}",
                h = String::from_utf8_lossy(&haystack)
            );
        }

        let mut haystack = vec![b'0'; 143];
        for i in 0..143 {
            haystack[i] = b'1';
            assert_eq!(
                position_by_byte::<4>(&haystack, b'1'),
                haystack.iter().position(|x| *x == b'1'),
                "result v. expected using {h}",
                h = String::from_utf8_lossy(&haystack)
            );
            haystack[i] = b'0';
        }
    }
}
