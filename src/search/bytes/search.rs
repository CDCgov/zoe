use crate::simd::SimdMaskFunctions;
use std::simd::{LaneCount, SupportedLaneCount, prelude::*};

/// Trait for splitting byte strings by byte(s).
pub(crate) trait ByteSplitIter<'a> {
    /// Lazily splits the haystack only on the `\r` and `\n` characters.
    fn lines_ascii<const N: usize>(&'a self) -> SplitByByte2<'a, N>
    where
        LaneCount<N>: SupportedLaneCount;
}

impl<'a, T: AsRef<[u8]>> ByteSplitIter<'a> for T {
    fn lines_ascii<const N: usize>(&'a self) -> SplitByByte2<'a, N>
    where
        LaneCount<N>: SupportedLaneCount, {
        let haystack = self.as_ref();
        SplitByByte2 {
            haystack,
            done: false,
            b1: b'\n',
            b2: b'\r',
        }
    }
}

pub(crate) struct SplitByByte2<'a, const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) haystack: &'a [u8],
    pub(crate) done:     bool,
    pub(crate) b1:       u8,
    pub(crate) b2:       u8,
}

impl<const N: usize> SplitByByte2<'_, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    pub fn remaining_len(&self) -> usize {
        self.haystack.len()
    }
}

/// Based on the implementation in for [`std::slice::Split`] but using
/// [`position_by_byte2`] internally.
impl<'a, const N: usize> Iterator for SplitByByte2<'a, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<&'a [u8]> {
        if self.done {
            return None;
        }

        if let Some(idx) = position_by_byte2(self.haystack, self.b1, self.b2) {
            let (slice_found, remaining) = (&self.haystack[..idx], &self.haystack[idx + 1..]);
            self.haystack = remaining;
            Some(slice_found)
        } else {
            self.done = true;
            Some(self.haystack)
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.done {
            (0, Some(0))
        } else {
            (1, Some(self.haystack.len() + 1))
        }
    }
}

/// Finds the index of the byte `b` in the `haystack`. The const parameter `N`
/// specifies the number of SIMD lanes used.
///
/// See [`position_simd`] for the SIMD byte search leveraged internally.
///
#[must_use]
#[inline]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn position_by_byte<const N: usize>(haystack: &[u8], b: u8) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    let (pre, mid, post) = haystack.as_simd();

    if let Some(p) = pre.iter().position(|x| *x == b) {
        Some(p)
    } else if let Some(p) = position_simd::<N, 4, _>(mid, |v| v.simd_eq(Simd::from_array([b; N]))) {
        Some(p + pre.len())
    } else {
        post.iter().position(|x| *x == b).map(|p| pre.len() + (mid.len() * N) + p)
    }
}

/// Finds the index of the bytes `b1` OR `b2` in the `haystack`. The const
/// parameter `N` specifies the number of SIMD lanes used.
///
/// See [`position_simd`] for the SIMD byte search leveraged internally.
///
#[must_use]
#[inline]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn position_by_byte2<const N: usize>(haystack: &[u8], b1: u8, b2: u8) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    let (pre, mid, post) = haystack.as_simd();

    if let Some(p) = pre.iter().position(|x| *x == b1 || *x == b2) {
        Some(p)
    } else if let Some(p) = position_simd::<N, 3, _>(mid, |v| {
        v.simd_eq(Simd::from_array([b1; N])) | v.simd_eq(Simd::from_array([b2; N]))
    }) {
        Some(p + pre.len())
    } else {
        post.iter()
            .position(|x| *x == b1 || *x == b2)
            .map(|p| pre.len() + mid.len() * N + p)
    }
}

/// Finds the index of the bytes `b1` OR `b2` or `b3` in the `haystack`. The const
/// parameter `N` specifies the number of SIMD lanes used.
///
/// See [`position_simd`] for the SIMD byte search leveraged internally.
///
#[must_use]
#[inline]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn position_by_byte3<const N: usize>(haystack: &[u8], b1: u8, b2: u8, b3: u8) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
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
/// The SIMD lanes `N`, unroll factor `UF` are const parameters that can be used
/// to adjust performance characteristics.
///
/// ## Acknowledgements
///
/// The unrolling algorithm inspired from previous work in the excellent [memchr crate](https://crates.io/crates/memchr).
///
#[allow(clippy::needless_range_loop)]
#[must_use]
#[inline]
pub fn position_simd<const N: usize, const UF: usize, P>(haystack: &[Simd<u8, N>], predicate: P) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount,
    P: Fn(Simd<u8, N>) -> Mask<i8, N>, {
    // unroll factor

    let chunk_size = N * UF;

    let chunks = haystack.chunks_exact(UF);
    let rem = chunks.remainder();
    let mut mask_buffer = [Mask::from_array([false; N]); UF];

    for (i, c) in chunks.enumerate() {
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

#[cfg(test)]
mod bench {
    use super::*;
    use std::{iter::once, sync::LazyLock};
    use test::Bencher;

    extern crate test;

    static LONG: LazyLock<Vec<u8>> = LazyLock::new(|| b"0".repeat(288).into_iter().chain(once(b'1')).collect());
    static SHORT: LazyLock<Vec<u8>> = LazyLock::new(|| b"0".repeat(11).into_iter().chain(once(b'1')).collect());

    #[bench]
    fn position_short_scalar(b: &mut Bencher) {
        b.iter(|| SHORT.iter().position(|x| *x == b'1'));
    }

    #[bench]
    fn position_long_scalar(b: &mut Bencher) {
        b.iter(|| LONG.iter().position(|x| *x == b'1'));
    }

    #[bench]
    fn position_short_simd(b: &mut Bencher) {
        b.iter(|| position_by_byte::<16>(&SHORT, b'1'));
    }

    #[bench]
    fn position_long_simd(b: &mut Bencher) {
        b.iter(|| position_by_byte::<16>(&LONG, b'1'));
    }
}
