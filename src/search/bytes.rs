#![allow(dead_code)]

use std::simd::{prelude::*, LaneCount, SupportedLaneCount};

pub(crate) struct SplitByByte2<'a, const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) haystack: &'a [u8],
    pub(crate) done:     bool,
    pub(crate) b1:       u8,
    pub(crate) b2:       u8,
}

impl<'a, const N: usize> SplitByByte2<'a, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    pub fn remaining_len(&self) -> usize {
        self.haystack.len()
    }
}

// Based on the implementation in for std::slice::Split
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

pub(crate) trait ByteSplitIter<'a> {
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

#[inline]
pub(crate) fn position_by_byte2<const N: usize>(haystack: &[u8], b1: u8, b2: u8) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    let (pre, mid, post) = haystack.as_simd();

    if let Some(p) = pre.iter().position(|x| *x == b1 || *x == b2) {
        Some(p)
    } else if let Some(p) = position_simd2(mid, b1, b2) {
        Some(p + pre.len())
    } else {
        post.iter()
            .position(|x| *x == b1 || *x == b2)
            .map(|p| pre.len() + mid.len() * N + p)
    }
}

#[inline]
pub(crate) fn position_by_byte<const N: usize>(haystack: &[u8], b: u8) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    let (pre, mid, post) = haystack.as_simd();

    if let Some(p) = pre.iter().position(|x| *x == b) {
        Some(p)
    } else if let Some(p) = position_simd(mid, b) {
        Some(p + pre.len())
    } else {
        post.iter().position(|x| *x == b).map(|p| pre.len() + mid.len() * N + p)
    }
}

#[inline]
fn position_simd<const N: usize>(haystack: &[Simd<u8, N>], b: u8) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    let needle = Simd::from_array([b; N]);
    for (i, &v) in haystack.iter().enumerate() {
        let mask = needle.simd_eq(v);
        if mask.any() {
            #[cfg(target_endian = "little")]
            return Some(i * N + mask.to_bitmask().trailing_zeros() as usize);

            #[cfg(target_endian = "big")]
            return Some(i * N + mask.to_bitmask().leading_zeros() as usize);
        }
    }
    None
}

#[inline]
fn position_simd2<const N: usize>(haystack: &[Simd<u8, N>], b1: u8, b2: u8) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    let n1 = Simd::from_array([b1; N]);
    let n2 = Simd::from([b2; N]);
    for (i, &v) in haystack.iter().enumerate() {
        let mask = n1.simd_eq(v) | n2.simd_eq(v);
        if mask.any() {
            #[cfg(target_endian = "little")]
            return Some(i * N + mask.to_bitmask().trailing_zeros() as usize);

            #[cfg(target_endian = "big")]
            return Some(i * N + mask.to_bitmask().leading_zeros() as usize);
        }
    }
    None
}

#[cfg(test)]
mod test {
    use crate::search::bytes::position_by_byte;

    #[test]
    fn byte_position_test() {
        let inputs =[b"".to_vec(),b"00000000000100000000001".to_vec(), b"1000000000000000".to_vec(),
        b"00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001".to_vec()];

        for haystack in inputs {
            assert_eq!(
                position_by_byte::<4>(&haystack, b'1'),
                haystack.iter().position(|x| *x == b'1'),
                "result v. expected using {h}",
                h = String::from_utf8_lossy(&haystack)
            );
        }
    }
}

#[cfg(test)]
mod bench {
    use super::*;
    use std::sync::LazyLock;
    use test::Bencher;

    extern crate test;

    static LONG: &[u8] = b"0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001";
    static SHORT: LazyLock<Vec<u8>> = LazyLock::new(|| b"000000000001".to_vec());

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
        b.iter(|| position_by_byte::<16>(LONG, b'1'));
    }
}
