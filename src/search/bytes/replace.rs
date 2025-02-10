use crate::simd::SimdByteFunctions;
use std::simd::{LaneCount, SupportedLaneCount};

/// Finds and replaces all instances of `needle` with the replacement byte.
#[inline]
pub fn replace_all_bytes<T>(v: &mut [T], needle: T, replacement: T)
where
    T: Copy + PartialEq, {
    for b in v.iter_mut() {
        if *b == needle {
            *b = replacement;
        }
    }
}

/// Finds and replace the specified byte value with another byte for some
/// mutable byte slice.
///
/// ## Limitations
///
/// Uses SIMD, but mostly useful just when the CPU target is below x86-64-v4,
/// otherwise the scalar code auto-vectorizes to the same performance (N=32).
/// For older targets, the fastest code may use N=16.
#[inline]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn replace_all_bytes_simd<const N: usize>(haystack: &mut [u8], needle: u8, replacement: u8)
where
    LaneCount<N>: SupportedLaneCount, {
    let (pre, mid, sfx) = haystack.as_simd_mut::<N>();
    replace_all_bytes(pre, needle, replacement);
    mid.iter_mut().for_each(|v| v.if_value_then_replace(needle, replacement));
    replace_all_bytes(sfx, needle, replacement);
}

#[cfg(test)]
mod bench {
    extern crate test;
    use super::*;
    use crate::generate::rand_sequence;
    use std::sync::LazyLock;
    use test::Bencher;

    static LONG: LazyLock<Vec<u8>> = LazyLock::new(|| rand_sequence(b"AGCT", 15_000, 42));
    static SHORT: LazyLock<Vec<u8>> = LazyLock::new(|| rand_sequence(b"AGCT", 150, 42));

    #[bench]
    fn find_replace_long_scalar(b: &mut Bencher) {
        b.iter(|| replace_all_bytes(&mut LONG.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_short_scalar(b: &mut Bencher) {
        b.iter(|| replace_all_bytes(&mut SHORT.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_long_simd32(b: &mut Bencher) {
        b.iter(|| replace_all_bytes_simd::<32>(&mut LONG.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_short_simd32(b: &mut Bencher) {
        b.iter(|| replace_all_bytes_simd::<32>(&mut SHORT.clone(), b'A', b'T'));
    }
}
