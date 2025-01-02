use crate::simd::SimdByteFunctions;
use std::simd::{prelude::*, LaneCount, SupportedLaneCount};

/// Calculates GC content using SIMD operations for better performance. Returns
/// count of G and C nucleotides in the sequence.
#[must_use]
pub fn gc_content_simd<const N: usize>(s: &[u8]) -> usize
where
    LaneCount<N>: SupportedLaneCount, {
    let g = Simd::splat(b'G');
    let c = Simd::splat(b'C');

    let (pre, mid, sfx) = s.as_simd::<N>();

    let count: usize = mid
        .iter()
        .map(|&w| {
            let v = w.to_ascii_uppercase();
            (v.simd_eq(g) | v.simd_eq(c)).to_bitmask().count_ones() as usize
        })
        .sum();

    count + gc_content(pre) + gc_content(sfx)
}

/// Calculates GC content using scalar operations. Returns count of G and C
/// nucleotides in the sequence.
#[must_use]
pub fn gc_content(s: &[u8]) -> usize {
    s.iter()
        .map(|&b| b.to_ascii_uppercase())
        .filter(|&b| b == b'G' || b == b'C')
        .count()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_gc_content() {
        for length in [16, 1200, 10000] {
            let s = crate::generate::rand_sequence(b"ATCG", length, 42);
            let (scalar, simd) = (gc_content(&s), gc_content_simd::<16>(&s));
            assert_eq!(scalar, simd, "when testing for length {length}");
        }
    }
}

#[cfg(test)]
mod bench {
    use super::*;
    use std::sync::LazyLock;
    use test::Bencher;

    extern crate test;

    const N: usize = 1600;
    static SEQ: LazyLock<Vec<u8>> = std::sync::LazyLock::new(|| crate::generate::rand_sequence(b"ATCG", N, 42));

    #[bench]
    fn gc_content_scalar(b: &mut Bencher) {
        b.iter(|| gc_content(&SEQ));
    }

    #[bench]
    fn gc_content_simd16(b: &mut Bencher) {
        b.iter(|| gc_content_simd::<16>(&SEQ));
    }
}
