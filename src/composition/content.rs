use crate::simd::SimdByteFunctions;
use std::simd::{prelude::*, LaneCount, SupportedLaneCount};

#[allow(non_snake_case, clippy::cast_abs_to_unsigned)]
pub(crate) fn gc_content_simd<const N: usize>(s: &[u8]) -> usize
where
    LaneCount<N>: SupportedLaneCount, {
    let G = Simd::splat(b'G');
    let C = Simd::splat(b'C');

    let (pre, mid, sfx) = s.as_simd::<N>();

    let count: usize = mid
        .iter()
        .map(|&w| {
            let v = w.to_ascii_uppercase();
            (v.simd_eq(G) | v.simd_eq(C)).to_int().reduce_sum().abs() as usize
        })
        .sum();

    count + gc_content(pre) + gc_content(sfx)
}

pub(crate) fn gc_content(s: &[u8]) -> usize {
    s.iter()
        .map(|&b| b.to_ascii_uppercase())
        .filter(|&b| b == b'G' || b == b'C')
        .count()
}
