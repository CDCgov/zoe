use crate::simd::SimdByteFunctions;

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
pub fn replace_all_bytes_simd<const N: usize>(haystack: &mut [u8], needle: u8, replacement: u8) {
    let (pre, mid, sfx) = haystack.as_simd_mut::<N>();
    replace_all_bytes(pre, needle, replacement);
    for v in mid.iter_mut() {
        v.if_value_then_replace(needle, replacement);
    }
    replace_all_bytes(sfx, needle, replacement);
}
