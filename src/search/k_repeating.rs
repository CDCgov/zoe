use crate::{iter_utils::SteppedWindows, search::position_by_byte_mapped_inner};
use std::simd::prelude::*;

/// Non-SIMD component of [`find_k_repeating`].
///
/// ## Limitations
///
/// Not optimized for small needles <21bp. Use [`find_k_repeating`] to more
/// appropriately handle this case.
#[inline]
#[must_use]
#[cfg(test)]
pub(super) fn find_k_repeating_scalar(haystack: &[u8], needle: u8, size: usize) -> Option<usize> {
    find_k_repeating_mapped_scalar(haystack, needle, size, |x| x)
}

/// Non-SIMD component of [`find_k_repeating_mapped`].
///
/// ## Limitations
///
/// Not optimized for small needles <21bp. Use [`find_k_repeating_mapped`] to
/// more appropriately handle this case.
#[inline]
#[must_use]
fn find_k_repeating_mapped_scalar<B>(haystack: &[u8], needle: u8, size: usize, map: B) -> Option<usize>
where
    B: Fn(u8) -> u8, {
    if size == 0 || size > haystack.len() {
        return None;
    }

    let mut idx = size - 1;
    let mut left_limit = 0usize;

    'outer: while idx < haystack.len() {
        if map(haystack[idx]) == needle {
            // Single match found: check backwards to see if it occurs `size` times
            for prev_idx in (left_limit..idx).rev() {
                if map(haystack[prev_idx]) != needle {
                    left_limit = idx + 1;
                    idx = prev_idx + size;
                    continue 'outer;
                }
            }

            // No mismatch found: return the starting index of the match
            return Some(idx + 1 - size);
        }

        // Mismatch found: advance idx to righmost position in leftmost possible match
        left_limit = idx + 1;
        idx += size;
    }

    None
}

/// This is the `size = 2` special case for the SIMD-acclerated portions of
/// [`find_k_repeating_mapped`].
#[inline]
#[must_use]
fn find_2_repeating_mapped_simd<const N: usize, S, B>(
    haystack: &[u8], needle: u8, simd_transform: S, byte_transform: B,
) -> Option<usize>
where
    S: Fn(Simd<u8, N>) -> Simd<u8, N>,
    B: Fn(u8) -> u8, {
    let nv = Simd::<u8, N>::splat(needle);

    // Get windows that are overlapping by 1
    let mut chunks = SteppedWindows::new(haystack, N - 1);

    for (i, chunk) in chunks.by_ref().copied().map(Simd::from_array).enumerate() {
        let mut mask = nv.simd_eq(simd_transform(chunk)).to_bitmask();
        mask = mask & (mask >> 1);

        if mask > 0 {
            return Some(i * (N - 1) + mask.trailing_zeros() as usize);
        }
    }

    let rem = chunks.final_partial_window();
    find_k_repeating_mapped_scalar(rem, needle, 2, byte_transform).map(|found| (haystack.len() - rem.len()) + found)
}

/// Given a haystack, needle and the number of times the needle repeats itself
/// in a row, find the starting index of the matching substring or return `None`
/// otherwise.
///
/// For small needle sizes (`size` < 21), SIMD acceleration is used.
///
/// For larger needle sizes (`size` ≥ 21), a scalar implementation is used,
/// enabling faster skipping through the string in a similar manner to the
/// Boyer-Moore algorithm.
///
/// ## Limitations
///
/// Not optimized for needles that are of comparable size to the haystack (e.g.,
/// over half the size of the haystack).
///
/// ## Parameters
///
/// `N` - The number of SIMD lanes to use for the search. This must be greater
/// than 2.
#[inline]
#[must_use]
#[allow(clippy::cast_possible_truncation)]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn find_k_repeating<const N: usize>(haystack: &[u8], needle: u8, size: usize) -> Option<usize> {
    find_k_repeating_inner::<N>(haystack, needle, size)
}

/// Given a haystack, needle and the number of times the needle repeats itself
/// in a row, find the starting index of the matching substring or return `None`
/// otherwise. The `haystack` is lazily mapped using `simd_transform` and
/// `byte_transform`.
///
/// For small needle sizes (`size` < 21), SIMD acceleration is used.
///
/// For larger needle sizes (`size` ≥ 21), a scalar implementation is used,
/// enabling faster skipping through the string in a similar manner to the
/// Boyer-Moore algorithm.
///
/// ## Limitations
///
/// Not optimized for needles that are of comparable size to the haystack (e.g.,
/// over half the size of the haystack).
///
/// ## Parameters
///
/// `N` - The number of SIMD lanes to use for the search. This must be greater
/// than 2.
#[inline]
#[must_use]
#[allow(clippy::cast_possible_truncation)]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn find_k_repeating_mapped<const N: usize, S, B>(
    haystack: &[u8], needle: u8, size: usize, simd_transform: S, byte_transform: B,
) -> Option<usize>
where
    S: Fn(Simd<u8, N>) -> Simd<u8, N>,
    B: Fn(u8) -> u8, {
    find_k_repeating_mapped_inner(haystack, needle, size, simd_transform, byte_transform)
}

/// Similar to [`find_k_repeating`], but without multiversioning (for use as a
/// helper function inside other multiversioned functions).
#[inline]
#[must_use]
pub(crate) fn find_k_repeating_inner<const N: usize>(haystack: &[u8], needle: u8, size: usize) -> Option<usize> {
    find_k_repeating_mapped_inner(haystack, needle, size, |v: Simd<u8, N>| v, |b| b)
}

/// Similar to [`find_k_repeating_mapped`], but without multiversioning (for use
/// as a helper function inside other multiversioned functions).
#[inline]
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub(crate) fn find_k_repeating_mapped_inner<const N: usize, S, B>(
    haystack: &[u8], needle: u8, size: usize, simd_transform: S, byte_transform: B,
) -> Option<usize>
where
    S: Fn(Simd<u8, N>) -> Simd<u8, N>,
    B: Fn(u8) -> u8, {
    const SCALAR_THRESHOLD: usize = 21;

    // TODO: fix when we have better const generics
    const {
        assert!(N > 2, "`find_k_repeating_simd` requires N > 2");
    };

    if size == 0 || size > haystack.len() {
        return None;
    } else if size == 1 {
        return position_by_byte_mapped_inner::<N, S, B>(haystack, needle, simd_transform, byte_transform);
    } else if size == 2 {
        return find_2_repeating_mapped_simd::<N, S, B>(haystack, needle, simd_transform, byte_transform);
    }

    // Fall back to Scalar algorithm when it isn't worth it
    if haystack.len() < N * 2 || size > SCALAR_THRESHOLD {
        return find_k_repeating_mapped_scalar(haystack, needle, size, byte_transform);
    }

    let nv = Simd::<u8, N>::splat(needle);
    // Get windows that are overlapping by 2
    let mut chunks = SteppedWindows::new(haystack, N - 2);
    let mut running_total = 0;

    for (i, chunk) in chunks.by_ref().copied().map(Simd::from_array).enumerate() {
        let mut mask = nv.simd_eq(simd_transform(chunk)).to_bitmask();
        mask = mask & (mask >> 1) & (mask >> 2);

        if mask == 0 {
            running_total = 0;
            continue;
        }

        let reset = mask & (1 << (N - 3)) == 0;
        let mut bit_offset = 0;
        if mask & 1 > 0 {
            let previous = running_total;
            bit_offset = mask.trailing_ones();
            mask >>= bit_offset;

            // CORRECTNESS: `trailing_ones()` cannot be larger than 64, which is
            // much smaller than the `usize::MAX` of hypothetical 16-bit `usize`
            // even though the value returned uses a `u32`.
            running_total += bit_offset as usize;

            if running_total >= size - 2 {
                return Some(i * (N - 2) - previous);
            }
        }

        while mask > 0 {
            let skip_count = mask.trailing_zeros();
            mask >>= skip_count;
            bit_offset += skip_count;

            // Restart running total
            running_total = mask.trailing_ones() as usize;
            mask >>= running_total;

            if running_total >= size - 2 {
                return Some(i * (N - 2) + bit_offset as usize);
            }

            bit_offset += running_total as u32;
        }

        if reset {
            running_total = 0;
        }
    }

    for (j, b) in chunks.final_partial_window().iter().copied().enumerate() {
        if byte_transform(b) == needle {
            running_total += 1;
        } else {
            running_total = 0;
        }

        if running_total == size {
            return Some(haystack.len() - chunks.final_partial_window().len() + j - (size - 1));
        }
    }

    None
}

#[cfg(test)]
mod test {
    use super::*;

    #[inline]
    #[must_use]
    pub fn find_k_repeating_naive(haystack: &[u8], needle: u8, size: usize) -> Option<usize> {
        let mut hits = 0usize;
        let offset = if size > 0 { size - 1 } else { return None };

        for (index, h) in haystack.iter().copied().enumerate() {
            if needle == h {
                hits += 1;
            } else {
                hits = 0;
            }

            if hits == size {
                return Some(index - offset);
            }
        }

        None
    }

    #[test]
    fn find_k_repeating_regressions() {
        let data = [
            (b"0000".to_vec(), b'0', 4),
            (b"caaaaaab".to_vec(), b'a', 6),
            (b"aaaaaaabaa".to_vec(), b'a', 8),
            (b"aabcaaaaa".to_vec(), b'a', 4),
            (b"1100000001".to_vec(), b'0', 3),
            (b"88000005000000".to_vec(), b'0', 8),
            (b"0000000000".to_vec(), b'0', 10),
            (b"AMMMMMMMMMMMMMMMMMMMMMMOMMML".to_vec(), b'M', 12),
            (b"aaaaaaabbbbbbbbbbbb".to_vec(), b'b', 11),
            (b"aababbabbbaaabbbbbbbbbbbbbbbbbbbbbbbbb".to_vec(), b'b', 24),
            (b"abbbbbbbbbbbbbbbbbbbbbbbaaaabbbbbbbbbb".to_vec(), b'b', 24),
            (b"abababababababababababb".to_vec(), b'b', 2),
        ];
        for (haystack, needle, size) in data {
            assert_eq!(
                find_k_repeating_naive(&haystack, needle, size),
                find_k_repeating::<8>(&haystack, needle, size)
            );
            assert_eq!(
                find_k_repeating_naive(&haystack, needle, size),
                find_k_repeating_scalar(&haystack, needle, size)
            );
            assert_eq!(
                find_k_repeating_naive(&haystack, needle, size),
                find_k_repeating::<8>(&haystack, needle, size)
            );
        }
    }

    #[test]
    fn find_k_repeating_units() {
        let mut s = b"aaaaaaab".repeat(15);
        s.extend(b"bb");

        assert_eq!(find_k_repeating_naive(&s, b'b', 3), find_k_repeating::<16>(&s, b'b', 3));
        assert_eq!(find_k_repeating_naive(&s, b'b', 3), find_k_repeating_scalar(&s, b'b', 3));
    }
}
