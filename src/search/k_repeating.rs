use crate::iter_utils::SteppedWindows;
use std::simd::{prelude::*, LaneCount, SupportedLaneCount};

/// Non-SIMD component of [`find_k_repeating`].
///
/// ### Limitations
///
/// Not optimized for small needles <21bp. Use [`find_k_repeating`] to more
/// appropriately handle this case.
#[inline]
#[must_use]
fn find_k_repeating_scalar(haystack: &[u8], needle: u8, size: usize) -> Option<usize> {
    if size == 0 || size > haystack.len() {
        return None;
    }

    let mut idx = size - 1;
    let mut left_limit = 0usize;

    'outer: while idx < haystack.len() {
        if haystack[idx] == needle {
            // Single match found: check backwards to see if it occurs `size` times
            for prev_idx in (left_limit..idx).rev() {
                if haystack[prev_idx] != needle {
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
/// [`find_k_repeating`].
#[inline]
#[must_use]
fn find_2_repeating_simd<const N: usize>(haystack: &[u8], needle: u8) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    let nv = Simd::splat(needle);
    let mut chunks = SteppedWindows::new(haystack, N - 1, N);

    for (i, chunk) in chunks.by_ref().map(Simd::from_slice).enumerate() {
        let mut mask = nv.simd_eq(chunk).to_bitmask();
        mask = mask & (mask >> 1);

        if mask > 0 {
            return Some(i * (N - 1) + mask.trailing_zeros() as usize);
        }
    }

    let rem = chunks.remainder();
    find_k_repeating_scalar(rem, needle, 2).map(|found| (haystack.len() - rem.len()) + found)
}

/// Given a haystack, needle and the number of times the needle repeats itself
/// in a row, find the starting index of the matching substring or return
/// `None` otherwise.
///
/// For small needle sizes (`size` < 21), SIMD acceleration is used. The
/// const parameter `N` is used to specify the number of SIMD lanes for the
/// search.
///
/// For larger needle sizes (`size` ≥ 21), a scalar implementation
/// is used, enabling faster skipping through the string in a similar manner
/// to the Boyer-Moore algorithm.
///
/// ### Limitations
///
/// Not optimized for needles that are of comparable size to the haystack (e.g.,
/// over half the size of the haystack).
///
#[allow(clippy::cast_possible_truncation)]
#[inline]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn find_k_repeating<const N: usize>(haystack: &[u8], needle: u8, size: usize) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    const SCALAR_THRESHOLD: usize = 21;

    // TODO: fix when we have better const generics
    const {
        assert!(N > 2, "`find_k_repeating_simd` requires N > 2");
    };

    if size == 0 || size > haystack.len() {
        return None;
    } else if size == 1 {
        return super::bytes::position_by_byte(haystack, needle);
    } else if size == 2 {
        return find_2_repeating_simd(haystack, needle);
    }

    // Fall back to Scalar algorithm when it isn't worth it
    if haystack.len() < N * 2 || size > SCALAR_THRESHOLD {
        return find_k_repeating_scalar(haystack, needle, size);
    }

    let nv = Simd::splat(needle);
    let mut chunks = SteppedWindows::new(haystack, N - 2, N);
    let mut running_total = 0;

    for (i, chunk) in chunks.by_ref().map(Simd::from_slice).enumerate() {
        let mut mask = nv.simd_eq(chunk).to_bitmask();
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

    for (j, b) in chunks.remainder().iter().copied().enumerate() {
        if b == needle {
            running_total += 1;
        } else {
            running_total = 0;
        }

        if running_total == size {
            return Some(haystack.len() - chunks.remainder().len() + j - (size - 1));
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

#[cfg(test)]
mod bench {
    use super::super::substring::substring_match_simd;
    use super::*;
    use std::sync::LazyLock;
    use test::Bencher;

    extern crate test;

    static WORST_REPEATING: LazyLock<Vec<u8>> = LazyLock::new(|| {
        let mut s = b"abb".to_vec();
        s.extend(b"aabb".repeat(398));
        s.extend(b"aabbb");
        s
    });

    static AVG_REPEATING: LazyLock<Vec<u8>> = LazyLock::new(|| {
        let mut s = b"aaaaaaab".repeat(199);
        s.extend(b"aaaaabbb");
        s
    });

    mod average {
        use super::*;

        #[bench]
        fn using_substring_match(b: &mut Bencher) {
            let needle = vec![b'b'; 3];
            b.iter(|| substring_match_simd::<16>(&AVG_REPEATING, &needle));
        }

        #[bench]
        fn composite(b: &mut Bencher) {
            b.iter(|| find_k_repeating::<16>(&AVG_REPEATING, b'b', 3));
        }

        #[bench]
        fn simd(b: &mut Bencher) {
            b.iter(|| find_k_repeating::<16>(&AVG_REPEATING, b'b', 3));
        }

        #[bench]
        fn scalar(b: &mut Bencher) {
            b.iter(|| find_k_repeating_scalar(&AVG_REPEATING, b'b', 3));
        }
    }

    mod worst {
        use super::*;

        #[bench]
        fn using_substring_match(b: &mut Bencher) {
            let needle = vec![b'b'; 3];
            b.iter(|| substring_match_simd::<16>(&WORST_REPEATING, &needle));
        }

        #[bench]
        fn composite(b: &mut Bencher) {
            b.iter(|| find_k_repeating::<16>(&WORST_REPEATING, b'b', 3));
        }

        #[bench]
        fn simd(b: &mut Bencher) {
            b.iter(|| find_k_repeating::<16>(&WORST_REPEATING, b'b', 3));
        }

        #[bench]
        fn scalar(b: &mut Bencher) {
            b.iter(|| find_k_repeating_scalar(&WORST_REPEATING, b'b', 3));
        }
    }
}
