use std::{
    iter::zip,
    simd::{prelude::*, LaneCount, SupportedLaneCount},
};

/// Provides the starting index of a matched "fuzzy substring" with however many
/// `differences_allowed` from the needle, otherwise [`None`] is returned.
///
/// ### Limitations
///
/// The naïve algorithm should not be used except for small haystacks.
///
#[must_use]
#[inline]
pub fn fuzzy_substring_match(haystack: &[u8], needle: &[u8], differences_allowed: u8) -> Option<usize> {
    if needle.len() > haystack.len() || needle.is_empty() || haystack.is_empty() {
        return None;
    }

    if differences_allowed as usize >= needle.len() {
        return Some(0);
    }

    fuzzy_substring_match_scalar(haystack, needle, differences_allowed)
}

/// A naïve fuzzy substring matching algorithm.
#[inline]
#[must_use]
pub(crate) fn fuzzy_substring_match_scalar(haystack: &[u8], needle: &[u8], differences_allowed: u8) -> Option<usize> {
    for (i, w) in haystack.windows(needle.len()).enumerate() {
        let mut differences = 0;
        for (&h, &n) in zip(w, needle) {
            if h != n {
                differences += 1;
                if differences > differences_allowed {
                    break;
                }
            }
        }

        if differences <= differences_allowed {
            return Some(i);
        }
    }

    None
}

/// Similar to [`fuzzy_substring_match`] but takes a const `N` parameter for the
/// number of SIMD lanes and a second const parameter `K` for the number of
/// differences allowed.
///
/// This algorithm implements a portable version of (1).
///
/// ### Limitations
///
/// Not optimized for either large needles or large differences allowed (`K`).
/// However, it should be quite suitable for short oligonucleotide searches with
/// 1 or 2 mismatches.
///
/// ### Citation
///
/// 1. Tamanna Chhabra, Sukhpal Singh Ghuman, and Jorma Tarhio (2023). "String
///    Searching with Mismatches Using Avx2 and Avx-512 Instructions." doi:
///    <http://dx.doi.org/10.2139/ssrn.4662323>
///
#[allow(clippy::cast_possible_truncation)]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn fuzzy_substring_match_simd<const N: usize, const K: usize>(haystack: &[u8], needle: &[u8]) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    // Scalar version takes u8
    const {
        assert!(K <= u8::MAX as usize);
    };

    if needle.len() > haystack.len() || needle.is_empty() || haystack.is_empty() {
        return None;
    }

    // Defer to the exact algorithm as needed.
    if K == 0 {
        return super::substring::substring_match_simd(haystack, needle);
    }

    // Since `needle.len() <= haystack.len()` we can just grab the first possible
    // match for trivial case.
    if K >= needle.len() {
        return Some(0);
    }

    // For short needles, this is slower. I suspect it has to do with branch
    // misprediction vs. the vanilla version.
    let minimum_simd_length = needle.len() + N - 1;
    if haystack.len() < std::cmp::max(2 * N, minimum_simd_length) {
        return fuzzy_substring_match_scalar(haystack, needle, K as u8);
    }

    let needles: Vec<_> = needle.iter().copied().map(|b| Simd::from_array([b; N])).collect();
    let mut i = 0;
    let (left, right) = needles.split_at(K);

    'outer: while i < haystack.len() - minimum_simd_length {
        // This is necessary because of incomplete `const_generic_exprs`
        let mut found = [u64::MAX; K];
        let mut found_final = u64::MAX;

        for (j, &b) in left.iter().enumerate() {
            let h = Simd::from_slice(&haystack[i + j..]);
            let c = h.simd_eq(b).to_bitmask();

            for k in (1..j).rev() {
                // j < left.len() < K
                found[k] &= found[k - 1] | c;
            }
            found[0] &= found_final | c;
            found_final &= c;
        }

        if found[K - 1] == 0 {
            i += N;
            continue 'outer;
        }

        for (j, &b) in right.iter().enumerate().map(|(j, v)| (j + K, v)) {
            let h = Simd::from_slice(&haystack[i + j..]);
            let c = h.simd_eq(b).to_bitmask();

            for k in (1..K).rev() {
                found[k] &= found[k - 1] | c;
            }
            found[0] &= found_final | c;
            found_final &= c;

            if found[K - 1] == 0 {
                i += N;
                continue 'outer;
            }
        }

        return Some(i + found[K - 1].trailing_zeros() as usize);
    }

    fuzzy_substring_match_scalar(&haystack[i..], needle, K as u8).map(|j| i + j)
}

#[cfg(test)]
mod test {
    #![allow(clippy::type_complexity)]

    use super::*;

    static DATA: [(&[u8], &[u8], u8, Option<usize>); 9] = [
        (b"aaaaaaaaaaaaaaaaaa", b"aaa", 0, Some(0)),
        (b"aaaaaaaaaaaaaaabbb", b"bbb", 0, Some(15)),
        (b"aaaaaabbbaaaaaaaaa", b"bbb", 0, Some(6)),
        (b"aaaaaababaaaaaaaaa", b"bbb", 1, Some(6)),
        (b"bababababababababa", b"bbb", 1, Some(0)),
        (b"aaaaaababaaaaaaaaa", b"bbb", 2, Some(4)),
        (b"aaaaaaaaaaaaaaaaaa", b"bbb", 2, None),
        (b"aaa", b"aaaaaaaaaaaaaaaaaa", 3, None),
        (b"aaaaaaaaaaaaaaaaaa", b"bbb", 3, Some(0)),
    ];

    #[test]
    fn test_fuzzy_naive() {
        for (haystack, needle, d, expected) in DATA {
            assert_eq!(fuzzy_substring_match(haystack, needle, d), expected);
        }
    }

    #[test]
    pub fn test_fuzzy_simd3() {
        for (haystack, needle, d, expected) in DATA {
            let computed = match d {
                1 => fuzzy_substring_match_simd::<4, 1>(haystack, needle),
                2 => fuzzy_substring_match_simd::<4, 2>(haystack, needle),
                3 => fuzzy_substring_match_simd::<4, 3>(haystack, needle),
                _ => fuzzy_substring_match_simd::<4, 0>(haystack, needle),
            };

            assert_eq!(
                computed,
                expected,
                "haystack: {}, needle: {} d: {}",
                String::from_utf8_lossy(haystack),
                String::from_utf8_lossy(needle),
                d
            );
        }
    }

    #[test]
    pub fn test_regression() {
        let data = [(b"ca".to_vec(), b"ab".to_vec()), (b"aa".to_vec(), b"ba".to_vec())];

        for (haystack, needle) in data {
            assert_eq!(
                fuzzy_substring_match(&haystack, &needle, 1),
                fuzzy_substring_match_simd::<4, 1>(&haystack, &needle)
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

    static L: LazyLock<Vec<u8>> = LazyLock::new(|| {
        let mut s = b"aaaaaaaab".repeat(15);
        s.extend(b"bbbbbbbbbbbbbbbb");
        s
    });

    static S: LazyLock<Vec<u8>> = LazyLock::new(|| {
        let mut s = b"aaaaaaaab".repeat(2);
        s.extend(b"bbb");
        s
    });

    mod short {
        use super::*;

        fn wrapper(haystack: &[u8], needle: &[u8], differences_allowed: u8) -> Option<usize> {
            fuzzy_substring_match(haystack, needle, differences_allowed)
        }

        #[bench]
        fn fuzzy_match_naive_scalar(b: &mut Bencher) {
            b.iter(|| fuzzy_substring_match(&S, b"bbb".as_slice(), 2));
        }

        #[bench]
        fn fuzzy_match_simd3(b: &mut Bencher) {
            b.iter(|| fuzzy_substring_match_simd::<16, 2>(&S, b"bbb".as_slice()));
        }

        #[bench]
        fn fuzzy_match_naive_scalar2(b: &mut Bencher) {
            b.iter(|| wrapper(&S, b"bbb".as_slice(), 2));
        }
    }

    mod long {
        use super::*;

        #[bench]
        fn fuzzy_match_naive_scalar(b: &mut Bencher) {
            b.iter(|| fuzzy_substring_match(&L, b"bbbbbbbbbbbbbbbb".as_slice(), 2));
        }

        #[bench]
        fn fuzzy_match_simd3(b: &mut Bencher) {
            b.iter(|| fuzzy_substring_match_simd::<16, 2>(&L, b"bbbbbbbbbbbbbbbb".as_slice()));
        }
    }
}
