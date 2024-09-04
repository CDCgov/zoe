use std::{
    ops::Range,
    simd::{prelude::*, LaneCount, SupportedLaneCount},
};

use super::inexact::fuzzy_substring_match_simd;

///
/// Trait for searching byte substrings.
///
/// [`Nucleotides`](crate::data::types::nucleotides::Nucleotides),
/// [`AminoAcids`](crate::data::types::amino_acids::AminoAcids), and
/// [`QualityScores`](crate::data::types::phred::QualityScores) are printable
/// subsets of ASCII underneath the hood and therefore be efficiently searched
/// in this manner.
///
pub trait ByteSubstring {
    /// Returns `true` is the substring is found and `false` if not. Searches in
    /// the forward direction.
    fn contains_substring(&self, needle: impl AsRef<[u8]>) -> bool;
    /// Returns the substring's index range if found or [`None`] if not. Searches in
    /// the forward direction.
    fn find_substring(&self, needle: impl AsRef<[u8]>) -> Option<Range<usize>>;
    /// Finds the `needle` in the `haystack` using inexact matching up to the
    /// `DIFFERENCES_ALLOWED`. See [`fuzzy_substring_match_simd`] for
    /// more details.
    fn find_fuzzy_substring<const DIFFERENCES_ALLOWED: usize>(&self, needle: impl AsRef<[u8]>) -> Option<Range<usize>>;
}

impl<T: AsRef<[u8]> + ?Sized> ByteSubstring for T {
    #[inline]
    #[must_use]
    fn contains_substring(&self, needle: impl AsRef<[u8]>) -> bool {
        let haystack = self.as_ref();
        let needle = needle.as_ref();

        substring_match_simd::<32>(haystack, needle).is_some()
    }

    #[inline]
    #[must_use]
    fn find_substring(&self, needle: impl AsRef<[u8]>) -> Option<Range<usize>> {
        let haystack = self.as_ref();
        let needle = needle.as_ref();
        substring_match_simd::<32>(haystack, needle).map(|s| s..s + needle.len())
    }

    #[inline]
    #[must_use]
    fn find_fuzzy_substring<const DIFFERENCES_ALLOWED: usize>(&self, needle: impl AsRef<[u8]>) -> Option<Range<usize>> {
        let haystack = self.as_ref();
        let needle = needle.as_ref();

        fuzzy_substring_match_simd::<32, DIFFERENCES_ALLOWED>(haystack, needle).map(|s| s..s + needle.len())
    }
}

/// Similar to [`ByteSubstring`] but requires the input byte string to be
/// mutable.
pub trait ByteSubstringMut {
    /// If the substring is found, the first instance is removed and index range
    /// of the removed substring is returned, otherwise [`None`]. Searches in
    /// the forward direction.
    fn remove_first_substring(&mut self, needle: impl AsRef<[u8]>) -> Option<Range<usize>>;

    /// Replace a single byte of the stored sequence. Please see the
    /// *retain* and *recode* functions for a more wholistic approach.
    fn replace_all_bytes(&mut self, needle: u8, replacement: u8);
}

impl<T: AsMut<Vec<u8>> + AsRef<[u8]>> ByteSubstringMut for T {
    #[inline]
    fn remove_first_substring(&mut self, needle: impl AsRef<[u8]>) -> Option<Range<usize>> {
        if let Some(r) = self.find_substring(needle) {
            self.as_mut().drain(r.clone());
            Some(r)
        } else {
            None
        }
    }

    #[inline]
    fn replace_all_bytes(&mut self, needle: u8, replacement: u8) {
        crate::search::replace_all_bytes(self.as_mut(), needle, replacement);
    }
}

/// Finds the `needle` byte substring in the `haystack`, returning the starting index or [`None`] otherwise.
///
/// ### Limitations
///
/// This is a naïve exact match implementation and should only be used for small byte strings.
///
#[inline]
#[must_use]
pub fn substring_match(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.len() > haystack.len() || needle.is_empty() || haystack.is_empty() {
        return None;
    }

    for (i, w) in haystack.windows(needle.len()).enumerate() {
        if needle == w {
            return Some(i);
        }
    }

    None
}

/// Returns the starting index of the matched substring or [`None`] otherwise.
/// The const parameter `N` is used to specify the number of SIMD lanes for the
/// search.
///
/// ### Limitations
///
/// The current version is not optimized for larger needles (> 30 bp).
///
/// ### Citation
///
/// 1. Muła, Wojciech (2018). "SIMD-friendly algorithms for substring searching".
///    Available at:
///    <http://0x80.pl/articles/simd-strfind.html#algorithm-1-generic-simd>.
///    Accessed September 3, 2024.
#[inline]
#[must_use]
pub fn substring_match_simd<const N: usize>(haystack: &[u8], needle: &[u8]) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    if needle.len() > haystack.len() || needle.is_empty() || haystack.is_empty() {
        return None;
    }

    let first = needle[0];
    if needle.len() == 1 {
        return super::position_by_byte(haystack, first);
    }

    let Some((n2_offset, last)) = needle.iter().copied().enumerate().rev().find(|(_, b)| *b != first) else {
        return find_k_repeating(haystack, first, needle.len());
    };

    let n1 = Simd::from_array([first; N]);
    let n2 = Simd::from_array([last; N]);

    // In order to verify the needle, we need to subtract it off. However, the last character in the vector counts.
    let chunks1 = haystack[..=(haystack.len() - needle.len())]
        .chunks_exact(N)
        .map(Simd::from_slice);
    let chunks2 = haystack[n2_offset..].chunks_exact(N).map(Simd::from_slice);
    let z = std::iter::zip(chunks1, chunks2);

    let mut i = 0;
    for (c1, c2) in z {
        let f1 = n1.simd_eq(c1);
        let f2 = n2.simd_eq(c2);

        let mut m = (f1 & f2).to_bitmask();

        while m > 0 {
            let bit_position = m.trailing_zeros() as usize;
            let candidate_index = i + bit_position;

            if &haystack[candidate_index..candidate_index + needle.len()] == needle {
                return Some(candidate_index);
            }
            m &= m - 1;
        }
        i += N;
    }

    if N <= 8 {
        substring_match(&haystack[i..], needle).map(|j| i + j)
    } else {
        substring_match_simd::<8>(&haystack[i..], needle).map(|j| i + j)
    }
}

/// Given a haystack, needle and the number of times the needle repeats itself
/// in a row, find the starting index of the matching substring or return
/// [`None`] otherwise.
#[inline]
#[must_use]
pub fn find_k_repeating(haystack: &[u8], needle: u8, size: usize) -> Option<usize> {
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

/// Similar to [`find_k_repeating`] but takes a const `N` parameter for the
/// number of SIMD lanes. This algorithm is a simplified version of (1).
///
/// ### Limitations
///
/// Not optimized for needles >9bp.
///
/// ### Citation
///
/// 1. Tamanna Chhabra, Sukhpal Singh Ghuman, and Jorma Tarhio (2023). "String
///    Searching with Mismatches Using Avx2 and Avx-512 Instructions." doi:
///    <http://dx.doi.org/10.2139/ssrn.4662323>
#[inline]
#[must_use]
pub fn find_k_repeating_simd<const N: usize>(haystack: &[u8], needle: u8, size: usize) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    if size == 0 || size > haystack.len() {
        return None;
    } else if size == 1 {
        return super::bytes::position_by_byte(haystack, needle);
    }

    let minimum_simd_length = size + N - 1;
    if haystack.len() < minimum_simd_length {
        return find_k_repeating(haystack, needle, size);
    } else if size == haystack.len() {
        let (p, m, s) = haystack.as_simd();
        if p.iter().any(|b| *b != needle)
            || s.iter().any(|b| *b != needle)
            || super::bytes::position_simd::<N, 4, _>(m, |v| v.simd_ne(Simd::splat(needle))).is_some()
        {
            return None;
        }
        return Some(0);
    }

    let nv = Simd::from_array([needle; N]);
    let mut i = 0;

    'outer: while i < haystack.len() - minimum_simd_length {
        let h1 = Simd::from_slice(&haystack[i..]);
        let h2 = Simd::from_slice(&haystack[i + 1..]);
        let mut found_final = (h1.simd_eq(nv) & h2.simd_eq(nv)).to_bitmask();

        if found_final == 0 {
            i += N;
            continue 'outer;
        }

        for j in 2..size {
            let h = Simd::from_slice(&haystack[i + j..]);

            let bitmask = h.simd_eq(nv).to_bitmask();
            found_final &= bitmask;

            if found_final == 0 {
                i += N;
                continue 'outer;
            }
        }

        return Some(i + found_final.trailing_zeros() as usize);
    }

    find_k_repeating(&haystack[i..], needle, size).map(|j| i + j)
}

/// Similar to [`find_k_repeating_simd`] but modifies the search for needles of
/// size 3+ to use shift and count operations.
#[inline]
#[must_use]
pub fn find_k_repeating_simd2<const N: usize>(haystack: &[u8], needle: u8, size: usize) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    if size == 0 || size > haystack.len() {
        return None;
    } else if size == 1 {
        return super::bytes::position_by_byte(haystack, needle);
    }

    let nv = Simd::splat(needle);
    let chunks1 = haystack.chunks_exact(N);
    let chunks2 = haystack[1..].chunks_exact(N);
    let mut running_total = 0u32;
    let mut i = 0;

    for (c1, c2) in chunks1.map(Simd::from_slice).zip(chunks2.map(Simd::from_slice)) {
        let mut mask = (nv.simd_eq(c1) & nv.simd_eq(c2)).to_bitmask();

        if mask == 0 {
            i += N;
            running_total = 0;
            continue;
        }

        let reset = mask & (1 << (N - 1)) == 0;
        let mut bit_offset = 0;
        if mask & 1 > 0 {
            let previous = running_total;
            bit_offset = mask.trailing_ones();
            mask >>= bit_offset;
            running_total += bit_offset;

            if running_total as usize >= size - 1 {
                return Some(i - previous as usize);
            }
        }

        while mask > 0 {
            let skip_count = mask.trailing_zeros();
            mask >>= skip_count;
            bit_offset += skip_count;

            // Restart running total
            running_total = mask.trailing_ones();
            mask >>= running_total;

            if running_total as usize >= size - 1 {
                return Some(i + bit_offset as usize);
            }

            bit_offset += running_total;
        }

        if reset {
            running_total = 0;
        }

        i += N;
    }

    for (j, b) in haystack[i..].iter().copied().enumerate() {
        if b == needle {
            running_total += 1;
        } else {
            running_total = 0;
        }

        if running_total as usize == size {
            return Some(i + j - (size - 1));
        }
    }

    None
}

/// Based on an original algorithm by W. D. Chettleburgh but modifies the verification stage to include some SIMD steps.
///
/// ### Limitations
///
/// Performs much better for longer needles than short ones.
#[inline]
#[must_use]
pub fn find_k_repeating_simd3<const N: usize>(haystack: &[u8], needle: u8, size: usize) -> Option<usize>
where
    LaneCount<N>: SupportedLaneCount, {
    // Handle edge-cases for inputs
    if size == 0 || size > haystack.len() {
        return None;
    }

    let mut idx = size - 1;
    let mut left_limit = 0;

    'outer: while idx < haystack.len() && idx < size + N {
        if haystack[idx] == needle {
            for prev_idx in (left_limit..idx).rev() {
                if haystack[prev_idx] != needle {
                    left_limit = idx + 1;
                    idx = prev_idx + size;
                    continue 'outer;
                }
            }
            return Some(idx + 1 - size);
        }

        left_limit = idx + 1;
        idx += size;
    }

    let nv = Simd::splat(needle);

    'outer2: while idx < haystack.len() {
        if haystack[idx] == needle {
            let c = haystack[..idx].rchunks_exact(N).map(Simd::from_slice);

            let mut offset = 1;
            for v in c {
                let ones = v.reverse().simd_eq(nv).to_bitmask().trailing_ones();

                offset += ones as usize;
                if offset >= size {
                    return Some(idx + 1 - size);
                }

                if (ones as usize) < N {
                    let mismatch_index = idx - offset;
                    idx = mismatch_index + size;
                    continue 'outer2;
                }
            }
        } else {
            idx += size;
        }
    }

    None
}

#[cfg(test)]
mod test {
    use super::*;

    static PAD: &[u8; 150] = &[b'a'; 150];
    static NEEDLE: &[u8; 5] = b"hello";

    #[test]
    fn substring_match_units() {
        let mut haystack = *PAD;
        assert_eq!(None, substring_match(&haystack, NEEDLE));
        for start in 0..haystack.len() - NEEDLE.len() {
            haystack = *PAD;
            haystack[start..start + NEEDLE.len()].copy_from_slice(NEEDLE);
            assert_eq!(Some(start), substring_match(&haystack, NEEDLE));
        }
    }

    #[test]
    fn substring_match_simd_units() {
        let mut haystack = *PAD;
        assert_eq!(None, substring_match_simd::<8>(&haystack, NEEDLE));
        for start in 0..haystack.len() - NEEDLE.len() {
            haystack = *PAD;
            haystack[start..start + NEEDLE.len()].copy_from_slice(NEEDLE);
            assert_eq!(Some(start), substring_match_simd::<32>(&haystack, NEEDLE));
        }
    }

    #[test]
    fn substring_match_regressions() {
        let data = [
            (b"aaaaabaadaa".to_vec(), b"baabbbb".to_vec()),
            (b"dcxxxaxxxx".to_vec(), b"axx".to_vec()),
        ];

        for (haystack, needle) in data {
            assert_eq!(
                substring_match(&haystack, &needle),
                substring_match_simd::<16>(&haystack, &needle)
            );
        }
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
        ];
        for (haystack, needle, size) in data.into_iter().skip(7) {
            assert_eq!(
                find_k_repeating(&haystack, needle, size),
                find_k_repeating_simd::<8>(&haystack, needle, size)
            );
            assert_eq!(
                find_k_repeating(&haystack, needle, size),
                find_k_repeating_simd2::<8>(&haystack, needle, size)
            );
            assert_eq!(
                find_k_repeating(&haystack, needle, size),
                find_k_repeating_simd3::<8>(&haystack, needle, size)
            );
        }
    }

    #[test]
    fn find_k_repeating_units() {
        let mut s = b"aaaaaaab".repeat(15);
        s.extend(b"bb");

        assert_eq!(find_k_repeating(&s, b'b', 3), find_k_repeating_simd::<16>(&s, b'b', 3));
        assert_eq!(find_k_repeating(&s, b'b', 3), find_k_repeating_simd2::<16>(&s, b'b', 3));
        assert_eq!(find_k_repeating(&s, b'b', 3), find_k_repeating_simd3::<16>(&s, b'b', 3));
    }
}

#[cfg(test)]
mod bench {
    use super::*;
    use std::sync::LazyLock;
    use test::Bencher;

    extern crate test;

    static WORST_REPEATING: LazyLock<Vec<u8>> = LazyLock::new(|| {
        let mut s = b"bbabba".repeat(266);
        s.extend(b"abbb");
        s
    });

    static AVG_REPEATING: LazyLock<Vec<u8>> = LazyLock::new(|| {
        let mut s = b"aaaaaaab".repeat(199);
        s.extend(b"aaaaabbb");
        s
    });

    mod find_k_repeating {
        use super::*;

        mod average {
            use super::*;

            #[bench]
            fn using_substring_match(b: &mut Bencher) {
                let needle = vec![b'b'; 3];
                b.iter(|| substring_match_simd::<16>(&AVG_REPEATING, &needle));
            }

            #[bench]
            fn scalar(b: &mut Bencher) {
                b.iter(|| find_k_repeating(&AVG_REPEATING, b'b', 3));
            }

            #[bench]
            fn simd(b: &mut Bencher) {
                b.iter(|| find_k_repeating_simd::<16>(&AVG_REPEATING, b'b', 3));
            }

            #[bench]
            fn simd2(b: &mut Bencher) {
                b.iter(|| find_k_repeating_simd2::<16>(&AVG_REPEATING, b'b', 3));
            }

            #[bench]
            fn simd3(b: &mut Bencher) {
                b.iter(|| find_k_repeating_simd3::<16>(&AVG_REPEATING, b'b', 3));
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
            fn scalar(b: &mut Bencher) {
                b.iter(|| find_k_repeating(&WORST_REPEATING, b'b', 3));
            }

            #[bench]
            fn simd(b: &mut Bencher) {
                b.iter(|| find_k_repeating_simd::<16>(&WORST_REPEATING, b'b', 3));
            }

            #[bench]
            fn simd2(b: &mut Bencher) {
                b.iter(|| find_k_repeating_simd2::<16>(&WORST_REPEATING, b'b', 3));
            }

            #[bench]
            fn simd3(b: &mut Bencher) {
                b.iter(|| find_k_repeating_simd3::<16>(&WORST_REPEATING, b'b', 3));
            }
        }
    }
}
