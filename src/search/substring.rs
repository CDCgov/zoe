use crate::search::RangeSearch;

use super::{inexact::fuzzy_substring_match_simd, k_repeating::find_k_repeating};
use std::{
    ops::Range,
    simd::{LaneCount, SupportedLaneCount, prelude::*},
};

/// Trait for searching byte substrings.
///
/// [`Nucleotides`](crate::data::types::nucleotides::Nucleotides),
/// [`AminoAcids`](crate::data::types::amino_acids::AminoAcids), and
/// [`QualityScores`](crate::data::types::phred::QualityScores) are printable
/// subsets of ASCII underneath the hood and therefore be efficiently searched
/// in this manner.
///
/// This trait is also compatible with [`RangeSearch`]. See [Restricting the
/// search range](crate::search#restricting-the-search-range) for more details.
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
    /// Find contiguous, repeating byte characters in a byte slice. See [`find_k_repeating`]
    /// for more details.
    fn find_repeating_byte(&self, needle: u8, size: usize) -> Option<Range<usize>>;
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

    #[inline]
    #[must_use]
    fn find_repeating_byte(&self, needle: u8, size: usize) -> Option<Range<usize>> {
        let haystack = self.as_ref();

        find_k_repeating::<32>(haystack, needle, size).map(|s| s..s + size)
    }
}

impl ByteSubstring for RangeSearch<'_> {
    #[inline]
    fn contains_substring(&self, needle: impl AsRef<[u8]>) -> bool {
        self.slice.contains_substring(needle)
    }

    #[inline]
    fn find_substring(&self, needle: impl AsRef<[u8]>) -> Option<Range<usize>> {
        self.slice.find_substring(needle).map(|r| self.adjust_to_context(&r))
    }

    #[inline]
    fn find_fuzzy_substring<const DIFFERENCES_ALLOWED: usize>(&self, needle: impl AsRef<[u8]>) -> Option<Range<usize>> {
        self.slice
            .find_fuzzy_substring::<DIFFERENCES_ALLOWED>(needle)
            .map(|r| self.adjust_to_context(&r))
    }

    #[inline]
    fn find_repeating_byte(&self, needle: u8, size: usize) -> Option<Range<usize>> {
        self.slice
            .find_repeating_byte(needle, size)
            .map(|r| self.adjust_to_context(&r))
    }
}

/// Similar to [`ByteSubstring`] but requires the input byte string to be
/// mutable.
pub trait ByteSubstringMut {
    /// If the substring is found, the first instance is removed and index range
    /// of the removed substring is returned, otherwise [`None`]. Searches in
    /// the forward direction.
    fn remove_first_substring(&mut self, needle: impl AsRef<[u8]>) -> Option<Range<usize>>;

    /// Replace a single byte of the stored sequence. See the [`Recode`]
    /// trait for a more wholistic approach.
    ///
    /// [`Recode`]: crate::data::Recode
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

/// Finds the `needle` byte substring in the `haystack`, returning the starting
/// index or [`None`] otherwise.
///
/// ### Limitations
///
/// This is a naïve exact match implementation and should only be used for small
/// byte strings.
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
/// 1. Muła, Wojciech (2018). "SIMD-friendly algorithms for substring searching."
///    Available at:
///    <http://0x80.pl/articles/simd-strfind.html#algorithm-1-generic-simd>.
///    Accessed September 3, 2024.
#[inline]
#[must_use]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
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

    // In order to verify the needle, we need to subtract it off. However, the
    // last character in the vector counts.
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
}
