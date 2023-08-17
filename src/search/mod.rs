#![allow(clippy::comparison_chain, clippy::explicit_iter_loop)]

use super::simd::SimdByteFunctions;
use std::ops::Range;
use std::simd::{LaneCount, SupportedLaneCount};

/// Finds and replaces the specified value with another for some mutable slice.
#[inline]
pub(crate) fn find_and_replace<T>(v: &mut [T], needle: T, replacement: T)
where
    T: Copy + PartialEq, {
    for b in v.iter_mut() {
        if *b == needle {
            *b = replacement;
        }
    }
}

#[allow(dead_code)]
/// Finds and replace the specified byte value with another byte for some
/// mutable byte slice. Uses SIMD, but mostly useful just when the CPU target is
/// below x86-64-v4, otherwise the scalar code auto-vectorizes better.
#[inline]
pub(crate) fn find_and_replace_simd<const N: usize>(haystack: &mut [u8], needle: u8, replacement: u8)
where
    LaneCount<N>: SupportedLaneCount, {
    let (pre, mid, sfx) = haystack.as_simd_mut::<N>();
    find_and_replace(pre, needle, replacement);
    mid.iter_mut().for_each(|v| v.if_value_then_replace(needle, replacement));
    find_and_replace(sfx, needle, replacement);
}

pub trait Subsequence {
    fn contains_subsequence(&self, needle: impl AsRef<[u8]>) -> bool;
    fn find_subsequence(&self, needle: impl AsRef<[u8]>) -> Option<Range<usize>>;
}

impl<T: AsRef<[u8]> + ?Sized> Subsequence for T {
    fn contains_subsequence(&self, needle: impl AsRef<[u8]>) -> bool {
        let haystack = self.as_ref();
        let needle = needle.as_ref();

        if needle.len() < haystack.len() {
            if needle.is_empty() {
                return false;
            }

            // TO-DO: improve with SIMD or memmem
            haystack.windows(needle.len()).any(|w| w == needle)
        } else if needle.len() > haystack.len() {
            false
        } else {
            // Implies: needle.len() == haystack.len()
            haystack == needle
        }
    }

    fn find_subsequence(&self, needle: impl AsRef<[u8]>) -> Option<Range<usize>> {
        let haystack = self.as_ref();
        let needle = needle.as_ref();

        if needle.len() < haystack.len() {
            if needle.is_empty() {
                return None;
            }

            // TO-DO: improve with SIMD or memmem
            for (start, w) in haystack.windows(needle.len()).enumerate() {
                if w == needle {
                    return Some(start..(start + needle.len()));
                }
            }
            None
        } else if needle.len() > haystack.len() {
            None
        } else {
            // Implies: needle.len() == haystack.len()
            if haystack == needle {
                Some(0..haystack.len())
            } else {
                None
            }
        }
    }
}

pub trait VectorSubsequence {
    fn remove_first_subsequence(&mut self, needle: impl AsRef<[u8]>) -> Option<Range<usize>>;
}

impl<T: AsMut<Vec<u8>> + AsRef<[u8]>> VectorSubsequence for T {
    fn remove_first_subsequence(&mut self, needle: impl AsRef<[u8]>) -> Option<Range<usize>> {
        if let Some(r) = self.find_subsequence(needle) {
            self.as_mut().drain(r.clone());
            Some(r)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod bench {
    extern crate test;
    use super::*;
    use crate::generate::rand_sequence;
    use lazy_static::lazy_static;
    use test::Bencher;

    lazy_static! {
        static ref LONG: Vec<u8> = rand_sequence(b"AGCT", 15_000, 42);
    }
    lazy_static! {
        static ref SHORT: Vec<u8> = rand_sequence(b"AGCT", 150, 42);
    }

    #[bench]
    fn find_replace_long_scalar(b: &mut Bencher) {
        b.iter(|| find_and_replace(&mut LONG.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_short_scalar(b: &mut Bencher) {
        b.iter(|| find_and_replace(&mut SHORT.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_long_simd16(b: &mut Bencher) {
        b.iter(|| find_and_replace_simd::<16>(&mut LONG.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_short_simd16(b: &mut Bencher) {
        b.iter(|| find_and_replace_simd::<16>(&mut SHORT.clone(), b'A', b'T'));
    }
}
