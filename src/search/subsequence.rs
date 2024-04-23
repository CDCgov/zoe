#![allow(clippy::comparison_chain)]

use std::ops::Range;

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
