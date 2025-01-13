use crate::{data::validation::CheckSequence, prelude::*};
use std::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

impl AsRef<[u8]> for Nucleotides {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsRef<[u8]> for NucleotidesView<'_> {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl AsRef<[u8]> for NucleotidesViewMut<'_> {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl AsMut<[u8]> for Nucleotides {
    #[inline]
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl AsMut<[u8]> for NucleotidesViewMut<'_> {
    #[inline]
    fn as_mut(&mut self) -> &mut [u8] {
        self.0
    }
}

impl AsRef<Vec<u8>> for Nucleotides {
    #[inline]
    fn as_ref(&self) -> &Vec<u8> {
        &self.0
    }
}

impl AsMut<Vec<u8>> for Nucleotides {
    #[inline]
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }
}

impl From<String> for Nucleotides {
    #[inline]
    fn from(s: String) -> Self {
        Nucleotides(s.as_bytes().to_vec())
    }
}

impl From<Vec<u8>> for Nucleotides {
    #[inline]
    fn from(vec: Vec<u8>) -> Self {
        Nucleotides(vec)
    }
}

impl From<&[u8]> for Nucleotides {
    #[inline]
    fn from(bytes: &[u8]) -> Self {
        Nucleotides(bytes.to_vec())
    }
}

impl<'a> From<&'a [u8]> for NucleotidesView<'a> {
    #[inline]
    fn from(bytes: &'a [u8]) -> Self {
        NucleotidesView(bytes)
    }
}

// TODO: Implement From &mut [u8] for NucleotidesViewMut?? (and other seq-like types)

impl<const N: usize> From<&[u8; N]> for Nucleotides {
    #[inline]
    fn from(bytes: &[u8; N]) -> Self {
        Nucleotides(bytes.to_vec())
    }
}

impl<'a, const N: usize> From<&'a [u8; N]> for NucleotidesView<'a> {
    #[inline]
    fn from(bytes: &'a [u8; N]) -> Self {
        NucleotidesView(bytes)
    }
}

impl IntoIterator for Nucleotides {
    type Item = u8;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for NucleotidesView<'a> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for NucleotidesViewMut<'a> {
    type Item = &'a mut u8;
    type IntoIter = std::slice::IterMut<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<'a> IntoIterator for &'a Nucleotides {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a NucleotidesView<'_> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a NucleotidesViewMut<'_> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a mut Nucleotides {
    type Item = &'a mut u8;
    type IntoIter = std::slice::IterMut<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<'a> IntoIterator for &'a mut NucleotidesViewMut<'_> {
    type Item = &'a mut u8;
    type IntoIter = std::slice::IterMut<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl FromIterator<u8> for Nucleotides {
    #[inline]
    fn from_iter<T: IntoIterator<Item = u8>>(iterable: T) -> Self {
        Nucleotides(iterable.into_iter().collect())
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for Nucleotides {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for NucleotidesView<'_> {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for NucleotidesViewMut<'_> {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> IndexMut<I> for Nucleotides {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> IndexMut<I> for NucleotidesViewMut<'_> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl std::fmt::Display for Nucleotides {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.0.is_ascii_simd::<16>() {
            // SAFETY: we just checked it is ASCII using our fast SIMD function.
            // ASCII is valid UTF8.
            write!(f, "{}", unsafe { std::str::from_utf8_unchecked(&self.0) })
        } else {
            write!(f, "{}", String::from_utf8_lossy(&self.0))
        }
    }
}

impl std::fmt::Display for NucleotidesView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.0.is_ascii_simd::<16>() {
            // SAFETY: we just checked it is ASCII using our fast SIMD function.
            // ASCII is valid UTF8.
            write!(f, "{}", unsafe { std::str::from_utf8_unchecked(self.0) })
        } else {
            write!(f, "{}", String::from_utf8_lossy(self.0))
        }
    }
}

impl std::fmt::Display for NucleotidesViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.0.is_ascii_simd::<16>() {
            // SAFETY: we just checked it is ASCII using our fast SIMD function.
            // ASCII is valid UTF8.
            write!(f, "{}", unsafe { std::str::from_utf8_unchecked(self.0) })
        } else {
            write!(f, "{}", String::from_utf8_lossy(self.0))
        }
    }
}

impl std::fmt::Debug for Nucleotides {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl std::fmt::Debug for NucleotidesView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl std::fmt::Debug for NucleotidesViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}
