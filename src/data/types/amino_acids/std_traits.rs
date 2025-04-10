use crate::{data::validation::CheckSequence, prelude::*};
use std::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

impl AsRef<[u8]> for AminoAcids {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsRef<[u8]> for AminoAcidsView<'_> {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl AsRef<[u8]> for AminoAcidsViewMut<'_> {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl AsMut<[u8]> for AminoAcids {
    #[inline]
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl AsMut<[u8]> for AminoAcidsViewMut<'_> {
    #[inline]
    fn as_mut(&mut self) -> &mut [u8] {
        self.0
    }
}

impl AsRef<Vec<u8>> for AminoAcids {
    #[inline]
    fn as_ref(&self) -> &Vec<u8> {
        &self.0
    }
}

impl AsMut<Vec<u8>> for AminoAcids {
    #[inline]
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }
}

impl From<String> for AminoAcids {
    #[inline]
    fn from(s: String) -> Self {
        AminoAcids(s.into_bytes())
    }
}

impl From<Vec<u8>> for AminoAcids {
    #[inline]
    fn from(vec: Vec<u8>) -> Self {
        AminoAcids(vec)
    }
}

impl<const N: usize> From<[u8; N]> for AminoAcids {
    #[inline]
    fn from(bytes: [u8; N]) -> Self {
        AminoAcids(bytes.to_vec())
    }
}

impl From<&[u8]> for AminoAcids {
    #[inline]
    fn from(bytes: &[u8]) -> Self {
        AminoAcids(bytes.to_vec())
    }
}

impl<'a> From<&'a [u8]> for AminoAcidsView<'a> {
    #[inline]
    fn from(bytes: &'a [u8]) -> Self {
        AminoAcidsView(bytes)
    }
}

impl From<&mut [u8]> for AminoAcids {
    #[inline]
    fn from(bytes: &mut [u8]) -> Self {
        AminoAcids(bytes.to_vec())
    }
}

impl<'a> From<&'a mut [u8]> for AminoAcidsView<'a> {
    #[inline]
    fn from(bytes: &'a mut [u8]) -> Self {
        AminoAcidsView(bytes)
    }
}

impl<'a> From<&'a mut [u8]> for AminoAcidsViewMut<'a> {
    #[inline]
    fn from(bytes: &'a mut [u8]) -> Self {
        AminoAcidsViewMut(bytes)
    }
}

impl<const N: usize> From<&[u8; N]> for AminoAcids {
    #[inline]
    fn from(bytes: &[u8; N]) -> Self {
        AminoAcids(bytes.to_vec())
    }
}

impl<'a, const N: usize> From<&'a [u8; N]> for AminoAcidsView<'a> {
    #[inline]
    fn from(bytes: &'a [u8; N]) -> Self {
        AminoAcidsView(bytes)
    }
}

impl<const N: usize> From<&mut [u8; N]> for AminoAcids {
    #[inline]
    fn from(bytes: &mut [u8; N]) -> Self {
        AminoAcids(bytes.to_vec())
    }
}

impl<'a, const N: usize> From<&'a mut [u8; N]> for AminoAcidsView<'a> {
    #[inline]
    fn from(bytes: &'a mut [u8; N]) -> Self {
        AminoAcidsView(bytes)
    }
}

impl<'a, const N: usize> From<&'a mut [u8; N]> for AminoAcidsViewMut<'a> {
    #[inline]
    fn from(bytes: &'a mut [u8; N]) -> Self {
        AminoAcidsViewMut(bytes)
    }
}

impl IntoIterator for AminoAcids {
    type Item = u8;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for AminoAcidsView<'a> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for AminoAcidsViewMut<'a> {
    type Item = &'a mut u8;
    type IntoIter = std::slice::IterMut<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<'a> IntoIterator for &'a AminoAcids {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a AminoAcidsView<'_> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a AminoAcidsViewMut<'_> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a mut AminoAcids {
    type Item = &'a mut u8;
    type IntoIter = std::slice::IterMut<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<'a> IntoIterator for &'a mut AminoAcidsViewMut<'_> {
    type Item = &'a mut u8;
    type IntoIter = std::slice::IterMut<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl FromIterator<u8> for AminoAcids {
    #[inline]
    fn from_iter<T: IntoIterator<Item = u8>>(iterable: T) -> Self {
        AminoAcids(iterable.into_iter().collect())
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for AminoAcids {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for AminoAcidsView<'_> {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for AminoAcidsViewMut<'_> {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> IndexMut<I> for AminoAcids {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> IndexMut<I> for AminoAcidsViewMut<'_> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl std::fmt::Display for AminoAcids {
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

impl std::fmt::Display for AminoAcidsView<'_> {
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

impl std::fmt::Display for AminoAcidsViewMut<'_> {
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

impl std::fmt::Debug for AminoAcids {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl std::fmt::Debug for AminoAcidsView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl std::fmt::Debug for AminoAcidsViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}
