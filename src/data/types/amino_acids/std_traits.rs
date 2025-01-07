use crate::{data::vec_types::CheckSequence, prelude::*};

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
        AminoAcids(s.as_bytes().to_vec())
    }
}

impl From<Vec<u8>> for AminoAcids {
    #[inline]
    fn from(vec: Vec<u8>) -> Self {
        AminoAcids(vec)
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

// TODO: Implement From &mut [u8] for AminoAcidsViewMut?? (and other seq-like types)

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

impl std::ops::Index<usize> for AminoAcids {
    type Output = u8;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl std::ops::Index<usize> for AminoAcidsView<'_> {
    type Output = u8;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl std::ops::Index<usize> for AminoAcidsViewMut<'_> {
    type Output = u8;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl std::ops::IndexMut<usize> for AminoAcids {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl std::ops::IndexMut<usize> for AminoAcidsViewMut<'_> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl std::ops::Index<std::ops::Range<usize>> for AminoAcids {
    type Output = [u8];

    #[inline]
    fn index(&self, index: std::ops::Range<usize>) -> &[u8] {
        &self.0[index]
    }
}

impl std::ops::Index<std::ops::Range<usize>> for AminoAcidsView<'_> {
    type Output = [u8];

    #[inline]
    fn index(&self, index: std::ops::Range<usize>) -> &[u8] {
        &self.0[index]
    }
}

impl std::ops::Index<std::ops::Range<usize>> for AminoAcidsViewMut<'_> {
    type Output = [u8];

    #[inline]
    fn index(&self, index: std::ops::Range<usize>) -> &[u8] {
        &self.0[index]
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
