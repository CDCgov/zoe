use super::*;

// Conversion
impl AsRef<[u8]> for Nucleotides {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for Nucleotides {
    #[inline]
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
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
    fn from(s: String) -> Self {
        Nucleotides(s.as_bytes().to_vec())
    }
}

impl From<Vec<u8>> for Nucleotides {
    fn from(vec: Vec<u8>) -> Self {
        Nucleotides(vec)
    }
}

impl From<&[u8]> for Nucleotides {
    fn from(bytes: &[u8]) -> Self {
        Nucleotides(bytes.to_vec())
    }
}

impl<const N: usize> From<&[u8; N]> for Nucleotides {
    fn from(bytes: &[u8; N]) -> Self {
        Nucleotides(bytes.to_vec())
    }
}

impl IntoIterator for Nucleotides {
    type Item = u8;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a Nucleotides {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

#[allow(dead_code)]
impl Nucleotides {
    fn iter(&self) -> std::slice::Iter<'_, u8> {
        <&Self as IntoIterator>::into_iter(self)
    }
}

impl FromIterator<u8> for Nucleotides {
    fn from_iter<T: IntoIterator<Item = u8>>(iterable: T) -> Self {
        Nucleotides(iterable.into_iter().collect())
    }
}

// Indexing
impl std::ops::Index<usize> for Nucleotides {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl std::ops::IndexMut<usize> for Nucleotides {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl std::ops::Index<std::ops::Range<usize>> for Nucleotides {
    type Output = [u8];

    #[inline]
    fn index(&self, index: std::ops::Range<usize>) -> &[u8] {
        &self.0[index]
    }
}

// Display
impl std::fmt::Display for Nucleotides {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}
