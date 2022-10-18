#[derive(Debug, Clone, Default)]
pub struct Nucleotides(pub(crate) Vec<u8>);

impl Nucleotides {
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Nucleotides(Vec::new())
    }

    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[inline]
    #[must_use]
    pub fn reverse_complement(&self) -> Self {
        Self(reverse_complement(&self.0))
    }
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_slice()
    }

    #[inline]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        self.0.as_mut_slice()
    }

    #[inline]
    #[must_use]
    pub fn as_vec(&self) -> &Vec<u8> {
        &self.0
    }

    #[inline]
    pub fn as_mut_vec(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }

    #[inline]
    pub fn find_and_replace(&mut self, needle: u8, replacement: u8) {
        crate::data::vec_types::find_and_replace(&mut self.0, needle, replacement);
    }

    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
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

impl FromIterator<u8> for Nucleotides {
    fn from_iter<T: IntoIterator<Item = u8>>(iterable: T) -> Self {
        Nucleotides(iterable.into_iter().collect())
    }
}

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

#[inline]
#[must_use]
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    use crate::data::matrices::REV_COMP;
    bases
        .iter()
        .rev()
        .copied()
        .map(|x| REV_COMP[x as usize])
        .collect()
}

impl std::fmt::Display for Nucleotides {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}
