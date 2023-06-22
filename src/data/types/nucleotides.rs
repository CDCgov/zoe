use crate::data::vec_types::BiologicalSequence;
use std::slice::ChunksExact;

#[derive(Debug, Clone, Default)]
#[repr(transparent)]
pub struct Nucleotides(pub(crate) Vec<u8>);

impl Nucleotides {
    // Standard functions
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
    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_slice()
    }

    #[inline]
    #[must_use]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        self.0.as_mut_slice()
    }

    #[inline]
    #[must_use]
    pub fn as_vec(&self) -> &Vec<u8> {
        &self.0
    }

    #[inline]
    #[must_use]
    pub fn as_mut_vec(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }

    #[inline]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    // Manipulation
    #[inline]
    pub fn find_and_replace(&mut self, needle: u8, replacement: u8) {
        crate::data::vec_types::find_and_replace(&mut self.0, needle, replacement);
    }

    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
    }

    // Domain functions
    #[inline]
    #[must_use]
    pub fn reverse_complement(&self) -> Self {
        Self(reverse_complement(&self.0))
    }

    #[inline]
    #[must_use]
    pub fn translate(&self) -> Vec<u8> {
        translate_sequence(&self.0)
    }

    #[must_use]
    pub fn into_aa_iter(&self) -> TranslatedNucleotidesIter {
        TranslatedNucleotidesIter {
            codons:        self.0.chunks_exact(3),
            has_remainder: true,
        }
    }

    /// # Distance
    ///
    /// Calculates hamming distance between `self`and another sequence.
    ///
    /// # Example
    /// ```
    /// use zoe::data::types::nucleotides::Nucleotides;
    ///
    /// let s1: Nucleotides = b"ATGCATCGATCGATCGATCGATCGATCGATGC".into();
    /// let s2: Nucleotides = b"ATGCATnGATCGATCGATCGAnCGATCGATnC".into();
    ///
    /// assert!(3 == s1.distance_hamming(&s2));
    ///
    /// let s3: &[u8] = b"ATGCATnGATCGATCGATCGAnCGATCGATnC";
    /// assert!(3 == s1.distance_hamming(&s3));
    /// ```
    ///
    #[inline]
    #[must_use]
    pub fn distance_hamming<T: BiologicalSequence + MaybeNucleic>(&self, other_sequence: &T) -> usize {
        crate::distance::hamming_simd::<16>(&self.0, other_sequence.get_inner_ref())
    }
}

/// Marker trait to restrict usage to types that might possible contain nucleotides.
pub trait MaybeNucleic {}
impl MaybeNucleic for Nucleotides {}
impl MaybeNucleic for Vec<u8> {}
impl MaybeNucleic for &[u8] {}

#[inline]
#[must_use]
pub fn translate_sequence(s: &[u8]) -> Vec<u8> {
    use crate::data::matrices::GENETIC_CODE;

    let mut codons = s.chunks_exact(3);
    let mut aa_sequence = Vec::with_capacity(s.len() / 3 + 1);

    for codon in codons.by_ref() {
        aa_sequence.push(if is_partial_codon(codon) {
            b'~'
        } else {
            *GENETIC_CODE.get(&codon.to_ascii_uppercase()).unwrap_or(&b'X')
        });
    }

    let tail = codons.remainder();
    if is_partial_codon(tail) {
        aa_sequence.push(b'~');
    }

    aa_sequence
}

pub struct TranslatedNucleotidesIter<'a> {
    codons:        ChunksExact<'a, u8>,
    has_remainder: bool,
}

impl<'a> Iterator for TranslatedNucleotidesIter<'a> {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        use crate::data::matrices::GENETIC_CODE;

        if let Some(codon) = self.codons.next() {
            if is_partial_codon(codon) {
                Some(b'~')
            } else {
                Some(*GENETIC_CODE.get(&codon.to_ascii_uppercase()).unwrap_or(&b'X'))
            }
        } else if self.has_remainder && is_partial_codon(self.codons.remainder()) {
            self.has_remainder = false;
            Some(b'~')
        } else {
            None
        }
    }
}

#[inline]
#[must_use]
fn is_partial_codon(codon: &[u8]) -> bool {
    if codon.is_empty() {
        false
    } else if codon.len() < 3 {
        true
    } else {
        let count: u8 = codon.iter().map(|b| u8::from(*b == b'-' || *b == b'~' || *b == b'.')).sum();

        count == 1 || count == 2
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

impl std::ops::Index<std::ops::Range<usize>> for Nucleotides {
    type Output = [u8];

    #[inline]
    fn index(&self, index: std::ops::Range<usize>) -> &[u8] {
        &self.0[index]
    }
}

#[inline]
#[must_use]
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    use crate::data::matrices::REV_COMP;
    bases.iter().rev().copied().map(|x| REV_COMP[x as usize]).collect()
}

impl std::fmt::Display for Nucleotides {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_translate() {
        let s = Nucleotides(b"ATGTCAGATcccagagaaTGAgg".to_vec());
        assert_eq!(s.translate(), b"MSDPRE*~");
    }

    #[test]
    fn test_translate_collect() {
        use super::super::amino_acids::AminoAcids;
        let s = Nucleotides(b"ATGTCAGATcccagagaaTGAgg".to_vec());
        assert_eq!(s.into_aa_iter().collect::<AminoAcids>().0, b"MSDPRE*~".to_vec());
    }
}
