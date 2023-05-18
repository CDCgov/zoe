use crate::data::{err::DistanceError, vec_types::BiologicalSequence};

/// Amino Acid new type to help encourage type safety.
#[derive(Debug, Clone, Default)]
#[repr(transparent)]
pub struct AminoAcids(pub(crate) Vec<u8>);

impl AminoAcids {
    // Standard functions
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        AminoAcids(Vec::new())
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

    // Manipulation
    #[inline]
    pub fn find_and_replace(&mut self, needle: u8, replacement: u8) {
        crate::data::vec_types::find_and_replace(&mut self.0, needle, replacement);
    }

    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
    }

    /// # Distance
    ///
    /// Calculates hamming distance between `self`and another sequence.
    ///
    /// # Example
    /// ```
    /// use zoe::data::types::amino_acids::AminoAcids;
    ///
    /// let s1: AminoAcids = b"MANATEEMANATEEMANATEE".into();
    /// let s2: AminoAcids = b"MANGAEEMANATEEMANGAEE".into();
    ///
    /// assert!(4 == s1.distance_hamming(&s2));
    ///
    /// let s3: &[u8] = b"MANGAEEMANATEEMANGAEE";
    /// assert!(4 == s1.distance_hamming(&s3));
    /// ```
    ///
    #[inline]
    #[must_use]
    pub fn distance_hamming<T: BiologicalSequence + MaybeAmino>(&self, other_sequence: &T) -> usize {
        crate::distance::hamming_simd::<16>(&self.0, other_sequence.get_inner_ref())
    }

    /// Calculates physiochemical distance between `self`and another protein sequence. See: [`crate::distance::physiochemical`]
    ///
    /// # Citation
    /// For factor analysis used by the function:
    ///
    /// > Atchley et al. 2008. "Solving the protein sequence metric problem." Proc
    /// > Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Published 2005 Apr 25.
    ///
    ///
    /// # Errors
    /// If either argument is empty or the sequence characters are invalid, an error
    /// is thrown. See [`DistanceError`].
    ///
    /// # Example
    /// ```
    /// use zoe::data::types::amino_acids::AminoAcids;
    ///
    /// let s1: AminoAcids = b"MANATEEMANATEEMANATEE".into();
    /// let s2: AminoAcids = b"MANGAEEMANATEEMANGAEE".into();
    ///
    /// assert!( (s1.distance_physiochemical(&s2).unwrap() - 0.7643302).abs() < 0.01 );
    ///
    /// let s3: &[u8] = b"MANGAEEMANATEEMANGAEE";
    /// assert!( (s1.distance_physiochemical(&s3).unwrap() - 0.7643302).abs() < 0.01 );
    /// ```
    ///
    #[inline]
    pub fn distance_physiochemical<T: BiologicalSequence + MaybeAmino>(
        &self, other_sequence: &T,
    ) -> Result<f32, DistanceError> {
        crate::distance::physiochemical(&self.0, other_sequence.get_inner_ref())
    }
}

/// Marker trait to restrict usage to types that might possible contain nucleotides.
pub trait MaybeAmino {}
impl MaybeAmino for AminoAcids {}
impl MaybeAmino for Vec<u8> {}
impl MaybeAmino for &[u8] {}

impl FromIterator<u8> for AminoAcids {
    fn from_iter<T: IntoIterator<Item = u8>>(iter: T) -> Self {
        let mut v = Vec::new();
        for aa in iter {
            v.push(aa);
        }
        AminoAcids(v)
    }
}

impl From<Vec<u8>> for AminoAcids {
    fn from(vec: Vec<u8>) -> Self {
        AminoAcids(vec)
    }
}
impl From<&[u8]> for AminoAcids {
    fn from(bytes: &[u8]) -> Self {
        AminoAcids(bytes.to_vec())
    }
}

impl<const N: usize> From<&[u8; N]> for AminoAcids {
    fn from(bytes: &[u8; N]) -> Self {
        AminoAcids(bytes.to_vec())
    }
}
