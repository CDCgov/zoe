use crate::data::err::DistanceError;

/// Amino Acid new type to help encourage type safety.
#[derive(Clone, Default, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct AminoAcids(pub(crate) Vec<u8>);

impl AminoAcids {
    // Standard functions
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        AminoAcids(Vec::new())
    }

    /// Consumes a [`Vec<u8>`] and return [`AminoAcids`] without checking for
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_vec_unchecked(v: Vec<u8>) -> Self {
        AminoAcids(v)
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
    pub fn into_vec(self) -> Vec<u8> {
        self.0
    }

    // Manipulation

    /// Truncates the length of the sequence to the specified `new_length`. This
    /// is equivalent to C-terminus trimming up to and including the index.
    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
    }

    /// Cuts the N-terminus of the [`AminoAcids`] just prior to the new starting
    /// index (0-based). Be aware that this clones the internal buffer!
    #[inline]
    pub fn cut_to_start(&mut self, new_start: usize) {
        *self = Self(self.0.drain(new_start..).collect());
    }

    /// If the end has a stop codon, remove it. Takes and gives ownership for
    /// chaining.
    #[must_use]
    #[inline]
    pub fn chop_stop(mut self) -> Self {
        if let Some(&b'*') = self.0.last() {
            self.0.pop();
            self
        } else {
            self
        }
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
    pub fn distance_hamming<T: AsRef<[u8]> + MaybeAmino>(&self, other_sequence: &T) -> usize {
        crate::distance::hamming_simd::<16>(&self.0, other_sequence.as_ref())
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
    pub fn distance_physiochemical<T: AsRef<[u8]> + MaybeAmino>(&self, other_sequence: &T) -> Result<f32, DistanceError> {
        crate::distance::physiochemical(&self.0, other_sequence.as_ref())
    }

    // Associated functions

    /// Generate a random AA sequence of given `length` and using a random `seed`.  Contains only uppercase, unaligned, non-ambiguous IUPAC codes.
    #[must_use]
    pub fn generate_random_aa(length: usize, seed: u64) -> Self {
        AminoAcids(crate::generate::rand_sequence(
            crate::data::constants::alphas::AMINO_ACIDS_UNALIGNED_UC,
            length,
            seed,
        ))
    }
}

/// Marker trait to restrict usage to types that might possible contain nucleotides.
pub trait MaybeAmino {}
impl MaybeAmino for AminoAcids {}
impl MaybeAmino for Vec<u8> {}
impl MaybeAmino for &[u8] {}

impl AsRef<[u8]> for AminoAcids {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for AminoAcids {
    #[inline]
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
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

impl FromIterator<u8> for AminoAcids {
    fn from_iter<T: IntoIterator<Item = u8>>(iter: T) -> Self {
        let mut v = Vec::new();
        for aa in iter {
            v.push(aa);
        }
        AminoAcids(v)
    }
}

impl std::fmt::Display for AminoAcids {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}

impl std::fmt::Debug for AminoAcids {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}

impl From<String> for AminoAcids {
    fn from(s: String) -> Self {
        AminoAcids(s.as_bytes().to_vec())
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
