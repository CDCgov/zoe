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

    // Domain functions
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
