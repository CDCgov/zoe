use std::ops::Range;

/// Provides method for recoding a sequence based on a byte mapping.
pub trait Recode: AsMut<[u8]> {
    /// Recode the sequence data using a given byte mapping.
    #[inline]
    fn recode(&mut self, transformation_mapping: [u8; 256]) {
        for b in self.as_mut() {
            *b = transformation_mapping[*b as usize];
        }
    }

    /// Replaces the provided [`Range`] with the specified byte. If the range
    /// does not exist, the function does nothing.
    #[inline]
    fn mask_if_exists(&mut self, range: Range<usize>, replacement: u8) {
        if let Some(slice) = self.as_mut().get_mut(range) {
            for b in slice {
                *b = replacement;
            }
        }
    }
}

impl<T: AsMut<[u8]>> Recode for T {}
