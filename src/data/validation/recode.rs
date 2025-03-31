use crate::data::view_traits::SliceRange;

/// Provides method for recoding a sequence based on a byte mapping.
pub trait Recode: AsMut<[u8]> {
    /// Recode the sequence data using a given byte mapping.
    #[inline]
    fn recode(&mut self, transformation_mapping: &'static [u8; 256]) {
        for b in self.as_mut() {
            *b = transformation_mapping[*b as usize];
        }
    }

    /// Replaces the provided `range` with the specified byte. If the range does
    /// not exist, the function does nothing.
    #[inline]
    fn mask_if_exists<R: SliceRange>(&mut self, range: R, replacement: u8) {
        if let Some(slice) = self.as_mut().get_mut(range) {
            for b in slice {
                *b = replacement;
            }
        }
    }
}

impl<T: AsMut<[u8]>> Recode for T {}
