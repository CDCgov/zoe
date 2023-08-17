/// Utility trait for sequence validation using byte-wise mappings.
pub trait ValidateSequence {
    fn retain_by_validation(&mut self, validation_mapping: [bool; 256]);
    fn retain_by_recoding(&mut self, transformation_mapping: [u8; 256]);
}

impl ValidateSequence for Vec<u8> {
    /// Allows for filtering of biological sequences using a validation mapping.
    #[inline]
    fn retain_by_validation(&mut self, validation_mapping: [bool; 256]) {
        self.retain(|b| validation_mapping[*b as usize]);
    }

    /// Allows for filtering of biological sequences using byte re-encoding. The
    /// 0-byte is assumed to be an invalid pattern.
    #[inline]
    fn retain_by_recoding(&mut self, validation_mapping: [u8; 256]) {
        self.retain_mut(|b| {
            *b = validation_mapping[*b as usize];
            *b > 0
        });
    }
}
