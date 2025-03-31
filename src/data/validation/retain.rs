use crate::data::types::{amino_acids::AminoAcids, nucleotides::Nucleotides};

/// Provides methods for validating and transforming sequence data using byte
/// mappings
pub trait RetainSequence {
    /// Retains only bytes that are marked as valid in the validation mapping
    fn retain_by_validation(&mut self, validation_mapping: &'static [bool; 256]);
    /// Retains and transforms bytes using the transformation mapping, removing
    /// any that map to 0
    fn retain_by_recoding(&mut self, transformation_mapping: &'static [u8; 256]);
}

impl RetainSequence for Vec<u8> {
    /// Allows for filtering of biological sequences using a validation mapping.
    #[inline]
    fn retain_by_validation(&mut self, validation_mapping: &'static [bool; 256]) {
        self.retain(|b| validation_mapping[*b as usize]);
    }

    /// Allows for filtering of biological sequences using byte re-encoding. The
    /// 0-byte is assumed to be an invalid pattern.
    #[inline]
    fn retain_by_recoding(&mut self, transformation_mapping: &'static [u8; 256]) {
        self.retain_mut(|b| {
            *b = transformation_mapping[*b as usize];
            *b > 0
        });
    }
}

/// Provides pass through functionality for the standard library.
pub trait StdForSequences: AsMut<Vec<u8>> {
    #[inline]
    fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&u8) -> bool, {
        self.as_mut().retain(f);
    }
}
impl StdForSequences for Nucleotides {}
impl StdForSequences for AminoAcids {}
