use super::{AminoAcids, AminoAcidsView, AminoAcidsViewMut};

/// Getter trait for structures providing read access to amino acids
pub trait AminoAcidsReadable {
    /// Get the amino acids from the struct as a byte slice
    fn amino_acids_bytes(&self) -> &[u8];
}

/// Getter trait for structures providing mutable access to amino acids
pub trait AminoAcidsMutable: AminoAcidsReadable {
    /// Get the amino acids from the struct as a mutable byte slice
    fn amino_acids_mut_bytes(&mut self) -> &mut [u8];
}

impl AminoAcidsReadable for AminoAcids {
    #[inline]
    fn amino_acids_bytes(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl AminoAcidsReadable for AminoAcidsView<'_> {
    #[inline]
    fn amino_acids_bytes(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl AminoAcidsReadable for AminoAcidsViewMut<'_> {
    #[inline]
    fn amino_acids_bytes(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl AminoAcidsMutable for AminoAcids {
    #[inline]
    fn amino_acids_mut_bytes(&mut self) -> &mut [u8] {
        self.as_mut_bytes()
    }
}

impl AminoAcidsMutable for AminoAcidsViewMut<'_> {
    #[inline]
    fn amino_acids_mut_bytes(&mut self) -> &mut [u8] {
        self.as_mut_bytes()
    }
}
