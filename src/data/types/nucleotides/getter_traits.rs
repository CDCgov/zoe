use super::{FastQ, FastQView, FastQViewMut, Nucleotides, NucleotidesView, NucleotidesViewMut};

/// Getter trait for structures providing read access to nucleotides
pub trait NucleotidesReadable {
    /// Get the nucleotides from the struct as a byte slice
    fn nucleotide_bytes(&self) -> &[u8];
}

/// Getter trait for structures providing mutable access to nucleotides
pub trait NucleotidesMutable: NucleotidesReadable {
    /// Get the nucleotides from the struct as a mutable byte slice
    fn nucleotide_mut_bytes(&mut self) -> &mut [u8];
}

impl NucleotidesReadable for Nucleotides {
    #[inline]
    fn nucleotide_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl NucleotidesReadable for NucleotidesView<'_> {
    #[inline]
    fn nucleotide_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl NucleotidesReadable for NucleotidesViewMut<'_> {
    #[inline]
    fn nucleotide_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl NucleotidesMutable for Nucleotides {
    #[inline]
    fn nucleotide_mut_bytes(&mut self) -> &mut [u8] {
        self.as_mut()
    }
}

impl NucleotidesMutable for NucleotidesViewMut<'_> {
    #[inline]
    fn nucleotide_mut_bytes(&mut self) -> &mut [u8] {
        self.as_mut()
    }
}

impl NucleotidesReadable for FastQ {
    #[inline]
    fn nucleotide_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl NucleotidesReadable for FastQView<'_> {
    #[inline]
    fn nucleotide_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl NucleotidesReadable for FastQViewMut<'_> {
    #[inline]
    fn nucleotide_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl NucleotidesMutable for FastQ {
    #[inline]
    fn nucleotide_mut_bytes(&mut self) -> &mut [u8] {
        self.sequence.as_mut()
    }
}

impl NucleotidesMutable for FastQViewMut<'_> {
    #[inline]
    fn nucleotide_mut_bytes(&mut self) -> &mut [u8] {
        self.sequence.as_mut()
    }
}
