use crate::{
    data::{
        mappings::{
            ANY_TO_DNA_CANONICAL_UPPER, IS_IUPAC_BASE, IS_UNALIGNED_IUPAC_BASE, IUPAC_TO_DNA_CANONICAL,
            IUPAC_TO_DNA_CANONICAL_UPPER, TO_DNA_UC, TO_UNALIGNED_DNA_UC,
        },
        types::nucleotides::NucleotidesMutable,
        validation::recode::Recode,
        validation::retain::RetainSequence,
    },
    prelude::*,
};
use std::ops::Range;

pub trait ToDNA: AsRef<[u8]> {
    fn filter_to_dna(&self) -> Nucleotides {
        let mut n = Nucleotides(self.as_ref().to_vec());
        n.retain_dna_uc();
        n
    }
}
impl ToDNA for String {}
impl ToDNA for Vec<u8> {}
impl ToDNA for &[u8] {}

/// Provides DNA-specific methods for recoding a sequence.
pub trait RecodeNucleotides: NucleotidesMutable {
    /// Recodes the stored sequence to an uppercase canonical (ACTG + N) one.
    /// Any non-canonical base becomes N, including gaps.
    #[inline]
    fn recode_any_to_actgn_uc(&mut self) {
        self.nucleotide_mut_bytes().recode(ANY_TO_DNA_CANONICAL_UPPER);
    }

    /// Recodes the stored sequence of valid IUPAC codes to a canonical (ACTG +
    /// N) sequence. Ambiguous bases become N while non-IUPAC bytes and gaps are
    /// left unchanged.
    #[inline]
    fn recode_iupac_to_actgn(&mut self) {
        self.nucleotide_mut_bytes().recode(IUPAC_TO_DNA_CANONICAL);
    }

    /// Recodes the stored sequence of valid IUPAC codes to an uppercase
    /// canonical (ACTG + N) sequence. Ambiguous bases become N while non-IUPAC
    /// bytes and gaps are left unchanged.
    #[inline]
    fn recode_iupac_to_actgn_uc(&mut self) {
        self.nucleotide_mut_bytes().recode(IUPAC_TO_DNA_CANONICAL_UPPER);
    }

    /// Masks the provided [`Range`] with `N`. If the range does not exist, the
    /// function does nothing.
    fn mask_if_exists(&mut self, range: Range<usize>) {
        self.nucleotide_mut_bytes().mask_if_exists(range, b'N');
    }
}
impl<T: NucleotidesMutable> RecodeNucleotides for T {}

pub trait RetainNucleotides: AsMut<Vec<u8>> {
    /// Only retains valid [IUPAC bases](https://www.bioinformatics.org/sms/iupac.html).
    /// Gaps are retained.
    #[inline]
    fn retain_iupac_bases(&mut self) {
        self.as_mut().retain_by_validation(IS_IUPAC_BASE);
    }

    /// Only retains valid [IUPAC bases](https://www.bioinformatics.org/sms/iupac.html).
    /// Gaps are NOT retained.
    #[inline]
    fn retain_unaligned_bases(&mut self) {
        self.as_mut().retain_by_validation(IS_UNALIGNED_IUPAC_BASE);
    }

    /// Only retains valid [IUPAC](https://www.bioinformatics.org/sms/iupac.html) DNA and converts to uppercase.
    /// Gaps are retained.
    #[inline]
    fn retain_dna_uc(&mut self) {
        self.as_mut().retain_by_recoding(TO_DNA_UC);
    }

    /// Only retains valid, unaligned [IUPAC](https://www.bioinformatics.org/sms/iupac.html) DNA and converts to uppercase.
    /// Gaps are NOT retained.
    #[inline]
    fn retain_unaligned_dna_uc(&mut self) {
        self.as_mut().retain_by_recoding(TO_UNALIGNED_DNA_UC);
    }
}

impl RetainNucleotides for Nucleotides {}
impl RetainNucleotides for Vec<u8> {}
