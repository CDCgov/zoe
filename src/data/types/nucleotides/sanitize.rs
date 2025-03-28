use crate::{
    data::{
        mappings::{
            ANY_TO_DNA_ACGTN_UC, ANY_TO_DNA_IUPAC_WITH_GAPS, ANY_TO_DNA_IUPAC_WITH_GAPS_UC, IS_DNA_ACGTN, IS_DNA_ACGTN_UC,
            IS_DNA_IUPAC, IS_DNA_IUPAC_UC, IS_DNA_IUPAC_WITH_GAPS, IUPAC_TO_DNA_ACGTN, IUPAC_TO_DNA_ACGTN_UC,
            TO_DNA_IUPAC_UC, TO_DNA_IUPAC_WITH_GAPS_UC,
        },
        types::nucleotides::NucleotidesMutable,
        validation::{recode::Recode, retain::RetainSequence},
        view_traits::SliceRange,
    },
    prelude::*,
};

use super::NucleotidesReadable;

pub trait ToDNA: Into<Nucleotides> {
    /// Filters to unaligned, IUPAC DNA.
    fn filter_to_dna(self) -> Nucleotides {
        let mut n = self.into();
        n.retain_iupac_uc();
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
    fn recode_any_to_acgtn_uc(&mut self) {
        self.nucleotide_mut_bytes().recode(ANY_TO_DNA_ACGTN_UC);
    }
    /// Recodes the stored sequence to a seqeunce containing only valid IUPAC codes.
    /// Non-IUPAC bytes are replaced with N, case is left unchanged and gaps are
    /// left unchanged.
    #[inline]
    fn recode_any_to_iupac_with_gaps(&mut self) {
        self.nucleotide_mut_bytes().recode(ANY_TO_DNA_IUPAC_WITH_GAPS);
    }

    /// Recodes the stored sequence to a seqeunce containing only valid IUPAC codes.
    /// Non-IUPAC bytes are replaced with N, lowercase is changed to uppercase and gaps are
    /// left unchanged.
    #[inline]
    fn recode_any_to_iupac_with_gaps_uc(&mut self) {
        self.nucleotide_mut_bytes().recode(ANY_TO_DNA_IUPAC_WITH_GAPS_UC);
    }

    /// Recodes the stored sequence of valid IUPAC codes to a canonical (ACTG +
    /// N) sequence. Ambiguous bases become N while non-IUPAC bytes and gaps are
    /// left unchanged.
    #[inline]
    fn recode_iupac_to_acgtn(&mut self) {
        self.nucleotide_mut_bytes().recode(IUPAC_TO_DNA_ACGTN);
    }

    /// Recodes the stored sequence of valid IUPAC codes to an uppercase
    /// canonical (ACTG + N) sequence. Ambiguous bases become N while non-IUPAC
    /// bytes and gaps are left unchanged.
    #[inline]
    fn recode_iupac_to_acgtn_uc(&mut self) {
        self.nucleotide_mut_bytes().recode(IUPAC_TO_DNA_ACGTN_UC);
    }

    /// Masks the provided `range` with `N`. If the range does not exist, the
    /// function does nothing.
    #[inline]
    fn mask_if_exists<R: SliceRange>(&mut self, range: R) {
        self.nucleotide_mut_bytes().mask_if_exists(range, b'N');
    }
}

impl<T: NucleotidesMutable> RecodeNucleotides for T {}

pub trait RetainNucleotides: AsMut<Vec<u8>> {
    /// Only retains valid [IUPAC bases](https://www.bioinformatics.org/sms/iupac.html).
    /// Gaps are retained.
    #[inline]
    fn retain_iupac_with_gaps(&mut self) {
        self.as_mut().retain_by_validation(IS_DNA_IUPAC_WITH_GAPS);
    }

    /// Only retains valid [IUPAC](https://www.bioinformatics.org/sms/iupac.html) DNA and converts to uppercase.
    /// Gaps are retained.
    #[inline]
    fn retain_iupac_with_gaps_uc(&mut self) {
        self.as_mut().retain_by_recoding(TO_DNA_IUPAC_WITH_GAPS_UC);
    }

    /// Only retains valid [IUPAC bases](https://www.bioinformatics.org/sms/iupac.html).
    /// Gaps are NOT retained.
    #[inline]
    fn retain_iupac(&mut self) {
        self.as_mut().retain_by_validation(IS_DNA_IUPAC);
    }

    /// Only retains valid, unaligned [IUPAC](https://www.bioinformatics.org/sms/iupac.html) DNA and converts to uppercase.
    /// Gaps are NOT retained.
    #[inline]
    fn retain_iupac_uc(&mut self) {
        self.as_mut().retain_by_recoding(TO_DNA_IUPAC_UC);
    }
}

impl RetainNucleotides for Nucleotides {}
impl RetainNucleotides for Vec<u8> {}

pub trait CheckNucleotides: NucleotidesReadable {
    /// Checks if nucleotide sequence only contains valid IUPAC bases, without
    /// gaps
    #[inline]
    fn is_iupac(&self) -> bool {
        self.nucleotide_bytes().iter().all(|&b| IS_DNA_IUPAC[b as usize])
    }

    /// Checks if nucleotide sequence only contains valid uppercase IUPAC bases,
    /// without gaps
    #[inline]
    fn is_iupac_uc(&self) -> bool {
        self.nucleotide_bytes().iter().all(|&b| IS_DNA_IUPAC_UC[b as usize])
    }

    /// Checks if the nucleotide sequence only contains `A`, `C`, `G`, `T`, `N` (no gaps).
    #[inline]
    fn is_acgtn(&self) -> bool {
        self.nucleotide_bytes().iter().all(|&b| IS_DNA_ACGTN[b as usize])
    }

    /// Checks if the nucleotide sequence only contains uppercase `A`, `C`, `G`, `T`, `N` (no gaps).
    #[inline]
    fn is_acgtn_uc(&self) -> bool {
        self.nucleotide_bytes().iter().all(|&b| IS_DNA_ACGTN_UC[b as usize])
    }
}

impl<T: NucleotidesReadable> CheckNucleotides for T {}
