use crate::{
    data::{
        mappings::dna::*,
        types::nucleotides::NucleotidesMutable,
        validation::{recode::Recode, retain::RetainSequence},
        view_traits::SliceRange,
    },
    prelude::*,
};

use super::NucleotidesReadable;

pub trait ToDNA: Into<Nucleotides> {
    /// Filters and recodes to uppercase IUPAC with corrected gaps.
    fn filter_to_dna(self) -> Nucleotides {
        let mut n = self.into();
        n.retain_and_recode_dna(RefineDNAStrat::IupacCorrectGapsUc);
        n
    }

    /// Filters and recodes to uppercase IUPAC without gaps.
    fn filter_to_dna_uanligned(self) -> Nucleotides {
        let mut n = self.into();
        n.retain_and_recode_dna(RefineDNAStrat::IupacNoGapsUc);
        n
    }

    /// Recodes to uppercase IUPAC with corrected gaps in-place. Data that cannot
    /// be recoded becomes `N`.
    fn recode_to_dna(self) -> Nucleotides {
        let mut n = self.into();
        n.recode_dna_aligned();
        n
    }
}
impl ToDNA for String {}
impl ToDNA for Vec<u8> {}
impl ToDNA for &[u8] {}

/// Enumeration for DNA recoding strategies.
pub enum RecodeDNAStrat {
    /// Converts any valid IUPAC DNA to ACGTN (preserving case), allowing gaps
    IupacToAcgtnWithGaps,
    /// Converts any valid IUPAC DNA to uppercase ACGTN, allowing gaps
    IupacToAcgtnWithGapsUpper,

    /// Converts any byte to uppercase ACGTN with N as catch-all, even gaps
    AnyToAcgtnNoGapsUpper,
    /// Converts any byte to uppercase ACGTN with N as catch-all, allowing gaps
    AnyToAcgtnWithGapsUpper,

    /// Converts to IUPAC nomenclature without case changes, allowing gaps
    AnyToIupacWithGaps,
    /// Converts to uppercase IUPAC nomenclature, allowing gaps
    AnyToIupacWithGapsUpper,
    /// Converts to uppercase IUPAC nomenclature, correcting non-standard gaps
    AnyToIupacCorrectGapsUpper,
}

impl RecodeDNAStrat {
    /// Returns the corresponding mapping array for the selected recoding strategy
    #[inline]
    const fn mapping(&self) -> &'static [u8; 256] {
        match self {
            RecodeDNAStrat::IupacToAcgtnWithGaps => &IUPAC_TO_DNA_ACGTN_WITH_GAPS,
            RecodeDNAStrat::IupacToAcgtnWithGapsUpper => &IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC,

            RecodeDNAStrat::AnyToAcgtnNoGapsUpper => &ANY_TO_DNA_ACGTN_NO_GAPS_UC,
            RecodeDNAStrat::AnyToAcgtnWithGapsUpper => &ANY_TO_DNA_ACGTN_WITH_GAPS_UC,

            RecodeDNAStrat::AnyToIupacWithGaps => &ANY_TO_DNA_IUPAC_WITH_GAPS,
            RecodeDNAStrat::AnyToIupacWithGapsUpper => &ANY_TO_DNA_IUPAC_WITH_GAPS_UC,
            RecodeDNAStrat::AnyToIupacCorrectGapsUpper => &ANY_TO_DNA_IUPAC_CORRECT_GAPS_UC,
        }
    }
}

/// Enumeration for DNA validation and retention strategies.
pub enum IsValidDNA {
    /// For valid IUPAC bases without gaps
    IupacNoGaps,
    /// For uppercase IUPAC bases without gaps
    IupacNoGapsUc,
    /// For valid IUPAC bases with gaps
    IupacWithGaps,
    /// For uppercase IUPAC bases with gaps
    IupacWithGapsUc,

    /// For ACGTN bases without gaps
    AcgtnNoGaps,
    /// For uppercase ACGTN bases without gaps
    AcgtnNoGapsUc,
    /// For uppercase ACGTN bases with standard gaps
    AcgtnStdGapsUc,
}

impl IsValidDNA {
    #[inline]
    const fn validation_mapping(&self) -> &'static [bool; 256] {
        match self {
            IsValidDNA::IupacNoGaps => &IS_DNA_IUPAC_NO_GAPS,
            IsValidDNA::IupacNoGapsUc => &IS_DNA_IUPAC_NO_GAPS_UC,
            IsValidDNA::IupacWithGaps => &IS_DNA_IUPAC_WITH_GAPS,
            IsValidDNA::IupacWithGapsUc => &IS_DNA_IUPAC_WITH_GAPS_UC,

            IsValidDNA::AcgtnNoGaps => &IS_DNA_ACGTN_NO_GAPS,
            IsValidDNA::AcgtnNoGapsUc => &IS_DNA_ACGTN_NO_GAPS_UC,
            IsValidDNA::AcgtnStdGapsUc => &IS_DNA_ACGTN_STD_GAPS_UC,
        }
    }
}

/// DNA retention strategies. Data that cannot be recoded is not retained.
pub enum RefineDNAStrat {
    /// Retains and recodes to uppercase IUPAC bases without gaps
    IupacNoGapsUc,
    /// Retains and recodes to uppercase IUPAC bases with gaps
    IupacWithGapsUc,
    /// Retains and recodes to uppercase IUPAC bases with corrected gaps
    IupacCorrectGapsUc,

    /// Retains and recodes to uppercase ACGTN bases without gaps
    AcgtnNoGapsUc,
    /// Retains and recodes to uppercase ACGTN bases with gaps
    AcgtnWithGapsUc,
    /// Retains and recodes to uppercase ACGTN bases with standard gaps
    AcgtnStdGapsUc,
}

impl RefineDNAStrat {
    /// Returns the corresponding recoding array for the selected retention and recoding strategy
    #[inline]
    const fn recoding_mapping(&self) -> &'static [u8; 256] {
        match self {
            RefineDNAStrat::IupacNoGapsUc => &TO_DNA_IUPAC_NO_GAPS_UC,
            RefineDNAStrat::IupacWithGapsUc => &TO_DNA_IUPAC_WITH_GAPS_UC,
            RefineDNAStrat::IupacCorrectGapsUc => &TO_DNA_IUPAC_CORRECT_GAPS_UC,

            RefineDNAStrat::AcgtnNoGapsUc => &TO_DNA_ACGTN_NO_GAPS_UC,
            RefineDNAStrat::AcgtnWithGapsUc => &TO_DNA_ACGTN_WITH_GAPS_UC,
            RefineDNAStrat::AcgtnStdGapsUc => &TO_DNA_ACGTN_STD_GAPS_UC,
        }
    }
}

/// Provides DNA-specific methods for recoding a sequence. Data that cannot be
/// recoded becomes `N`. See [`RecodeDNAStrategy`] for recoding strategies.
pub trait RecodeNucleotides: NucleotidesMutable {
    /// Recodes the stored sequences according to the strategy in
    /// [`RecodeDNAStrategy`]. Data that cannot be recoded becomes `N`.
    #[inline]
    fn recode_dna(&mut self, strategy: RecodeDNAStrat) {
        self.nucleotide_mut_bytes().recode(strategy.mapping());
    }

    /// Recodes the stored sequence using
    /// [`RecodeDNAStrategy::AnyToAcgtnNoGapsUpper`], which is the preferred
    /// strategy for read data. Data that cannot be recoded becomes `N`.
    #[inline]
    fn recode_dna_reads(&mut self) {
        self.nucleotide_mut_bytes()
            .recode(RecodeDNAStrat::AnyToAcgtnNoGapsUpper.mapping());
    }

    /// Recodes the stored sequence using
    /// [`RecodeDNAStrategy::AnyToIupacCorrectGapsUpper`], which is the
    /// preferred strategy for aligned, multiple sequence alignment data. Data
    /// that cannot be recoded becomes `N`.
    #[inline]
    fn recode_dna_aligned(&mut self) {
        self.nucleotide_mut_bytes()
            .recode(RecodeDNAStrat::AnyToIupacCorrectGapsUpper.mapping());
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
    /// Retains nucleotides according if they are valid according to the
    /// retention strategy.
    #[inline]
    fn retain_dna(&mut self, strategy: IsValidDNA) {
        self.as_mut().retain_by_validation(strategy.validation_mapping());
    }

    /// Retains and recodes nucleotides according to the specified retention strategy.
    #[inline]
    fn retain_and_recode_dna(&mut self, strategy: RefineDNAStrat) {
        self.as_mut().retain_by_recoding(strategy.recoding_mapping());
    }
}

impl RetainNucleotides for Nucleotides {}
impl RetainNucleotides for Vec<u8> {}

pub trait CheckNucleotides: NucleotidesReadable {
    /// Checks if nucleotide sequence only contains valid IUPAC bases, without
    /// gaps
    #[inline]
    fn is_iupac_no_gaps(&self) -> bool {
        self.nucleotide_bytes().iter().all(|&b| IS_DNA_IUPAC_NO_GAPS[b as usize])
    }

    /// Checks if nucleotide sequence only contains valid uppercase IUPAC bases,
    /// without gaps
    #[inline]
    fn is_iupac_no_gaps_uc(&self) -> bool {
        self.nucleotide_bytes().iter().all(|&b| IS_DNA_IUPAC_NO_GAPS_UC[b as usize])
    }

    /// Checks if the nucleotide sequence only contains `A`, `C`, `G`, `T`, `N` (no gaps).
    #[inline]
    fn is_acgtn(&self) -> bool {
        self.nucleotide_bytes().iter().all(|&b| IS_DNA_ACGTN_NO_GAPS[b as usize])
    }

    /// Checks if the nucleotide sequence only contains uppercase `A`, `C`, `G`, `T`, `N` (no gaps).
    #[inline]
    fn is_acgtn_uc(&self) -> bool {
        self.nucleotide_bytes().iter().all(|&b| IS_DNA_ACGTN_NO_GAPS_UC[b as usize])
    }
}

impl<T: NucleotidesReadable> CheckNucleotides for T {}
