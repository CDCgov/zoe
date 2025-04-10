use crate::{
    DEFAULT_SIMD_LANES,
    data::{
        mappings::dna::*,
        types::nucleotides::{NucleotidesMutable, NucleotidesReadable},
        validation::{recode::Recode, retain::RetainSequence},
        view_traits::SliceRange,
    },
    prelude::*,
};
use std::simd::prelude::*;

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

/// Enumeration for DNA recoding strategies. In all strategies, U is recoded to
/// T.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[non_exhaustive]
pub enum RecodeDNAStrat {
    /// Converts any valid IUPAC DNA to ACGTN, preserving case and gaps
    IupacToAcgtnWithGaps,
    /// Converts any valid IUPAC DNA to uppercase ACGTN, preserving gaps
    IupacToAcgtnWithGapsUpper,

    /// Converts any byte to uppercase ACGTN with N as catch-all, even gaps
    AnyToAcgtnNoGapsUpper,
    /// Converts any byte to uppercase ACGTN with N as catch-all, preserving
    /// gaps
    AnyToAcgtnWithGapsUpper,

    /// Converts to IUPAC nomenclature, preserving case and gaps. Capital N is
    /// used for non-IUPAC symbols
    AnyToIupacWithGaps,
    /// Converts to uppercase IUPAC nomenclature, preserving gaps. Capital N is
    /// used for non-IUPAC symbols
    AnyToIupacWithGapsUpper,
    /// Converts to uppercase IUPAC nomenclature, correcting non-standard gaps
    /// (`:` and `~`) to `-`. Gaps represented by `.` are preserved. Capital N
    /// is used for non-IUPAC symbols
    AnyToIupacCorrectGapsUpper,
}

impl RecodeDNAStrat {
    /// Returns the corresponding mapping array for the selected recoding strategy
    #[inline]
    const fn mapping(self) -> &'static [u8; 256] {
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
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[non_exhaustive]
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
    pub(crate) const fn mapping(self) -> &'static [bool; 256] {
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

    #[inline]
    pub(crate) const fn is_valid(self, index: u8) -> bool {
        self.mapping()[index as usize]
    }
}

/// DNA retention strategies. Data that cannot be recoded is not retained.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[non_exhaustive]
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
    /// Retains and recodes to uppercase ACGTN bases with standard gaps (`-`)
    AcgtnStdGapsUc,
}

impl RefineDNAStrat {
    /// Returns the corresponding recoding array for the selected retention and recoding strategy
    #[inline]
    const fn mapping(self) -> &'static [u8; 256] {
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
/// recoded becomes `N`. See [`RecodeDNAStrat`] for recoding strategies.
pub trait RecodeNucleotides: NucleotidesMutable {
    /// Recodes the stored sequences according to the strategy in
    /// [`RecodeDNAStrat`]. Data that cannot be recoded becomes `N`.
    #[inline]
    fn recode_dna(&mut self, strategy: RecodeDNAStrat) {
        self.nucleotide_mut_bytes().recode(strategy.mapping());
    }

    /// Recodes the stored sequence using
    /// [`RecodeDNAStrat::AnyToAcgtnNoGapsUpper`], which is the preferred
    /// strategy for read data. Data that cannot be recoded becomes `N`.
    #[inline]
    fn recode_dna_reads(&mut self) {
        self.nucleotide_mut_bytes()
            .recode(RecodeDNAStrat::AnyToAcgtnNoGapsUpper.mapping());
    }

    /// Recodes the stored sequence using
    /// [`RecodeDNAStrat::AnyToIupacCorrectGapsUpper`], which is the preferred
    /// strategy for aligned, multiple sequence alignment data. Data that cannot
    /// be recoded becomes `N`.
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
        self.as_mut().retain_by_validation(strategy.mapping());
    }

    /// Retains and recodes nucleotides according to the specified retention strategy.
    #[inline]
    fn retain_and_recode_dna(&mut self, strategy: RefineDNAStrat) {
        self.as_mut().retain_by_recoding(strategy.mapping());
    }
}

impl RetainNucleotides for Nucleotides {}
impl RetainNucleotides for Vec<u8> {}

pub trait CheckNucleotides {
    /// Checks if nucleotide sequence is valid according to the specified
    /// validation strategy.
    ///
    /// If calling with `strategy` as [`AcgtnNoGapsUc`], consider using
    /// [`is_acgtn_uc`] instead, which is SIMD accelerated.
    ///
    /// [`AcgtnNoGapsUc`]: IsValidDNA::AcgtnNoGapsUc
    /// [`is_acgtn_uc`]: CheckNucleotides::is_acgtn_uc
    fn is_valid_dna(&self, strategy: IsValidDNA) -> bool;

    /// Checks if the nucleotide sequence only contains uppercase `A`, `C`, `G`,
    /// `T`, `N` (no gaps). This version is SIMD accelerated.
    ///
    /// The output will be identical to [`is_valid_dna`] with `strategy` as
    /// [`AcgtnNoGapsUc`].
    ///
    /// [`AcgtnNoGapsUc`]: IsValidDNA::AcgtnNoGapsUc
    /// [`is_valid_dna`]: CheckNucleotides::is_valid_dna
    fn is_acgtn_uc(&self) -> bool;
}

impl<T: NucleotidesReadable> CheckNucleotides for T {
    #[inline]
    fn is_valid_dna(&self, strategy: IsValidDNA) -> bool {
        self.nucleotide_bytes().iter().all(|&b| strategy.is_valid(b))
    }

    #[inline]
    fn is_acgtn_uc(&self) -> bool {
        is_acgtn_uc_simd(self.nucleotide_bytes())
    }
}

impl CheckNucleotides for &[u8] {
    #[inline]
    fn is_valid_dna(&self, strategy: IsValidDNA) -> bool {
        self.iter().all(|&b| strategy.is_valid(b))
    }

    #[inline]
    fn is_acgtn_uc(&self) -> bool {
        is_acgtn_uc_simd(self)
    }
}

impl CheckNucleotides for Vec<u8> {
    #[inline]
    fn is_valid_dna(&self, strategy: IsValidDNA) -> bool {
        self.iter().all(|&b| strategy.is_valid(b))
    }

    #[inline]
    fn is_acgtn_uc(&self) -> bool {
        is_acgtn_uc_simd(self)
    }
}

#[inline]
#[allow(non_snake_case)]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
fn is_acgtn_uc_simd(s: &[u8]) -> bool {
    let (left, middle, right) = s.as_simd::<{ DEFAULT_SIMD_LANES }>();

    let (A, G, C, T, N) = (
        Simd::splat(b'A'),
        Simd::splat(b'G'),
        Simd::splat(b'C'),
        Simd::splat(b'T'),
        Simd::splat(b'N'),
    );

    for v in middle {
        let valid = v.simd_eq(A) | v.simd_eq(G) | v.simd_eq(C) | v.simd_eq(T) | v.simd_eq(N);
        if !valid.all() {
            return false;
        }
    }

    left.iter().chain(right).all(|&b| IsValidDNA::AcgtnNoGapsUc.is_valid(b))
}
