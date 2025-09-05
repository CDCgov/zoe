use super::ByteIndexMap;
use crate::data::alphas::{AA_IUPAC_NO_GAPS_UC, AA_IUPAC_NO_GAPS_UC_X};

/// Conversion from IUPAC amino acids to indices, including `X` as a catch-all.
pub const AA_PROFILE_MAP: ByteIndexMap<21> = ByteIndexMap::new_ignoring_case(*AA_IUPAC_NO_GAPS_UC_X, b'X');

/// Conversion from unambiguous IUPAC amino acids to indices. Converting any
/// other symbol will result in Alanine.
pub const AA_UNAMBIG_PROFILE_MAP: ByteIndexMap<20> = ByteIndexMap::new_ignoring_case(*AA_IUPAC_NO_GAPS_UC, b'A');

/// Conversion from IUPAC amino acids, stop codons `*`, and ambiguous amino
/// acids `BJZ` to indices, including `X` as a catch-all.
///
/// See [IUPAC Standards](crate::data#iupac-standards) for more discussion.
pub const AA_ALL_AMBIG_PROFILE_MAP_WITH_STOP: ByteIndexMap<25> =
    ByteIndexMap::new_ignoring_case(*b"ACDEFGHIKLMNPQRSTVWY*BJZX", b'X');
