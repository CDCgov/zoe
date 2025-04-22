use super::ByteIndexMap;
use crate::data::alphas::AA_IUPAC_NO_GAPS_UC;

/// Conversion from IUPAC amino acids to indices. Alanine is used as a
/// catch-all.
pub const AA_UNAMBIG_PROFILE_MAP: ByteIndexMap<20> = ByteIndexMap::new_ignoring_case(*AA_IUPAC_NO_GAPS_UC, b'A');
