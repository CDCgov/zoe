use super::*;
use crate::data::alphas::*;

//
// For if: u8 -> bool
//

/// A boolean mapping of valid IUPAC nucleotide codes without gaps.
pub(crate) const IS_DNA_IUPAC_NO_GAPS: [bool; 256] = make_is_alpha_mapping(DNA_IUPAC_NO_GAPS);

/// A boolean mapping of uppercase IUPAC nucleotide codes without gaps.
pub(crate) const IS_DNA_IUPAC_NO_GAPS_UC: [bool; 256] = make_is_alpha_mapping(DNA_IUPAC_NO_GAPS_UC);

/// A boolean mapping of all valid IUPAC nucleotide codes inluding gaps.
pub(crate) const IS_DNA_IUPAC_WITH_GAPS: [bool; 256] = make_is_alpha_mapping(DNA_IUPAC_WITH_GAPS);

/// A boolean mapping of uppercase IUPAC nucleotide codes including.
pub(crate) const IS_DNA_IUPAC_WITH_GAPS_UC: [bool; 256] = make_is_alpha_mapping(DNA_IUPAC_WITH_GAPS_UC);

/// A boolean mapping of strictly canonical nucleotides + `n`/`N`.
pub(crate) const IS_DNA_ACGTN_NO_GAPS: [bool; 256] = make_is_alpha_mapping(DNA_ACGTN_NO_GAPS);

/// A boolean mapping of uppercase `ACGTN` bytes.
pub(crate) const IS_DNA_ACGTN_NO_GAPS_UC: [bool; 256] = make_is_alpha_mapping(DNA_ACGTN_NO_GAPS_UC);

/// A boolean mapping of uppercase `ACGTN-` bytes.
pub(crate) const IS_DNA_ACGTN_STD_GAPS_UC: [bool; 256] = make_is_alpha_mapping(DNA_ACGTN_STD_GAPS_UC);

//
// For retain: u8 -> u8 else 0
//

/// Used to convert nucleotide sequences to unaligned, IUPAC-validated,
/// uppercase DNA. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_IUPAC_NO_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgturyswkmbdhvnACGTURYSWKMBDHVN",
    b"ACGTTRYSWKMBDHVNACGTTRYSWKMBDHVN",
    0,
);

/// Used to convert nucleotide sequences to valid IUPAC DNA (uppercase). The
/// 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_IUPAC_WITH_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-",
    b"ACGTTRYSWKMBDHVNACGTTRYSWKMBDHVN.-",
    0,
);

/// Used to convert nucleotide sequences to valid IUPAC DNA (uppercase), also
/// correct non-standard gaps. The 0-byte is used for filtering out unwanted
/// patterns.
pub(crate) const TO_DNA_IUPAC_CORRECT_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-:~",
    b"ACGTTRYSWKMBDHVNACGTTRYSWKMBDHVN.---",
    0,
);

/// Used to convert nucleotide sequences to only ACGTN (uppercase) without gaps.
/// The 0-byte is used for filtering out unwanted patterns.
#[rustfmt::skip]
pub(crate) const TO_DNA_ACGTN_NO_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgtnACGTN",
    b"ACGTNACGTN",
    0,
);

/// Used to convert nucleotide sequences to only uppercase ACGTN without gaps.
/// The 0-byte is used for filtering out unwanted patterns.
#[rustfmt::skip]
pub(crate) const TO_DNA_ACGTN_WITH_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgtnACGTN",
    b"ACGTNACGTN",
    0,
);

/// Used to convert nucleotide sequences to only uppercase ACGTN without gaps.
/// The 0-byte is used for filtering out unwanted patterns.
#[rustfmt::skip]
pub(crate) const TO_DNA_ACGTN_STD_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgtnACGTN-",
    b"ACGTNACGTN-",
    0,
);

//
// For recode: u8 -> u8 else self
//

/// Used to convert any valid IUPAC DNA to uppercase ACGTN.
pub(crate) const IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC: [u8; 256] = make_mapping_otherwise_self(
    b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-",
    b"ACGTTNNNNNNNNNNNACGTTNNNNNNNNNNN.-",
);

/// Used to convert any valid IUPAC DNA to ACGTN.
pub(crate) const IUPAC_TO_DNA_ACGTN_WITH_GAPS: [u8; 256] = make_mapping_otherwise_self(
    b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-",
    b"acgttnnnnnnnnnnnACGTTNNNNNNNNNNN.-",
);

/// Maps bytes to themselves but valid IUPAC nucleotides to their reverse
/// complement.
pub(crate) const TO_REVERSE_COMPLEMENT: [u8; 256] = make_mapping_otherwise_self(
    b"gcatrykmbvdhuGCATRYKMBVDHU",
    b"cgtayrmkvbhdaCGTAYRMKVBHDA",
);

//
// For recode: u8 -> u8 else default
//

/// Used to convert any byte to uppercase ACGTN. N is used as a catch-all.
#[rustfmt::skip]
pub(crate) const ANY_TO_DNA_ACGTN_NO_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgtunACGTUN",
    b"ACGTTNACGTTN",
    b'N'
);

/// Used to convert any byte to uppercase ACGTN. N is used as a catch-all. Allow
/// gaps.
#[rustfmt::skip]
pub(crate) const ANY_TO_DNA_ACGTN_WITH_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgtunACGTUN.-",
    b"ACGTTNACGTTN.-",
    b'N'
);

/// Used to convert nucleotide sequences to IUPAC nomenclature (no case change),
/// while allowing for gaps and replacing unknown characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_WITH_GAPS: [u8; 256] = make_mapping_with_default(
    b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-",
    b"acgttryswkmbdhvnACGTTRYSWKMBDHVN.-",
    b'N',
);

/// Used to convert nucleotide sequences to upprecase IUPAC nomenclature,
/// while allowing for gaps and replacing unknown characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_WITH_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-",
    b"ACGTTRYSWKMBDHVNACGTTRYSWKMBDHVN.-",
    b'N',
);

/// Used to convert nucleotide sequences to upprecase IUPAC nomenclature, while
/// allowing for and correcting non-standard gaps and replacing unknown
/// characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_CORRECT_GAPS_UC: [u8; 256] = make_mapping_with_default(
    b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-:~",
    b"ACGTTRYSWKMBDHVNACGTTRYSWKMBDHVN.---",
    b'N',
);

//
// For ByteIndexMap: u8 -> usize
//

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T, 4: N}.
/// N is used as a catch-all. U is treated as T.
pub const DNA_PROFILE_MAP: ByteIndexMap<5> =
    ByteIndexMap::new_ignoring_case(*b"ACGTN", b'N').add_synonym_ignoring_case(b'U', b'T');

/// Used by [`ThreeBitKmerEncoder`] to convert any byte to `u8` indices where
/// {3: N, 4: A, 5: C, 6: G, 7: T}. N is used as a catch-all. U is treated as T.
///
/// [`ThreeBitKmerEncoder`]: crate::kmer::ThreeBitKmerEncoder
pub(crate) const THREE_BIT_MAPPING: ByteIndexMap<5> = ByteIndexMap::new_ignoring_case(*b"NACGT", b'N')
    .add_synonym_ignoring_case(b'U', b'T')
    .update_starting_index(3);
