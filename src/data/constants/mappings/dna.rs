use super::*;
use crate::data::{alphas::*, array_types::make_lowercase};

//
// Typically used by:   RetainSequence::retain_by_validation, std::slice::Iter::all
// Mapping:             u8 -> bool
//

/// A boolean mapping of valid IUPAC nucleotide codes without gaps.
pub(crate) const IS_DNA_IUPAC_NO_GAPS: [bool; 256] = make_is_alpha_mapping(DNA_IUPAC_NO_GAPS);

/// A boolean mapping of uppercase IUPAC nucleotide codes without gaps.
pub(crate) const IS_DNA_IUPAC_NO_GAPS_UC: [bool; 256] = make_is_alpha_mapping(DNA_IUPAC_NO_GAPS_UC);

/// A boolean mapping of all valid IUPAC nucleotide codes inluding gaps.
pub(crate) const IS_DNA_IUPAC_WITH_GAPS: [bool; 256] = make_is_alpha_mapping(DNA_IUPAC_WITH_GAPS);

/// A boolean mapping of uppercase IUPAC nucleotide codes including gaps.
pub(crate) const IS_DNA_IUPAC_WITH_GAPS_UC: [bool; 256] = make_is_alpha_mapping(DNA_IUPAC_WITH_GAPS_UC);

/// A boolean mapping of strictly canonical nucleotides + `n`/`N` without gaps.
pub(crate) const IS_DNA_ACGTN_NO_GAPS: [bool; 256] = make_is_alpha_mapping(DNA_ACGTN_NO_GAPS);

/// A boolean mapping of uppercase `ACGTN` bytes without gaps.
pub(crate) const IS_DNA_ACGTN_NO_GAPS_UC: [bool; 256] = make_is_alpha_mapping(DNA_ACGTN_NO_GAPS_UC);

/// A boolean mapping of uppercase `ACGTN-` bytes.
pub(crate) const IS_DNA_ACGTN_STD_GAPS_UC: [bool; 256] = make_is_alpha_mapping(DNA_ACGTN_STD_GAPS_UC);

/// A boolean mapping of `ACGT` bytes without gaps.
pub(crate) const IS_DNA_ACGT_NO_GAPS: [bool; 256] = make_is_alpha_mapping(b"acgtACGT");

/// A boolean mapping of uppercase `ACGT` bytes without gaps.
pub(crate) const IS_DNA_ACGT_NO_GAPS_UC: [bool; 256] = make_is_alpha_mapping(b"ACGT");

//
// Typically used by:   RetainSequence::retain_by_recoding
// Mapping:             u8 -> u8 else 0
//

/// Converts to uppercase IUPAC DNA without gaps. The 0-byte is used for
/// filtering out unwanted patterns.
pub const RETAIN_DNA_IUPAC_NO_GAPS_UC: ByteMap = ByteMap::all(0)
    .preserve(DNA_IUPAC_NO_GAPS_UC)
    .map(&make_lowercase(DNA_IUPAC_NO_GAPS_UC), DNA_IUPAC_NO_GAPS_UC)
    .map(b"Uu", b"TT");

/// Converts to uppercase IUPAC DNA with gaps preserved. The 0-byte is used for
/// filtering out unwanted patterns.
pub const RETAIN_DNA_IUPAC_WITH_GAPS_UC: ByteMap = ByteMap::all(0)
    .preserve(DNA_IUPAC_WITH_GAPS_UC)
    .map(&make_lowercase(DNA_IUPAC_WITH_GAPS_UC), DNA_IUPAC_WITH_GAPS_UC)
    .map(b"Uu", b"TT");

/// Converts to uppercase `ACGTN` without gaps. The 0-byte is used for filtering
/// out unwanted patterns.
pub const RETAIN_DNA_ACGTN_NO_GAPS_UC: ByteMap = ByteMap::all(0)
    .preserve(DNA_ACGTN_NO_GAPS_UC)
    .map(&make_lowercase(DNA_ACGTN_NO_GAPS_UC), DNA_ACGTN_NO_GAPS_UC)
    .map(b"Uu", b"TT");

/// Converts to uppercase `ACGTN` with gaps preserved. The 0-byte is used for
/// filtering out unwanted patterns.
pub const RETAIN_DNA_ACGTN_WITH_GAPS_UC: ByteMap = ByteMap::all(0)
    .preserve(DNA_ACGTN_NO_GAPS_UC)
    .map(&make_lowercase(DNA_ACGTN_NO_GAPS_UC), DNA_ACGTN_NO_GAPS_UC)
    .map(b"Uu", b"TT")
    .preserve(b".-");

/// Converts to uppercase IUPAC DNA with standard gaps preserved and
/// non-standard gaps corrected. The 0-byte is used for filtering out unwanted
/// patterns.
pub(crate) const RETAIN_DNA_IUPAC_CORRECT_GAPS_UC: ByteMap = ByteMap::all(0)
    .preserve(DNA_IUPAC_WITH_GAPS_UC)
    .map(&make_lowercase(DNA_IUPAC_WITH_GAPS_UC), DNA_IUPAC_WITH_GAPS_UC)
    .map(b"Uu", b"TT")
    .map(b":~", b"--");

/// Converts to uppercase `ACGT` without gaps. The 0-byte is used for filtering
/// out unwanted patterns.
pub(crate) const RETAIN_DNA_ACGT_NO_GAPS_UC: ByteMap = ByteMap::all(0).map(b"ACGTacgt", b"ACGTACGT").map(b"Uu", b"TT");

/// Converts to uppercase `ACGTN` with all gaps standardized to `-`. The 0-byte
/// is used for filtering out unwanted patterns.
pub(crate) const RETAIN_DNA_ACGTN_STD_GAPS_UC: ByteMap = ByteMap::all(0)
    .preserve(DNA_ACGTN_NO_GAPS_UC)
    .map(&make_lowercase(DNA_ACGTN_NO_GAPS_UC), DNA_ACGTN_NO_GAPS_UC)
    .map(b"Uu", b"TT")
    .map(b".-:~", b"----");

//
// Typically used by:   Recode::recode
// Mapping:             u8 -> u8 else self
//

/// Converts any IUPAC DNA to uppercase `ACGTN`.
pub(crate) const IUPAC_TO_DNA_ACGTN_UC: ByteMap = ByteMap::identity()
    .map(&make_lowercase(DNA_IUPAC_NO_GAPS_UC), DNA_IUPAC_NO_GAPS_UC)
    .map(b"Uu", b"TT")
    .map_to_one(b"RYSWKMBDHVN", b'N')
    .map_to_one(b"ryswkmbdhvn", b'N');

/// Converts any IUPAC DNA to `ACGTN`.
pub(crate) const IUPAC_TO_DNA_ACGTN: ByteMap = ByteMap::identity()
    .map(b"Uu", b"Tt")
    .map_to_one(b"RYSWKMBDHVN", b'N')
    .map_to_one(b"ryswkmbdhvn", b'n');

/// Maps bytes to themselves but valid IUPAC nucleotides to their complements.
#[rustfmt::skip]
pub(crate) const TO_COMPLEMENT: ByteMap = ByteMap::identity()
    .map(
        b"GCATRYKMBVDHUgcatrykmbvdhu",
        b"CGTAYRMKVBHDAcgtayrmkvbhda"
    );

//
// Typically used by:   Recode::recode
// Mapping:             u8 -> u8 else default
//

/// Converts to uppercase `ACGTN` without gaps. `N` is used for invalid data.
pub const RECODE_TO_DNA_ACGTN_NO_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_ACGTN_NO_GAPS_UC)
    .map(&make_lowercase(DNA_ACGTN_NO_GAPS_UC), DNA_ACGTN_NO_GAPS_UC)
    .map(b"Uu", b"TT");

/// Converts to uppercase `ACGTN` with gaps preserved. `N` is used for invalid
/// data.
pub const RECODE_TO_DNA_ACGTN_WITH_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_ACGTN_NO_GAPS_UC)
    .map(&make_lowercase(DNA_ACGTN_NO_GAPS_UC), DNA_ACGTN_NO_GAPS_UC)
    .map(b"Uu", b"TT")
    .preserve(b".-");

/// Converts to `ACGTN` without gaps while preserving case. Capital `N` is used
/// for invalid data.
pub const RECODE_TO_DNA_ACGTN_NO_GAPS: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_ACGTN_NO_GAPS_UC)
    .preserve(&make_lowercase(DNA_ACGTN_NO_GAPS_UC))
    .map(b"Uu", b"Tt");

/// Converts to `ACGTN` with gaps and case preserved. Capital `N` is used for
/// invalid data.
pub const RECODE_TO_DNA_ACGTN_WITH_GAPS: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_ACGTN_NO_GAPS_UC)
    .preserve(&make_lowercase(DNA_ACGTN_NO_GAPS_UC))
    .map(b"Uu", b"Tt")
    .preserve(b".-");

/// Converts to uppercase IUPAC DNA without gaps. `N` is used for invalid data.
pub const RECODE_TO_DNA_IUPAC_NO_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_IUPAC_NO_GAPS_UC)
    .map(&make_lowercase(DNA_IUPAC_NO_GAPS_UC), DNA_IUPAC_NO_GAPS_UC)
    .map(b"Uu", b"TT");

/// Converts to uppercase IUPAC DNA with gaps preserved. `N` is used for invalid
/// data.
pub const RECODE_TO_DNA_IUPAC_WITH_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_IUPAC_WITH_GAPS_UC)
    .map(&make_lowercase(DNA_IUPAC_WITH_GAPS_UC), DNA_IUPAC_WITH_GAPS_UC)
    .map(b"Uu", b"TT");

/// Converts to IUPAC DNA without gaps while preserving case. Capital `N` is
/// used for invalid data.
pub const RECODE_TO_DNA_IUPAC_NO_GAPS: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_IUPAC_NO_GAPS_UC)
    .preserve(&make_lowercase(DNA_IUPAC_NO_GAPS_UC))
    .map(b"Uu", b"Tt");

/// Converts to IUPAC DNA with gaps and case preserved. Capital `N` is used for
/// invalid data.
pub const RECODE_TO_DNA_IUPAC_WITH_GAPS: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_IUPAC_WITH_GAPS_UC)
    .preserve(&make_lowercase(DNA_IUPAC_WITH_GAPS_UC))
    .map(b"Uu", b"Tt");

/// Converts any byte to RNA while preserving case, replacing `T` with `U`.
/// Capital `N` is used for invalid data.
pub(crate) const TO_RNA: ByteMap = ByteMap::all(b'N').preserve(b"ACGUNacgun").map(b"Tt", b"Uu");

/// Converts any byte to uppercase RNA, replacing `T` with `U`. `N` is used for
/// invalid data.
pub(crate) const TO_RNA_UC: ByteMap = ByteMap::all(b'N').map(b"ACGUT", b"ACGUU").map(b"acgut", b"ACGUU");

/// Converts to upprecase IUPAC DNA, allowing for and correcting non-standard
/// gaps. `N` is used for invalid data.
pub(crate) const RECODE_TO_DNA_IUPAC_CORRECT_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .preserve(DNA_IUPAC_WITH_GAPS_UC)
    .map(&make_lowercase(DNA_IUPAC_WITH_GAPS_UC), DNA_IUPAC_WITH_GAPS_UC)
    .map(b"Uu", b"TT")
    .map(b":~", b"--");

//
// Typically used by:   ByteIndexMap
// Mapping:             u8 -> usize
//

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T, 4: N}.
/// N is used as a catch-all. U is treated as T.
pub const DNA_PROFILE_MAP: ByteIndexMap<5> =
    ByteIndexMap::new_ignoring_case(*b"ACGTN", b'N').add_synonym_ignore_case(b'U', b'T');

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T}.
/// A is used as a catch-all. U is treated as T.
pub const DNA_UNAMBIG_PROFILE_MAP: ByteIndexMap<4> =
    ByteIndexMap::new_ignoring_case(*b"ACGT", b'A').add_synonym_ignore_case(b'U', b'T');

/// Used by [`ThreeBitKmerEncoder`] to convert any byte to `u8` indices where
/// {3: N, 4: A, 5: C, 6: G, 7: T}. N is used as a catch-all. U is treated as T.
///
/// [`ThreeBitKmerEncoder`]:
///     crate::kmer::encoders::three_bit::ThreeBitKmerEncoder
pub(crate) const THREE_BIT_MAPPING: ByteMap = ByteMap::all(3)
    .map(b"ACG", &[4, 5, 6])
    .map(b"acg", &[4, 5, 6])
    .map_to_one(b"TUtu", 7);

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T, 4:
/// N, 5: Gap, 6: Other IUPAC, 7: Invalid}. U is treated as T.
pub(crate) const DNA_COUNT_PROFILE_MAP: ByteMap = ByteMap::all(7)
    .map_to_one(DNA_IUPAC_WITH_GAPS, 6)
    .map(b"ACGTUN-", &[0, 1, 2, 3, 3, 4, 5])
    .map(b"acgtun", &[0, 1, 2, 3, 3, 4]);
