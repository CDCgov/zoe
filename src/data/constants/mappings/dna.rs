use super::*;
use crate::data::alphas::*;

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

/// Used to convert nucleotide sequences to valid, uppercase IUPAC DNA without
/// gaps. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_IUPAC_NO_GAPS_UC: ByteMap = ByteMap::all(0)
    .map_ignore_case(DNA_IUPAC_NO_GAPS_UC, DNA_IUPAC_NO_GAPS_UC)
    .map_ignore_case(b"U", b"T");

/// Used to convert nucleotide sequences to valid, uppercase IUPAC DNA, allowing
/// gaps. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_IUPAC_WITH_GAPS_UC: ByteMap = ByteMap::all(0)
    .map_ignore_case(DNA_IUPAC_WITH_GAPS_UC, DNA_IUPAC_WITH_GAPS_UC)
    .map_ignore_case(b"U", b"T");

/// Used to convert nucleotide sequences to valid, uppercase IUPAC DNA, also
/// correcting non-standard gaps. The 0-byte is used for filtering out unwanted
/// patterns.
pub(crate) const TO_DNA_IUPAC_CORRECT_GAPS_UC: ByteMap = ByteMap::all(0)
    .map_ignore_case(DNA_IUPAC_WITH_GAPS_UC, DNA_IUPAC_WITH_GAPS_UC)
    .map_ignore_case(b"U", b"T")
    .map(b":~", b"--");

/// Used to convert nucleotide sequences to uppercase ACGTN without gaps. The
/// 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_ACGTN_NO_GAPS_UC: ByteMap = ByteMap::all(0)
    .map_ignore_case(DNA_ACGTN_NO_GAPS_UC, DNA_ACGTN_NO_GAPS_UC)
    .map_ignore_case(b"U", b"T");

/// Used to convert nucleotide sequences to uppercase ACGT without gaps. The
/// 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_ACGT_NO_GAPS_UC: ByteMap =
    ByteMap::all(0).map_ignore_case(b"ACGT", b"ACGT").map_ignore_case(b"U", b"T");

/// Used to convert nucleotide sequences to uppercase ACGTN with gaps. The
/// 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_ACGTN_WITH_GAPS_UC: ByteMap = ByteMap::all(0)
    .map_ignore_case(DNA_ACGTN_NO_GAPS_UC, DNA_ACGTN_NO_GAPS_UC)
    .map_ignore_case(b"U", b"T")
    .preserve(b".-");

/// Used to convert nucleotide sequences to uppercase ACGTN and standardizing
/// gaps to `-`. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_ACGTN_STD_GAPS_UC: ByteMap = ByteMap::all(0)
    .map_ignore_case(DNA_ACGTN_NO_GAPS_UC, DNA_ACGTN_NO_GAPS_UC)
    .map_ignore_case(b"U", b"T")
    .map(b".-:~", b"----");

//
// Typically used by:   Recode::recode
// Mapping:             u8 -> u8 else self
//

/// Used to convert any valid IUPAC DNA to uppercase ACGTN.
pub(crate) const IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC: ByteMap = ByteMap::identity()
    .map_ignore_case(DNA_IUPAC_NO_GAPS_UC, DNA_IUPAC_NO_GAPS_UC)
    .map_ignore_case(b"U", b"T")
    .many_to_one_ignore_case(b"RYSWKMBDHVN", b'N');

/// Used to convert any valid IUPAC DNA to ACGTN.
pub(crate) const IUPAC_TO_DNA_ACGTN_WITH_GAPS: ByteMap = ByteMap::identity()
    .preserve_both_cases(DNA_IUPAC_NO_GAPS_UC)
    .map(b"Uu", b"Tt")
    .many_to_one(b"RYSWKMBDHVN", b'N')
    .many_to_one(b"ryswkmbdhvn", b'n');

/// Maps bytes to themselves but valid IUPAC nucleotides to their reverse
/// complement.
#[rustfmt::skip]
pub(crate) const TO_REVERSE_COMPLEMENT: ByteMap = ByteMap::identity()
    .map(
        b"GCATRYKMBVDHUgcatrykmbvdhu",
        b"CGTAYRMKVBHDAcgtayrmkvbhda"
    );

//
// Typically used by:   Recode::recode
// Mapping:             u8 -> u8 else default
//

/// Maps bytes to themselves but converts `T` to `U` for converting DNA to RNA.
pub(crate) const TO_RNA: ByteMap = ByteMap::all(b'N').preserve_both_cases(b"ACGUN").map(b"Tt", b"Uu");

/// Used to convert any byte to uppercase RNA. Replaces `T` with `U`.
pub(crate) const TO_RNA_UC: ByteMap = ByteMap::all(b'N').map_ignore_case(b"ACGUT", b"ACGUU");

/// Used to convert any byte to uppercase ACGTN. N is used as a catch-all.
pub(crate) const ANY_TO_DNA_ACGTN_NO_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .map_ignore_case(DNA_ACGTN_NO_GAPS_UC, DNA_ACGTN_NO_GAPS_UC)
    .map_ignore_case(b"U", b"T");

/// Used to convert any byte to uppercase ACGTN. N is used as a catch-all. Allow
/// gaps.
pub(crate) const ANY_TO_DNA_ACGTN_WITH_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .map_ignore_case(DNA_ACGTN_NO_GAPS_UC, DNA_ACGTN_NO_GAPS_UC)
    .map_ignore_case(b"U", b"T")
    .preserve(b".-");

/// Used to convert nucleotide sequences to IUPAC nomenclature (no case change),
/// while allowing for gaps and replacing unknown characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_WITH_GAPS: ByteMap = ByteMap::all(b'N')
    .preserve_both_cases(DNA_IUPAC_WITH_GAPS_UC)
    .map(b"Uu", b"Tt");

/// Used to convert nucleotide sequences to upprecase IUPAC nomenclature,
/// while allowing for gaps and replacing unknown characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_WITH_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .map_ignore_case(DNA_IUPAC_WITH_GAPS_UC, DNA_IUPAC_WITH_GAPS_UC)
    .map_ignore_case(b"U", b"T");

/// Used to convert nucleotide sequences to upprecase IUPAC nomenclature, while
/// allowing for and correcting non-standard gaps and replacing unknown
/// characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_CORRECT_GAPS_UC: ByteMap = ByteMap::all(b'N')
    .map_ignore_case(DNA_IUPAC_WITH_GAPS_UC, DNA_IUPAC_WITH_GAPS_UC)
    .map_ignore_case(b"U", b"T")
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
    .map_ignore_case(b"ACG", &[4, 5, 6])
    .many_to_one_ignore_case(b"TU", 7);

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T, 4:
/// N, 5: Gap, 6: Other IUPAC, 7: Invalid}. U is treated as T.
pub(crate) const DNA_COUNT_PROFILE_MAP: ByteMap = ByteMap::all(7)
    .many_to_one_ignore_case(DNA_IUPAC_WITH_GAPS_UC, 6)
    .map_ignore_case(b"ACGTUN-", &[0, 1, 2, 3, 3, 4, 5]);
