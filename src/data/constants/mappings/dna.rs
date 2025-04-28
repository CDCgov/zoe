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

impl ByteMapBuilder {
    /// Maps `ACGTRYSWKMBDHVN` to themselves, as well as `U` to `T`.
    ///
    /// This does not handle lowercase unless [`handle_case`] was called with
    /// [`FromCase::Any`].
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    const fn map_dna_iupac(self) -> Self {
        self.map_to_self(b"ACGTRYSWKMBDHVN").map(b"U", b"T")
    }

    /// Maps `ACGTN` to themselves, as well as `U` to `T`.
    ///
    /// This does not handle lowercase unless [`handle_case`] was called with
    /// [`FromCase::Any`].
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    const fn map_acgtn(self) -> Self {
        self.map_to_self(b"ACGTN").map(b"U", b"T")
    }

    /// Maps `ACGT` to themselves, as well as `RYSWKMBDHVN` to `N`.
    ///
    /// This does not handle lowercase unless [`handle_case`] was called with
    /// [`FromCase::Any`].
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    const fn map_iupac_to_acgtn(self) -> Self {
        self.map_to_self(b"ACGT").map(b"U", b"T").map_many(b"RYSWKMBDHVN", b'N')
    }

    /// Maps `.` and `-` to themselves.
    const fn with_gaps(self) -> Self {
        self.map_to_self(b".-")
    }
}

//
// Typically used by:   RetainSequence::retain_by_recoding
// Mapping:             u8 -> u8 else 0
//

/// Used to convert nucleotide sequences to valid, uppercase IUPAC DNA without
/// gaps. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_IUPAC_NO_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_dna_iupac()
    .default_to_byte(0)
    .build();

/// Used to convert nucleotide sequences to valid, uppercase IUPAC DNA, allowing
/// gaps. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_IUPAC_WITH_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_dna_iupac()
    .with_gaps()
    .default_to_byte(0)
    .build();

/// Used to convert nucleotide sequences to valid, uppercase IUPAC DNA, also
/// correcting non-standard gaps. The 0-byte is used for filtering out unwanted
/// patterns.
pub(crate) const TO_DNA_IUPAC_CORRECT_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_dna_iupac()
    .map(b".-:~", b".---")
    .default_to_byte(0)
    .build();

/// Used to convert nucleotide sequences to uppercase ACGTN without gaps. The
/// 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_ACGTN_NO_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_acgtn()
    .default_to_byte(0)
    .build();

/// Used to convert nucleotide sequences to uppercase ACGT without gaps. The
/// 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_ACGT_NO_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_to_self(b"ACGT")
    .map(b"U", b"T")
    .default_to_byte(0)
    .build();

/// Used to convert nucleotide sequences to uppercase ACGTN with gaps. The
/// 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_ACGTN_WITH_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_acgtn()
    .with_gaps()
    .default_to_byte(0)
    .build();

/// Used to convert nucleotide sequences to uppercase ACGTN and standardizing
/// gaps to `-`. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_ACGTN_STD_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_acgtn()
    .map(b".-:~", b"----")
    .default_to_byte(0)
    .build();

//
// Typically used by:   Recode::recode
// Mapping:             u8 -> u8 else self
//

/// Used to convert any valid IUPAC DNA to uppercase ACGTN.
pub(crate) const IUPAC_TO_DNA_ACGTN_WITH_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_iupac_to_acgtn()
    .with_gaps()
    .build();

/// Used to convert any valid IUPAC DNA to ACGTN.
pub(crate) const IUPAC_TO_DNA_ACGTN_WITH_GAPS: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Preserve)
    .map_iupac_to_acgtn()
    .with_gaps()
    .build();

/// Maps bytes to themselves but valid IUPAC nucleotides to their reverse
/// complement.
#[rustfmt::skip]
pub(crate) const TO_REVERSE_COMPLEMENT: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Preserve)
    .map(
        b"GCATRYKMBVDHU",
        b"CGTAYRMKVBHDA"
    )
    .build();

//
// Typically used by:   Recode::recode
// Mapping:             u8 -> u8 else default
//

/// Maps bytes to themselves but converts `T` to `U` for converting DNA to RNA.
pub(crate) const TO_RNA: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Preserve)
    .map_to_self(b"ACGUN")
    .map(b"T", b"U")
    .default_to_byte(b'N')
    .build();

/// Used to convert any byte to uppercase RNA. Replaces `T` with `U`.
pub(crate) const TO_RNA_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_to_self(b"ACGUN")
    .map(b"T", b"U")
    .default_to_byte(b'N')
    .build();

/// Used to convert any byte to uppercase ACGTN. N is used as a catch-all.
pub(crate) const ANY_TO_DNA_ACGTN_NO_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_acgtn()
    .default_to_byte(b'N')
    .build();

/// Used to convert any byte to uppercase ACGTN. N is used as a catch-all. Allow
/// gaps.
pub(crate) const ANY_TO_DNA_ACGTN_WITH_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_acgtn()
    .with_gaps()
    .default_to_byte(b'N')
    .build();

/// Used to convert nucleotide sequences to IUPAC nomenclature (no case change),
/// while allowing for gaps and replacing unknown characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_WITH_GAPS: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Preserve)
    .map_dna_iupac()
    .with_gaps()
    .default_to_byte(b'N')
    .build();

/// Used to convert nucleotide sequences to upprecase IUPAC nomenclature,
/// while allowing for gaps and replacing unknown characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_WITH_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_dna_iupac()
    .with_gaps()
    .default_to_byte(b'N')
    .build();

/// Used to convert nucleotide sequences to upprecase IUPAC nomenclature, while
/// allowing for and correcting non-standard gaps and replacing unknown
/// characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_CORRECT_GAPS_UC: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Upper)
    .map_dna_iupac()
    .map(b".-:~", b".---")
    .default_to_byte(b'N')
    .build();

//
// Typically used by:   ByteIndexMap
// Mapping:             u8 -> usize
//

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T, 4: N}.
/// N is used as a catch-all. U is treated as T.
pub const DNA_PROFILE_MAP: ByteIndexMap<5> =
    ByteIndexMap::new_ignoring_case(*b"ACGTN", b'N').add_synonym_ignoring_case(b'U', b'T');

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T}.
/// A is used as a catch-all. U is treated as T.
pub const DNA_UNAMBIG_PROFILE_MAP: ByteIndexMap<4> =
    ByteIndexMap::new_ignoring_case(*b"ACGT", b'A').add_synonym_ignoring_case(b'U', b'T');

/// Used by [`ThreeBitKmerEncoder`] to convert any byte to `u8` indices where
/// {3: N, 4: A, 5: C, 6: G, 7: T}. N is used as a catch-all. U is treated as T.
///
/// [`ThreeBitKmerEncoder`]:
///     crate::kmer::encoders::three_bit::ThreeBitKmerEncoder
pub(crate) const THREE_BIT_MAPPING: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Exact)
    .map(b"ACG", &[4, 5, 6])
    .map_many(b"TU", 7)
    .default_to_byte(3)
    .build();

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T, 4:
/// N, 5: Gap, 6: Other IUPAC, 7: Invalid}. U is treated as T.
pub(crate) const DNA_COUNT_PROFILE_MAP: ByteMap = ByteMapBuilder::new()
    .handle_case(FromCase::Any, ToCase::Exact)
    .map(b"ACG", &[0, 1, 2])
    .map_many(b"TU", 3)
    .map(b"N-", &[4, 5])
    .map_many(b".RYSWKMBDHV", 6)
    .default_to_byte(7)
    .build();
