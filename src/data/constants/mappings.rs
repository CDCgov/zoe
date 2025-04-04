use crate::data::{
    alphas::{DNA_ACGTN, DNA_ACGTN_UC, DNA_IUPAC, DNA_IUPAC_UC, DNA_IUPAC_WITH_GAPS},
    array_types::{self, arr_max},
};
use std::ops::Index;

use self::array_types::position;

/// Maps bytes to themselves but valid IUPAC nucleotides to their reverse
/// complement.
#[allow(clippy::cast_possible_truncation)]
pub(crate) const TO_REVERSE_COMPLEMENT: [u8; 256] = {
    let mut comp = [0u8; 256];
    let mut i: usize = 0;
    while i < comp.len() {
        // usize cannot be greater than bound
        comp[i] = i as u8;
        i += 1;
    }

    let forward = b"gcatrykmbvdhuGCATRYKMBVDHU";
    let reverse: &[u8; 26] = b"cgtayrmkvbhdaCGTAYRMKVBHDA";
    let mut b: usize = 0;
    while b < forward.len() {
        comp[forward[b] as usize] = reverse[b];
        b += 1;
    }

    comp
};

/// A boolean mapping of all valid IUPAC nucleotide codes. Useful for sequence
/// filtering.
pub(crate) const IS_DNA_IUPAC_WITH_GAPS: [bool; 256] = {
    let mut v = [false; 256];
    let mut i = 0;

    while i < DNA_IUPAC_WITH_GAPS.len() {
        v[DNA_IUPAC_WITH_GAPS[i] as usize] = true;
        i += 1;
    }
    v
};

/// A boolean mapping of valid, unaligned IUPAC nucleotide codes. Useful for
/// sequence filtering.
pub(crate) const IS_DNA_IUPAC: [bool; 256] = {
    let mut v = [false; 256];
    let mut i = 0;

    while i < DNA_IUPAC.len() {
        v[DNA_IUPAC[i] as usize] = true;
        i += 1;
    }
    v
};

/// A boolean mapping of uppercase IUPAC nucleotide codes.
pub(crate) const IS_DNA_IUPAC_UC: [bool; 256] = {
    let mut v = [false; 256];
    let mut i = 0;

    while i < DNA_IUPAC_UC.len() {
        v[DNA_IUPAC_UC[i] as usize] = true;
        i += 1;
    }
    v
};

/// A boolean mapping of canonical nucleotides + `n`/`N`.
pub(crate) const IS_DNA_ACGTN: [bool; 256] = {
    let mut v = [false; 256];
    let mut i = 0;

    while i < DNA_ACGTN.len() {
        v[DNA_ACGTN[i] as usize] = true;
        i += 1;
    }
    v
};

/// A boolean mapping of uppercase canonical nucleotides + `N`.
pub(crate) const IS_DNA_ACGTN_UC: [bool; 256] = {
    let mut v = [false; 256];
    let mut i = 0;

    while i < DNA_ACGTN_UC.len() {
        v[DNA_ACGTN_UC[i] as usize] = true;
        i += 1;
    }
    v
};

/// Used to convert nucleotide sequences to unaligned, IUPAC-validated,
/// uppercase DNA. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_IUPAC_UC: [u8; 256] = {
    const FROM_BYTE: &[u8; 32] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN";
    const DEST_BYTE: &[u8; 32] = b"ACGTTRYSWKMBDHVNACGTTRYSWKMBDHVN";

    let mut v = [0u8; 256];
    let mut i = 0;

    while i < FROM_BYTE.len() {
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

/// Used to convert nucleotide sequences to valid IUPAC DNA (uppercase). The
/// 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_DNA_IUPAC_WITH_GAPS_UC: [u8; 256] = {
    const FROM_BYTE: &[u8; 34] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-";
    const DEST_BYTE: &[u8; 34] = b"ACGTTRYSWKMBDHVNACGTTRYSWKMBDHVN.-";

    let mut v = [0u8; 256];
    let mut i = 0;

    while i < FROM_BYTE.len() {
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

/// Used to convert any valid IUPAC DNA to uppercase ACGTN.
#[allow(clippy::cast_possible_truncation)]
pub(crate) const IUPAC_TO_DNA_ACGTN_UC: [u8; 256] = {
    const FROM_BYTE: &[u8; 32] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN";
    const DEST_BYTE: &[u8; 32] = b"ACGTTNNNNNNNNNNNACGTTNNNNNNNNNNN";

    let mut v = [0u8; 256];
    let mut i = 0;
    while i < v.len() {
        v[i] = i as u8;
        i += 1;
    }

    i = 0;
    while i < FROM_BYTE.len() {
        // reserve 0 for invalid state
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

/// Used to convert any valid IUPAC DNA to ACGTN.
#[allow(clippy::cast_possible_truncation)]
pub(crate) const IUPAC_TO_DNA_ACGTN: [u8; 256] = {
    const FROM_BYTE: &[u8; 32] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN";
    const DEST_BYTE: &[u8; 32] = b"acgttnnnnnnnnnnnACGTTNNNNNNNNNNN";

    let mut v = [0u8; 256];
    let mut i = 0;
    while i < v.len() {
        v[i] = i as u8;
        i += 1;
    }

    i = 0;
    while i < FROM_BYTE.len() {
        // reserve 0 for invalid state
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

/// Used to convert any byte to uppercase ACGTN. N is used as a catch-all.
pub(crate) const ANY_TO_DNA_ACGTN_UC: [u8; 256] = {
    const FROM_BYTE: &[u8; 12] = b"acgtunACGTUN";
    const DEST_BYTE: &[u8; 12] = b"ACGTTNACGTTN";

    // Use N as the catch-all
    let mut v = [b'N'; 256];
    let mut i = 0;

    while i < FROM_BYTE.len() {
        // reserve 0 for invalid state
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

/// Used to convert nucleotide sequences to IUPAC nomenclature (no case change),
/// while allowing for gaps and replacing unknown characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_WITH_GAPS: [u8; 256] = {
    const FROM_BYTE: &[u8; 34] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-";
    const DEST_BYTE: &[u8; 34] = b"acgttryswkmbdhvnACGTTRYSWKMBDHVN.-";

    let mut v = [b'N'; 256];
    let mut i = 0;

    while i < FROM_BYTE.len() {
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

/// Used to convert nucleotide sequences to upprecase IUPAC nomenclature,
/// while allowing for gaps and replacing unknown characters with N.
pub(crate) const ANY_TO_DNA_IUPAC_WITH_GAPS_UC: [u8; 256] = {
    const FROM_BYTE: &[u8; 34] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN.-";
    const DEST_BYTE: &[u8; 34] = b"ACGTTRYSWKMBDHVNACGTTRYSWKMBDHVN.-";

    let mut v = [b'N'; 256];
    let mut i = 0;

    while i < FROM_BYTE.len() {
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

/// Represents a mapping between bytes and indices. For example, this could be a
/// map from DNA bases to profile indices, such as [`DNA_PROFILE_MAP`].
///
/// ## Type Parameters
/// * `KEYS` - The number of bytes being mapped (such as 5 for DNA including
///   *N*)
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct ByteIndexMap<const KEYS: usize> {
    pub(crate) index_map: [u8; 256],
    pub(crate) byte_keys: [u8; KEYS],
}

impl<const S: usize> ByteIndexMap<S> {
    /// Create a new [`ByteIndexMap`] struct to represent a mapping between
    /// bytes and indices. For example, this could be a map from DNA bases to
    /// profile indices. Any byte that is not specified in `byte_keys` is mapped
    /// to the same thing as `catch_all`.
    ///
    /// ## Panics
    /// No duplicates can be present in `byte_keys`. `catch_all` must be present
    /// in `byte_keys`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new(byte_keys: [u8; S], catch_all: u8) -> Self {
        assert!(array_types::is_unique(&byte_keys));

        let catch_all_index =
            position(&byte_keys, catch_all).expect("The catch_all must be present in the byte_keys.") as u8;
        let mut out = ByteIndexMap {
            index_map: [catch_all_index; 256],
            byte_keys,
        };

        let mut i = 0;
        while i < byte_keys.len() {
            // Truncation will not occur because i cannot exceed
            // byte_keys.len(), and index must contain unique u8 values
            out.set_byte(byte_keys[i], i as u8);

            i += 1;
        }
        out
    }

    /// Create a new [`ByteIndexMap`] struct to represent a mapping between
    /// bytes and indices. For example, this could be a map from DNA bases to
    /// profile indices. Both `byte_keys` and `catch_all` ignore case.
    ///
    /// ## Panics
    /// No duplicates can be present in `byte_keys`. `catch_all` must be present
    /// in `byte_keys`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new_ignoring_case(mut byte_keys: [u8; S], catch_all: u8) -> Self {
        byte_keys = array_types::make_uppercase(&byte_keys);
        assert!(array_types::is_unique(&byte_keys));

        let catch_all_index = position(&byte_keys, catch_all.to_ascii_uppercase())
            .expect("The catch_all must be present in the byte_keys.") as u8;
        let mut out = ByteIndexMap {
            index_map: [catch_all_index; 256],
            byte_keys,
        };

        let mut i = 0;
        while i < byte_keys.len() {
            // Truncation will not occur because i cannot exceed
            // byte_keys.len(), and byte_keys must contain unique u8 values
            out.set_byte_ignoring_case(byte_keys[i], i as u8);

            i += 1;
        }
        out
    }

    /// Increment all values in the map to support a different starting index.
    const fn update_starting_index(mut self, starting_index: u8) -> Self {
        assert!(arr_max(&self.index_map).unwrap() < u8::MAX - starting_index);
        let mut i = 0;
        while i < self.index_map.len() {
            self.index_map[i] += starting_index;
            i += 1;
        }
        self
    }

    /// Set the index for a byte.
    #[inline]
    const fn set_byte(&mut self, byte: u8, index: u8) {
        self.index_map[byte as usize] = index;
    }

    /// Set the index for a byte, ignoring case.
    #[inline]
    const fn set_byte_ignoring_case(&mut self, byte: u8, index: u8) {
        self.index_map[byte.to_ascii_lowercase() as usize] = index;
        self.index_map[byte.to_ascii_uppercase() as usize] = index;
    }

    /// Change the [`ByteIndexMap`] so that `new_byte` maps to the same thing as
    /// `byte_key`.
    #[inline]
    #[must_use]
    pub const fn add_synonym(mut self, new_key: u8, previous_key: u8) -> Self {
        self.set_byte(new_key, self.copy_index(previous_key));
        self
    }

    /// Change the [`ByteIndexMap`] so that `new_byte` maps to the same thing as
    /// `byte_key`, ignoring case.
    #[inline]
    #[must_use]
    pub const fn add_synonym_ignoring_case(mut self, new_key: u8, previous_key: u8) -> Self {
        self.set_byte_ignoring_case(new_key, self.copy_index(previous_key));
        self
    }

    /// Get the length of `byte_keys`.
    #[inline]
    #[must_use]
    #[allow(clippy::len_without_is_empty)]
    pub const fn len(&self) -> usize {
        self.byte_keys.len()
    }

    /// Convert a base `b` into an index.
    #[inline]
    #[must_use]
    pub const fn to_index(&self, b: u8) -> usize {
        self.index_map[b as usize] as usize
    }

    #[inline]
    #[must_use]
    const fn copy_index(&self, b: u8) -> u8 {
        self.index_map[b as usize]
    }
}

impl<const S: usize> Index<u8> for ByteIndexMap<S> {
    type Output = u8;

    #[inline]
    fn index(&self, index: u8) -> &u8 {
        &self.index_map[index as usize]
    }
}

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

macro_rules! fill_std_gc {
    ($( $key: expr => $val: expr ),*) => {{
        let mut map = StdGeneticCode::new_raw_table();
        $( StdGeneticCode::insert(&mut map, $key, $val); )*
        map
   }}
}

/// A custom, perfect hashmap for translating resolvable codons to amino acids
/// under the _standard_ genetic code..
///
/// If a codon involving ambiguous IUPAC letters translates to the same amino
/// acid in all cases, then it is also included in the hashmap. Stop codons are
/// translated to `*`, and the codons `...`, `---`, and `NNN` are translated to
/// `.`, `-`, and `X` respectively.
#[derive(Debug)]
pub struct StdGeneticCode;
static STD_GENETIC_CODE: [u32; StdGeneticCode::TABLE_SIZE] = fill_std_gc!(
    b"TAA"=>b'*', b"TAG"=>b'*', b"TAR"=>b'*', b"TGA"=>b'*', b"TRA"=>b'*', b"GCA"=>b'A', b"GCB"=>b'A', b"GCC"=>b'A', b"GCD"=>b'A', b"GCG"=>b'A', b"GCH"=>b'A',
    b"GCK"=>b'A', b"GCM"=>b'A', b"GCN"=>b'A', b"GCR"=>b'A', b"GCS"=>b'A', b"GCT"=>b'A', b"GCV"=>b'A', b"GCW"=>b'A', b"GCY"=>b'A', b"TGC"=>b'C', b"TGT"=>b'C',
    b"TGY"=>b'C', b"GAC"=>b'D', b"GAT"=>b'D', b"GAY"=>b'D', b"GAA"=>b'E', b"GAG"=>b'E', b"GAR"=>b'E', b"TTC"=>b'F', b"TTT"=>b'F', b"TTY"=>b'F', b"GGA"=>b'G',
    b"GGB"=>b'G', b"GGC"=>b'G', b"GGD"=>b'G', b"GGG"=>b'G', b"GGH"=>b'G', b"GGK"=>b'G', b"GGM"=>b'G', b"GGN"=>b'G', b"GGR"=>b'G', b"GGS"=>b'G', b"GGT"=>b'G',
    b"GGV"=>b'G', b"GGW"=>b'G', b"GGY"=>b'G', b"CAC"=>b'H', b"CAT"=>b'H', b"CAY"=>b'H', b"ATA"=>b'I', b"ATC"=>b'I', b"ATH"=>b'I', b"ATM"=>b'I', b"ATT"=>b'I',
    b"ATW"=>b'I', b"ATY"=>b'I', b"AAA"=>b'K', b"AAG"=>b'K', b"AAR"=>b'K', b"CTA"=>b'L', b"CTB"=>b'L', b"CTC"=>b'L', b"CTD"=>b'L', b"CTG"=>b'L', b"CTH"=>b'L',
    b"CTK"=>b'L', b"CTM"=>b'L', b"CTN"=>b'L', b"CTR"=>b'L', b"CTS"=>b'L', b"CTT"=>b'L', b"CTV"=>b'L', b"CTW"=>b'L', b"CTY"=>b'L', b"TTA"=>b'L', b"TTG"=>b'L',
    b"TTR"=>b'L', b"YTA"=>b'L', b"YTG"=>b'L', b"YTR"=>b'L', b"ATG"=>b'M', b"AAC"=>b'N', b"AAT"=>b'N', b"AAY"=>b'N', b"CCA"=>b'P', b"CCB"=>b'P', b"CCC"=>b'P',
    b"CCD"=>b'P', b"CCG"=>b'P', b"CCH"=>b'P', b"CCK"=>b'P', b"CCM"=>b'P', b"CCN"=>b'P', b"CCR"=>b'P', b"CCS"=>b'P', b"CCT"=>b'P', b"CCV"=>b'P', b"CCW"=>b'P',
    b"CCY"=>b'P', b"CAA"=>b'Q', b"CAG"=>b'Q', b"CAR"=>b'Q', b"AGA"=>b'R', b"AGG"=>b'R', b"AGR"=>b'R', b"CGA"=>b'R', b"CGB"=>b'R', b"CGC"=>b'R', b"CGD"=>b'R',
    b"CGG"=>b'R', b"CGH"=>b'R', b"CGK"=>b'R', b"CGM"=>b'R', b"CGN"=>b'R', b"CGR"=>b'R', b"CGS"=>b'R', b"CGT"=>b'R', b"CGV"=>b'R', b"CGW"=>b'R', b"CGY"=>b'R',
    b"MGA"=>b'R', b"MGG"=>b'R', b"MGR"=>b'R', b"AGC"=>b'S', b"AGT"=>b'S', b"AGY"=>b'S', b"TCA"=>b'S', b"TCB"=>b'S', b"TCC"=>b'S', b"TCD"=>b'S', b"TCG"=>b'S',
    b"TCH"=>b'S', b"TCK"=>b'S', b"TCM"=>b'S', b"TCN"=>b'S', b"TCR"=>b'S', b"TCS"=>b'S', b"TCT"=>b'S', b"TCV"=>b'S', b"TCW"=>b'S', b"TCY"=>b'S', b"ACA"=>b'T',
    b"ACB"=>b'T', b"ACC"=>b'T', b"ACD"=>b'T', b"ACG"=>b'T', b"ACH"=>b'T', b"ACK"=>b'T', b"ACM"=>b'T', b"ACN"=>b'T', b"ACR"=>b'T', b"ACS"=>b'T', b"ACT"=>b'T',
    b"ACV"=>b'T', b"ACW"=>b'T', b"ACY"=>b'T', b"GTA"=>b'V', b"GTB"=>b'V', b"GTC"=>b'V', b"GTD"=>b'V', b"GTG"=>b'V', b"GTH"=>b'V', b"GTK"=>b'V', b"GTM"=>b'V',
    b"GTN"=>b'V', b"GTR"=>b'V', b"GTS"=>b'V', b"GTT"=>b'V', b"GTV"=>b'V', b"GTW"=>b'V', b"GTY"=>b'V', b"TGG"=>b'W', b"TAC"=>b'Y', b"TAT"=>b'Y', b"TAY"=>b'Y',
    b"..."=>b'.', b"---"=>b'-', b"NNN"=>b'X',

    b"UAA"=>b'*', b"UAG"=>b'*', b"UAR"=>b'*', b"UGA"=>b'*', b"URA"=>b'*', b"GCU"=>b'A', b"UGC"=>b'C', b"UGU"=>b'C', b"UGY"=>b'C', b"GAU"=>b'D', b"UUC"=>b'F',
    b"UUU"=>b'F', b"UUY"=>b'F', b"GGU"=>b'G', b"CAU"=>b'H', b"AUA"=>b'I', b"AUC"=>b'I', b"AUH"=>b'I', b"AUM"=>b'I', b"AUU"=>b'I', b"AUW"=>b'I', b"AUY"=>b'I',
    b"CUA"=>b'L', b"CUB"=>b'L', b"CUC"=>b'L', b"CUD"=>b'L', b"CUG"=>b'L', b"CUH"=>b'L', b"CUK"=>b'L', b"CUM"=>b'L', b"CUN"=>b'L', b"CUR"=>b'L', b"CUS"=>b'L',
    b"CUU"=>b'L', b"CUV"=>b'L', b"CUW"=>b'L', b"CUY"=>b'L', b"UUA"=>b'L', b"UUG"=>b'L', b"UUR"=>b'L', b"YUA"=>b'L', b"YUG"=>b'L', b"YUR"=>b'L', b"AUG"=>b'M',
    b"AAU"=>b'N', b"CCU"=>b'P', b"CGU"=>b'R', b"AGU"=>b'S', b"UCA"=>b'S', b"UCB"=>b'S', b"UCC"=>b'S', b"UCD"=>b'S', b"UCG"=>b'S', b"UCH"=>b'S', b"UCK"=>b'S',
    b"UCM"=>b'S', b"UCN"=>b'S', b"UCR"=>b'S', b"UCS"=>b'S', b"UCU"=>b'S', b"UCV"=>b'S', b"UCW"=>b'S', b"UCY"=>b'S', b"ACU"=>b'T', b"GUA"=>b'V', b"GUB"=>b'V',
    b"GUC"=>b'V', b"GUD"=>b'V', b"GUG"=>b'V', b"GUH"=>b'V', b"GUK"=>b'V', b"GUM"=>b'V', b"GUN"=>b'V', b"GUR"=>b'V', b"GUS"=>b'V', b"GUU"=>b'V', b"GUV"=>b'V',
    b"GUW"=>b'V', b"GUY"=>b'V', b"UGG"=>b'W', b"UAC"=>b'Y', b"UAU"=>b'Y', b"UAY"=>b'Y'
);

impl StdGeneticCode {
    const TABLE_SIZE: usize = 2048;
    const MASK: u32 = 0xFF_FF_FF;
    const TABLE: &[u32; Self::TABLE_SIZE] = &STD_GENETIC_CODE;

    /// Retrieve the amino acid corresponding to a codon, or return `None` if
    /// the translation cannot be resolved. To instead return 'X', see
    /// [`translate_codon`](StdGeneticCode::translate_codon).
    ///
    /// ## Panics
    ///
    /// Panics if codon has fewer than three elements. Any additional elements
    /// past the first three are ignored.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    #[inline]
    pub const fn get(codon: &[u8]) -> Option<u8> {
        let key = Self::to_upper_u32(codon);
        let found = StdGeneticCode::TABLE[Self::get_index(key)];

        let key_found = found & Self::MASK;
        let aa_found = found >> 24;

        if key_found == key {
            // The first packed byte does not truncate
            Some(aa_found as u8)
        } else {
            None
        }
    }

    /// Retrieve the amino acid corresponding to a codon, or return `'X'` if the
    /// translation cannot be resolved. To instead return an `Option`, see
    /// [`get`](StdGeneticCode::get).
    ///
    /// ## Panics
    ///
    /// Panics if codon has fewer than three elements. Any additional elements
    /// past the first three are ignored.
    #[inline]
    #[must_use]
    pub fn translate_codon(codon: &[u8]) -> u8 {
        Self::get(codon).unwrap_or(b'X')
    }

    /// Determine whether the provided codon is a stop codon. Any additional
    /// elements past the first three are ignored and partial codons are
    /// considered `false`.
    #[must_use]
    #[inline]
    pub fn is_stop_codon(c: &[u8]) -> bool {
        c.len() > 2
            && matches!(
                &[
                    c[0].to_ascii_uppercase(),
                    c[1].to_ascii_uppercase(),
                    c[2].to_ascii_uppercase()
                ],
                b"TAA" | b"TAG" | b"TAR" | b"TGA" | b"TRA" | b"UAA" | b"UAG" | b"UAR" | b"UGA" | b"URA"
            )
    }

    #[must_use]
    const fn new_raw_table() -> [u32; Self::TABLE_SIZE] {
        let mut raw_table = [0u32; Self::TABLE_SIZE];
        // Initializing the table with 0 makes it appear that key 0 is in the
        // hashmap with value 0. Since no valid codons have an index of 0 we
        // must invalidate the key at index 0 to guard against false positives
        // (consider a trivial codon input b"\0\0\0" as an example). We do this
        // by setting the key to 1 at index 0.
        //
        // This change leaves one more issue. A false entry with key of 1 must
        // not map to index 0 under our current hash function. This assertion
        // ensures any future change to the hash function will not violate this
        // invariant.
        assert!(Self::get_index(1) != 0);
        raw_table[0] = 1;
        raw_table
    }

    #[inline]
    const fn insert(raw_table: &mut [u32; Self::TABLE_SIZE], codon: &[u8], aa: u8) {
        let num = Self::to_upper_u32(codon);
        let index = Self::get_index(num);

        assert!(raw_table[index] == 0, "A collision occurred between codons.",);

        raw_table[index] = ((aa as u32) << 24) | num;
    }

    #[must_use]
    #[inline]
    const fn get_index(num: u32) -> usize {
        let hash = num ^ (num >> 9) ^ (num >> 12);
        // The compiler seems to be smart enough to optimize this
        hash as usize % Self::TABLE_SIZE
    }

    #[must_use]
    #[inline]
    const fn to_upper_u32(codon: &[u8]) -> u32 {
        (codon[2].to_ascii_uppercase() as u32)
            | ((codon[1].to_ascii_uppercase() as u32) << 8)
            | ((codon[0].to_ascii_uppercase() as u32) << 16)
    }
}

#[cfg(test)]
mod test {
    use crate::data::{
        alphas::{DNA_ACGTN, DNA_ACGTN_UC, DNA_IUPAC, DNA_IUPAC_UC, DNA_IUPAC_WITH_GAPS},
        mappings::{
            ANY_TO_DNA_ACGTN_UC, ANY_TO_DNA_IUPAC_WITH_GAPS, ANY_TO_DNA_IUPAC_WITH_GAPS_UC, IS_DNA_IUPAC_WITH_GAPS,
            IUPAC_TO_DNA_ACGTN, IUPAC_TO_DNA_ACGTN_UC, StdGeneticCode,
        },
    };

    use super::{
        ByteIndexMap, DNA_PROFILE_MAP, IS_DNA_ACGTN, IS_DNA_ACGTN_UC, IS_DNA_IUPAC, IS_DNA_IUPAC_UC, TO_DNA_IUPAC_UC,
        TO_DNA_IUPAC_WITH_GAPS_UC,
    };

    #[test]
    fn gc_self_test() {
        #[rustfmt::skip]
        let gc = [
            (b"TAA", b'*'), (b"TAG", b'*'), (b"TAR", b'*'), (b"TGA", b'*'), (b"TRA", b'*'), (b"GCA", b'A'), (b"GCB", b'A'), (b"GCC", b'A'), (b"GCD", b'A'), (b"GCG", b'A'), (b"GCH", b'A'),
            (b"GCK", b'A'), (b"GCM", b'A'), (b"GCN", b'A'), (b"GCR", b'A'), (b"GCS", b'A'), (b"GCT", b'A'), (b"GCV", b'A'), (b"GCW", b'A'), (b"GCY", b'A'), (b"TGC", b'C'), (b"TGT", b'C'),
            (b"TGY", b'C'), (b"GAC", b'D'), (b"GAT", b'D'), (b"GAY", b'D'), (b"GAA", b'E'), (b"GAG", b'E'), (b"GAR", b'E'), (b"TTC", b'F'), (b"TTT", b'F'), (b"TTY", b'F'), (b"GGA", b'G'),
            (b"GGB", b'G'), (b"GGC", b'G'), (b"GGD", b'G'), (b"GGG", b'G'), (b"GGH", b'G'), (b"GGK", b'G'), (b"GGM", b'G'), (b"GGN", b'G'), (b"GGR", b'G'), (b"GGS", b'G'), (b"GGT", b'G'),
            (b"GGV", b'G'), (b"GGW", b'G'), (b"GGY", b'G'), (b"CAC", b'H'), (b"CAT", b'H'), (b"CAY", b'H'), (b"ATA", b'I'), (b"ATC", b'I'), (b"ATH", b'I'), (b"ATM", b'I'), (b"ATT", b'I'),
            (b"ATW", b'I'), (b"ATY", b'I'), (b"AAA", b'K'), (b"AAG", b'K'), (b"AAR", b'K'), (b"CTA", b'L'), (b"CTB", b'L'), (b"CTC", b'L'), (b"CTD", b'L'), (b"CTG", b'L'), (b"CTH", b'L'),
            (b"CTK", b'L'), (b"CTM", b'L'), (b"CTN", b'L'), (b"CTR", b'L'), (b"CTS", b'L'), (b"CTT", b'L'), (b"CTV", b'L'), (b"CTW", b'L'), (b"CTY", b'L'), (b"TTA", b'L'), (b"TTG", b'L'),
            (b"TTR", b'L'), (b"YTA", b'L'), (b"YTG", b'L'), (b"YTR", b'L'), (b"ATG", b'M'), (b"AAC", b'N'), (b"AAT", b'N'), (b"AAY", b'N'), (b"CCA", b'P'), (b"CCB", b'P'), (b"CCC", b'P'),
            (b"CCD", b'P'), (b"CCG", b'P'), (b"CCH", b'P'), (b"CCK", b'P'), (b"CCM", b'P'), (b"CCN", b'P'), (b"CCR", b'P'), (b"CCS", b'P'), (b"CCT", b'P'), (b"CCV", b'P'), (b"CCW", b'P'),
            (b"CCY", b'P'), (b"CAA", b'Q'), (b"CAG", b'Q'), (b"CAR", b'Q'), (b"AGA", b'R'), (b"AGG", b'R'), (b"AGR", b'R'), (b"CGA", b'R'), (b"CGB", b'R'), (b"CGC", b'R'), (b"CGD", b'R'),
            (b"CGG", b'R'), (b"CGH", b'R'), (b"CGK", b'R'), (b"CGM", b'R'), (b"CGN", b'R'), (b"CGR", b'R'), (b"CGS", b'R'), (b"CGT", b'R'), (b"CGV", b'R'), (b"CGW", b'R'), (b"CGY", b'R'),
            (b"MGA", b'R'), (b"MGG", b'R'), (b"MGR", b'R'), (b"AGC", b'S'), (b"AGT", b'S'), (b"AGY", b'S'), (b"TCA", b'S'), (b"TCB", b'S'), (b"TCC", b'S'), (b"TCD", b'S'), (b"TCG", b'S'),
            (b"TCH", b'S'), (b"TCK", b'S'), (b"TCM", b'S'), (b"TCN", b'S'), (b"TCR", b'S'), (b"TCS", b'S'), (b"TCT", b'S'), (b"TCV", b'S'), (b"TCW", b'S'), (b"TCY", b'S'), (b"ACA", b'T'),
            (b"ACB", b'T'), (b"ACC", b'T'), (b"ACD", b'T'), (b"ACG", b'T'), (b"ACH", b'T'), (b"ACK", b'T'), (b"ACM", b'T'), (b"ACN", b'T'), (b"ACR", b'T'), (b"ACS", b'T'), (b"ACT", b'T'),
            (b"ACV", b'T'), (b"ACW", b'T'), (b"ACY", b'T'), (b"GTA", b'V'), (b"GTB", b'V'), (b"GTC", b'V'), (b"GTD", b'V'), (b"GTG", b'V'), (b"GTH", b'V'), (b"GTK", b'V'), (b"GTM", b'V'),
            (b"GTN", b'V'), (b"GTR", b'V'), (b"GTS", b'V'), (b"GTT", b'V'), (b"GTV", b'V'), (b"GTW", b'V'), (b"GTY", b'V'), (b"TGG", b'W'), (b"TAC", b'Y'), (b"TAT", b'Y'), (b"TAY", b'Y'),
            (b"...", b'.'), (b"---", b'-'), (b"NNN", b'X'),

            (b"UAA", b'*'), (b"UAG", b'*'), (b"UAR", b'*'), (b"UGA", b'*'), (b"URA", b'*'), (b"GCU", b'A'), (b"UGC", b'C'), (b"UGU", b'C'), (b"UGY", b'C'), (b"GAU", b'D'), (b"UUC", b'F'),
            (b"UUU", b'F'), (b"UUY", b'F'), (b"GGU", b'G'), (b"CAU", b'H'), (b"AUA", b'I'), (b"AUC", b'I'), (b"AUH", b'I'), (b"AUM", b'I'), (b"AUU", b'I'), (b"AUW", b'I'), (b"AUY", b'I'),
            (b"CUA", b'L'), (b"CUB", b'L'), (b"CUC", b'L'), (b"CUD", b'L'), (b"CUG", b'L'), (b"CUH", b'L'), (b"CUK", b'L'), (b"CUM", b'L'), (b"CUN", b'L'), (b"CUR", b'L'), (b"CUS", b'L'),
            (b"CUU", b'L'), (b"CUV", b'L'), (b"CUW", b'L'), (b"CUY", b'L'), (b"UUA", b'L'), (b"UUG", b'L'), (b"UUR", b'L'), (b"YUA", b'L'), (b"YUG", b'L'), (b"YUR", b'L'), (b"AUG", b'M'),
            (b"AAU", b'N'), (b"CCU", b'P'), (b"CGU", b'R'), (b"AGU", b'S'), (b"UCA", b'S'), (b"UCB", b'S'), (b"UCC", b'S'), (b"UCD", b'S'), (b"UCG", b'S'), (b"UCH", b'S'), (b"UCK", b'S'),
            (b"UCM", b'S'), (b"UCN", b'S'), (b"UCR", b'S'), (b"UCS", b'S'), (b"UCU", b'S'), (b"UCV", b'S'), (b"UCW", b'S'), (b"UCY", b'S'), (b"ACU", b'T'), (b"GUA", b'V'), (b"GUB", b'V'),
            (b"GUC", b'V'), (b"GUD", b'V'), (b"GUG", b'V'), (b"GUH", b'V'), (b"GUK", b'V'), (b"GUM", b'V'), (b"GUN", b'V'), (b"GUR", b'V'), (b"GUS", b'V'), (b"GUU", b'V'), (b"GUV", b'V'),
            (b"GUW", b'V'), (b"GUY", b'V'), (b"UGG", b'W'), (b"UAC", b'Y'), (b"UAU", b'Y'), (b"UAY", b'Y')
        ];

        for (codon, aa) in gc {
            assert_eq!(aa, StdGeneticCode::get(codon).unwrap());
            assert_eq!(aa, StdGeneticCode::translate_codon(codon));
            assert_eq!(aa == b'*', StdGeneticCode::is_stop_codon(codon));
        }
    }

    #[test]
    fn test_dna_map() {
        for i in 0..=255 {
            match i {
                b'A' | b'a' => assert!(DNA_PROFILE_MAP.to_index(i) == 0),
                b'C' | b'c' => assert!(DNA_PROFILE_MAP.to_index(i) == 1),
                b'G' | b'g' => assert!(DNA_PROFILE_MAP.to_index(i) == 2),
                b'T' | b't' | b'U' | b'u' => assert!(DNA_PROFILE_MAP.to_index(i) == 3),
                _ => assert!(DNA_PROFILE_MAP.to_index(i) == 4),
            }
        }
    }

    #[test]
    fn test_dna_map_ignores_case() {
        const MAP1: ByteIndexMap<5> = ByteIndexMap::new_ignoring_case(*b"acgtn", b'N').add_synonym_ignoring_case(b'u', b'T');
        const MAP2: ByteIndexMap<5> = ByteIndexMap::new_ignoring_case(*b"AcGtN", b'n').add_synonym_ignoring_case(b'U', b't');
        assert_eq!(DNA_PROFILE_MAP, MAP1);
        assert_eq!(DNA_PROFILE_MAP, MAP2);
    }

    #[test]
    fn test_def_a_map() {
        const DEF_A_MAP: ByteIndexMap<4> =
            ByteIndexMap::new_ignoring_case(*b"ACGT", b'A').add_synonym_ignoring_case(b'U', b'T');
        for i in 0..=255 {
            match i {
                b'C' | b'c' => assert!(DEF_A_MAP.to_index(i) == 1),
                b'G' | b'g' => assert!(DEF_A_MAP.to_index(i) == 2),
                b'T' | b't' | b'U' | b'u' => assert!(DEF_A_MAP.to_index(i) == 3),
                _ => assert!(DEF_A_MAP.to_index(i) == 0),
            }
        }
    }

    #[test]
    fn test_case_sensitive() {
        const DNA_MAP: ByteIndexMap<10> = ByteIndexMap::new(*b"ACGTNacgtn", b'N')
            .add_synonym(b'U', b'T')
            .add_synonym(b'u', b't');
        for i in 0..=255 {
            match i {
                b'A' => assert!(DNA_MAP.to_index(i) == 0),
                b'C' => assert!(DNA_MAP.to_index(i) == 1),
                b'G' => assert!(DNA_MAP.to_index(i) == 2),
                b'T' | b'U' => assert!(DNA_MAP.to_index(i) == 3),
                b'a' => assert!(DNA_MAP.to_index(i) == 5),
                b'c' => assert!(DNA_MAP.to_index(i) == 6),
                b'g' => assert!(DNA_MAP.to_index(i) == 7),
                b't' | b'u' => assert!(DNA_MAP.to_index(i) == 8),
                b'n' => assert!(DNA_MAP.to_index(i) == 9),
                _ => assert!(DNA_MAP.to_index(i) == 4),
            }
        }
    }

    #[test]
    #[should_panic = "assertion failed: array_types::is_unique(&byte_keys)"]
    fn test_duplicate() {
        let _ = ByteIndexMap::new(*b"ACGTNA", b'N');
    }

    #[test]
    #[should_panic = "assertion failed: array_types::is_unique(&byte_keys)"]
    fn test_duplicate_nocase() {
        let _ = ByteIndexMap::new_ignoring_case(*b"ACGTNA", b'N');
    }

    #[test]
    #[should_panic = "The catch_all must be present in the byte_keys."]
    fn test_missing_catch_all() {
        let _ = ByteIndexMap::new(*b"ACGT", b'N');
    }

    #[test]
    #[should_panic = "The catch_all must be present in the byte_keys."]
    fn test_missing_catch_all_nocase() {
        let _ = ByteIndexMap::new_ignoring_case(*b"ACGT", b'N');
    }

    const VALIDATOR_PAIRS: [(&[u8], [bool; 256]); 5] = [
        (DNA_IUPAC_WITH_GAPS, IS_DNA_IUPAC_WITH_GAPS),
        (DNA_IUPAC, IS_DNA_IUPAC),
        (DNA_IUPAC_UC, IS_DNA_IUPAC_UC),
        (DNA_ACGTN, IS_DNA_ACGTN),
        (DNA_ACGTN_UC, IS_DNA_ACGTN_UC),
    ];

    #[test]
    fn test_alphabet_validators() {
        for (alphabet, validator) in VALIDATOR_PAIRS {
            for c in u8::MIN..=u8::MAX {
                assert_eq!(alphabet.contains(&c), validator[c as usize]);
            }
        }
    }

    const TO_DNA_PAIRS: [(&[u8], [u8; 256]); 2] =
        [(DNA_IUPAC, TO_DNA_IUPAC_UC), (DNA_IUPAC_WITH_GAPS, TO_DNA_IUPAC_WITH_GAPS_UC)];

    #[test]
    fn test_to_dna() {
        for (alphabet, converter) in TO_DNA_PAIRS {
            for c in u8::MIN..=u8::MAX {
                let in_alphabet = alphabet.contains(&c);
                if in_alphabet {
                    if matches!(c, b'u' | b'U') {
                        assert_eq!(converter[c as usize], b'T');
                    } else if c.is_ascii_lowercase() {
                        assert_eq!(converter[c as usize], c.to_ascii_uppercase());
                    } else {
                        assert_eq!(converter[c as usize], c);
                    }
                } else {
                    assert_eq!(converter[c as usize], 0);
                }
            }
        }
    }

    #[test]
    fn test_iupac_to_dna_acgtn_uc() {
        for c in u8::MIN..=u8::MAX {
            if matches!(c, b'u' | b'U') {
                assert_eq!(IUPAC_TO_DNA_ACGTN_UC[c as usize], b'T');
            } else if DNA_ACGTN.contains(&c) {
                assert_eq!(IUPAC_TO_DNA_ACGTN_UC[c as usize], c.to_ascii_uppercase());
            } else if DNA_IUPAC.contains(&c) {
                assert_eq!(IUPAC_TO_DNA_ACGTN_UC[c as usize], b'N');
            } else {
                assert_eq!(IUPAC_TO_DNA_ACGTN_UC[c as usize], c);
            }
        }
    }

    #[test]
    fn test_iupac_to_dna_acgtn() {
        for c in u8::MIN..=u8::MAX {
            if c == b'u' {
                assert_eq!(IUPAC_TO_DNA_ACGTN[c as usize], b't');
            } else if c == b'U' {
                assert_eq!(IUPAC_TO_DNA_ACGTN[c as usize], b'T');
            } else if DNA_ACGTN.contains(&c) {
                assert_eq!(IUPAC_TO_DNA_ACGTN[c as usize], c);
            } else if DNA_IUPAC.contains(&c) {
                if c.is_ascii_lowercase() {
                    assert_eq!(IUPAC_TO_DNA_ACGTN[c as usize], b'n');
                } else {
                    assert_eq!(IUPAC_TO_DNA_ACGTN[c as usize], b'N');
                }
            } else {
                assert_eq!(IUPAC_TO_DNA_ACGTN[c as usize], c);
            }
        }
    }

    #[test]
    fn test_any_to_dna_acgtn_uc() {
        for c in u8::MIN..=u8::MAX {
            if matches!(c, b'u' | b'U') {
                assert_eq!(ANY_TO_DNA_ACGTN_UC[c as usize], b'T');
            } else if DNA_ACGTN.contains(&c) {
                assert_eq!(ANY_TO_DNA_ACGTN_UC[c as usize], c.to_ascii_uppercase());
            } else {
                assert_eq!(ANY_TO_DNA_ACGTN_UC[c as usize], b'N');
            }
        }
    }

    #[test]
    fn test_any_to_dna_iupac_with_gaps() {
        for c in u8::MIN..=u8::MAX {
            if c == b'u' {
                assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS[c as usize], b't');
            } else if c == b'U' {
                assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS[c as usize], b'T');
            } else if DNA_IUPAC_WITH_GAPS.contains(&c) {
                assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS[c as usize], c);
            } else {
                assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS[c as usize], b'N');
            }
        }
    }

    #[test]
    fn test_any_to_dna_iupac_with_gaps_uc() {
        for c in u8::MIN..=u8::MAX {
            if matches!(c, b'u' | b'U') {
                assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS_UC[c as usize], b'T');
            } else if DNA_IUPAC_WITH_GAPS.contains(&c) {
                assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS_UC[c as usize], c.to_ascii_uppercase());
            } else {
                assert_eq!(ANY_TO_DNA_IUPAC_WITH_GAPS_UC[c as usize], b'N');
            }
        }
    }
}
