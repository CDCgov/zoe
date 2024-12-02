use super::alphas::{NUCLEIC_IUPAC, NUCLEIC_IUPAC_UNALIGNED};
use std::{ops::Index, sync::LazyLock};

/// Maps bytes to themselves but valid IUPAC nucleotides to their reverse complement.
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

/// A boolean mapping of all valid IUPAC nucleotide codes. Useful for
/// sequence filtering.
pub(crate) const IS_IUPAC_BASE: [bool; 256] = {
    let mut v = [false; 256];
    let mut i = 0;

    while i < NUCLEIC_IUPAC.len() {
        v[NUCLEIC_IUPAC[i] as usize] = true;
        i += 1;
    }
    v
};

/// A boolean mapping of valid, unaligned IUPAC nucleotide codes. Useful for
/// sequence filtering.
pub(crate) const IS_UNALIGNED_IUPAC_BASE: [bool; 256] = {
    let mut v = [false; 256];
    let mut i = 0;

    while i < NUCLEIC_IUPAC_UNALIGNED.len() {
        v[NUCLEIC_IUPAC_UNALIGNED[i] as usize] = true;
        i += 1;
    }
    v
};

/// Used to convert nucleotide sequences to unaligned, IUPAC-validated,
/// uppercase DNA. The 0-byte is used for filtering out unwanted patterns.
pub(crate) const TO_UNALIGNED_DNA_UC: [u8; 256] = {
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
pub(crate) const TO_DNA_UC: [u8; 256] = {
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

#[allow(clippy::cast_possible_truncation)]
/// Used to convert any valid IUPAC DNA to uppercase ACGTN.
pub(crate) const IUPAC_TO_DNA_CANONICAL_UPPER: [u8; 256] = {
    const FROM_BYTE: &[u8; 32] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN";
    const DEST_BYTE: &[u8; 32] = b"ACGTTNNNNNNNNNNNACGTTNNNNNNNNNNN";

    let mut v = [0u8; 256];
    let mut i = 0;
    while i < v.len() {
        v[i] = i as u8;
        i += 1;
    }

    while i < FROM_BYTE.len() {
        // reserve 0 for invalid state
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

#[allow(clippy::cast_possible_truncation)]
/// Used to convert any valid IUPAC DNA to ACGTN.
pub(crate) const IUPAC_TO_DNA_CANONICAL: [u8; 256] = {
    const FROM_BYTE: &[u8; 32] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN";
    const DEST_BYTE: &[u8; 32] = b"acgttnnnnnnnnnnnACGTTNNNNNNNNNNN";

    let mut v = [0u8; 256];
    let mut i = 0;
    while i < v.len() {
        v[i] = i as u8;
        i += 1;
    }

    while i < FROM_BYTE.len() {
        // reserve 0 for invalid state
        v[FROM_BYTE[i] as usize] = DEST_BYTE[i];
        i += 1;
    }
    v
};

/// Used to convert any byte to uppercase ACGTN. N is used as a catch-all.
pub(crate) const ANY_TO_DNA_CANONICAL_UPPER: [u8; 256] = {
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

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct ResidueMapping<const S: usize> {
    pub(crate) mapping: [u8; 256],
    pub(crate) index:   [u8; S],
}

impl<const S: usize> ResidueMapping<S> {
    /// Create a new [`ResidueMapping`] struct to represent a mapping between
    /// residues and indices. For DNA bases, you may consider using
    /// [`DNA_RESIDUE_MAPPING`], which uniquely represents `ACGTN` with `N`
    /// being a catch-all.
    ///
    /// # Panics
    /// Uppercase ASCII is expected for the `index` byte string and the
    /// `catch_all` value. The `catch_all` value must be present in `index`, and
    /// no duplicates can be present in `index`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new(index: [u8; S], catch_all: u8) -> Self {
        assert!(catch_all.is_ascii_uppercase());

        let mut catch_all_index = None;
        let mut i = 0;
        while i < index.len() {
            assert!(index[i].is_ascii_uppercase());

            let mut j = i + 1;
            while j < index.len() {
                assert!(index[i] != index[j]);
                j += 1;
            }

            if index[i] == catch_all {
                catch_all_index = Some(i);
            }

            i += 1;
        }

        // If this unwrap fails, it means catch_all wasn't present in index
        let mut mapping = [catch_all_index.unwrap() as u8; 256];

        let mut i = 0;
        while i < index.len() {
            // Truncation will not occur because i cannot exceed index.len(),
            // and index must contain unique u8 values
            mapping[index[i] as usize] = i as u8;
            mapping[index[i].to_ascii_lowercase() as usize] = i as u8;

            i += 1;
        }
        ResidueMapping { mapping, index }
    }

    #[inline]
    #[must_use]
    #[allow(clippy::len_without_is_empty)]
    pub const fn len(&self) -> usize {
        self.index.len()
    }

    #[inline]
    #[must_use]
    pub const fn get_index(&self, b: u8) -> usize {
        self.mapping[b as usize] as usize
    }
}

impl<const S: usize> Index<u8> for ResidueMapping<S> {
    type Output = u8;

    #[inline]
    fn index(&self, index: u8) -> &u8 {
        &self.mapping[index as usize]
    }
}

/// Used to convert any byte to `u8` indices where {0: A, 1: C, 2: G, 3: T, 4: N}.
/// N is used as a catch-all.
pub const DNA_RESIDUE_MAPPING: ResidueMapping<5> = ResidueMapping::new(*b"ACGTN", b'N');

macro_rules! fill_map {
    ($( $key: expr => $val: expr ),*) => {{
        let mut map = CodonTranslator::new();
        $( map.insert( $key, $val); )*
        map
   }}
}

/// A custom, perfect hashmap for translating resolvable codons to amino acids.
///
/// If a codon involving ambiguous IUPAC letters translates to the same amino
/// acid in all cases, then it is also included in the hashmap. Stop codons are
/// translated to `*`, and the codons `...`, `---`, and `NNN` are translated to
/// `.`, `-`, and `X` respectively.
#[derive(Debug)]
pub struct CodonTranslator {
    table: Vec<u32>,
}

impl CodonTranslator {
    const TABLE_SIZE: usize = 2048;
    const MASK: u32 = 0xFF_FF_FF;

    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    #[inline]
    /// Retrieve the amino acid corresponding to a codon, or return None if the
    /// translation cannot be resolved.
    ///
    /// # Panics
    ///
    /// Panics if codon has fewer than three elements. Any additional elements
    /// past the first three are ignored.
    pub fn get(&self, codon: &[u8]) -> Option<u8> {
        let key = Self::to_upper_u32(codon);
        let found = self.table[Self::get_index(key)];

        let key_found = found & Self::MASK;
        let aa_found = found >> 24;

        if key_found == key {
            // The first packed byte does not truncate
            Some(aa_found as u8)
        } else {
            None
        }
    }

    #[must_use]
    fn new() -> Self {
        let mut raw_table = vec![0u32; Self::TABLE_SIZE];
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
        CodonTranslator { table: raw_table }
    }

    #[inline]
    fn insert(&mut self, codon: &[u8], aa: u8) {
        let num = Self::to_upper_u32(codon);
        let index = Self::get_index(num);
        assert!(
            self.table[index] == 0,
            "Codon {c} @ {index} already contains {aa2}, trying to insert {aa}",
            c = std::str::from_utf8(codon).unwrap(),
            aa2 = (self.table[index] >> 24) as u8 as char,
            aa = aa as char,
        );
        self.table[index] = u32::from(aa) << 24 | num;
    }

    #[must_use]
    #[inline]
    fn get_index(num: u32) -> usize {
        let hash = num ^ (num >> 9) ^ (num >> 12);
        // The compiler seems to be smart enough to optimize this
        hash as usize % Self::TABLE_SIZE
    }

    #[must_use]
    #[inline]
    fn to_upper_u32(codon: &[u8]) -> u32 {
        u32::from(codon[2].to_ascii_uppercase())
            | u32::from(codon[1].to_ascii_uppercase()) << 8
            | u32::from(codon[0].to_ascii_uppercase()) << 16
    }
}

/// A custom, perfect hashmap for the _standard_ Genetic Code where ambiguous
/// codons can be unambiguously translated.
pub(crate) static GENETIC_CODE: LazyLock<CodonTranslator> = LazyLock::new(|| {
    fill_map!(
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
    )
});

#[cfg(test)]
mod test {
    use super::GENETIC_CODE;

    #[test]
    fn gc_self_test() {
        #[rustfmt::skip]
        let gc = vec![
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
            assert_eq!(aa, GENETIC_CODE.get(codon).unwrap());
        }
    }
}
