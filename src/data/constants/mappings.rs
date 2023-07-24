#![allow(dead_code)]
use super::alphas::{NUCLEIC_IUPAC, NUCLEIC_IUPAC_UNALIGNED};
use lazy_static::lazy_static;
use std::collections::HashMap;
use std::simd::Simd;

// Nucleotide reverse complement
#[allow(clippy::cast_possible_truncation)]
pub(crate) const TO_REVERSE_COMPLEMENT: [u8; 256] = {
    let mut comp = [0u8; 256];
    let mut i: usize = 0;
    while i < 256 {
        // usize cannot be greater than bound
        comp[i] = i as u8;
        i += 1;
    }

    let forward = b"gcatrykmbvdhuGCATRYKMBVDHU";
    let reverse = b"cgtayrmkvbhdaCGTAYRMKVBHDA";
    let mut b: usize = 0;
    while b < 26 {
        comp[forward[b] as usize] = reverse[b];
        b += 1;
    }

    comp
};
pub(crate) const SIMD64_REVERSE_COMPLEMENT: Simd<u8, 64> = {
    let mut rc: [u8; 64] = [0; 64];
    let mut i: usize = 0;
    while i < rc.len() {
        rc[i] = TO_REVERSE_COMPLEMENT[64 + i];
        i += 1;
    }
    Simd::from_array(rc)
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

macro_rules! fill_map {
    ($( $key: expr => $val: expr ),*) => {{
        let mut map: HashMap<Vec<u8>, u8, ahash::RandomState> = HashMap::default();
        $( map.insert($key.to_vec(), $val); )*
        map
   }}
}

// Assumes the library will uppercase input before hashing.
lazy_static! {
    pub(crate) static ref GENETIC_CODE: HashMap<Vec<u8>, u8, ahash::RandomState> = fill_map!(
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
        b"..."=>b'.', b"---"=>b'-', b"NNN"=>b'X'
    );
}
