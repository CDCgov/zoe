#![allow(dead_code)]
/// IUPAC nucleotide codes, including lower and upper case, gaps, missing characters, and ambiguity codes.
pub(crate) const NUCLEIC_IUPAC: &[u8; 34] = b"acgturyswkmbdhvn-.ACGTURYSWKMBDHVN";
/// IUPAC nucleotides excluding gaps and missing character codes.
pub(crate) const NUCLEIC_IUPAC_UNALIGNED: &[u8; 32] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN";

/// Upper and lowercase English alphabet.
pub(crate) const ENGLISH: &[u8; 52] = b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

/// IUPAC Amino Acids in upper case.
pub(crate) const AMINO_ACIDS_UNALIGNED_UC: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// IUPAC Amino Acids.
pub(crate) const AMINO_ACIDS_IUPAC: &[u8; 42] = b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.";
/// DAIS-Ribosome style amino acid codes: IUPAC + gaps + X + partial codons `~`.
pub(crate) const AMINO_ACIDS: &[u8; 45] = b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.~Xx";
