#![allow(dead_code)]
/// IUPAC nucleotide codes, including lower and upper case, gaps, missing
/// characters, and ambiguity codes.
pub(crate) const DNA_IUPAC_WITH_GAPS: &[u8; 34] = b"acgturyswkmbdhvn-.ACGTURYSWKMBDHVN";
/// IUPAC nucleotides excluding gaps and missing character codes.
pub(crate) const DNA_IUPAC: &[u8; 32] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN";
/// Uppercase IUPAC nucleotide bases, including ambiguity codes.
pub(crate) const DNA_IUPAC_UC: &[u8; 16] = b"ACGTURYSWKMBDHVN";
/// Canonical nucleotide bases + `n` or `N`.
pub(crate) const DNA_ACGTN: &[u8; 10] = b"acgtnACGTN";
/// Canonical uppercase nucleotide bases + `N`.
pub(crate) const DNA_ACGTN_UC: &[u8; 5] = b"ACGTN";

/// Upper and lowercase English alphabet.
pub(crate) const ENGLISH: &[u8; 52] = b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

/// IUPAC Amino Acids in upper case.
pub(crate) const AA_IUPAC_UC: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// IUPAC Amino Acids.
pub(crate) const AA_IUPAC_WITH_GAPS: &[u8; 42] = b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.";
/// DAIS-Ribosome style amino acid codes: IUPAC + gaps + X + partial codons `~`.
pub(crate) const AA_DAIS_WITH_GAPS: &[u8; 45] = b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.~Xx";

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use super::*;

    const ALPHABETS: [&[u8]; 9] = [
        DNA_IUPAC_WITH_GAPS,
        DNA_IUPAC,
        DNA_IUPAC_UC,
        DNA_ACGTN,
        DNA_ACGTN_UC,
        ENGLISH,
        AA_IUPAC_UC,
        AA_IUPAC_WITH_GAPS,
        AA_DAIS_WITH_GAPS,
    ];

    const CAP_PAIRS: [(&[u8], &[u8]); 2] = [(DNA_ACGTN, DNA_ACGTN_UC), (DNA_IUPAC, DNA_IUPAC_UC)];

    const SUBSET_PAIRS: [(&[u8], &[u8]); 4] = [
        (DNA_ACGTN, DNA_IUPAC),
        (DNA_IUPAC, DNA_IUPAC_WITH_GAPS),
        (AA_IUPAC_UC, AA_IUPAC_WITH_GAPS),
        (AA_IUPAC_WITH_GAPS, AA_DAIS_WITH_GAPS),
    ];

    #[test]
    fn no_duplicates() {
        for alphabet in ALPHABETS {
            let mut set = HashSet::new();
            for c in alphabet {
                assert!(set.insert(c));
            }
        }
    }

    #[test]
    fn caps() {
        for (mixed, upper) in CAP_PAIRS {
            let mixed_set = mixed.iter().copied().collect::<HashSet<_>>();
            for c in upper {
                assert!(mixed_set.contains(c));

                if c.is_ascii_alphabetic() {
                    assert!(c.is_ascii_uppercase());
                    assert!(mixed_set.contains(&c.to_ascii_lowercase()));
                }
            }
        }
    }

    #[test]
    fn subsets() {
        for (smaller, larger) in SUBSET_PAIRS {
            let larger_set = larger.iter().copied().collect::<HashSet<_>>();
            for c in smaller {
                assert!(larger_set.contains(c));
            }
        }
    }
}
