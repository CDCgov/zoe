#![allow(dead_code)]

/// IUPAC nucleotides without gaps.
pub(crate) const DNA_IUPAC: &[u8; 32] = b"acgturyswkmbdhvnACGTURYSWKMBDHVN";
/// Uppercase IUPAC nucleotides without gaps.
pub(crate) const DNA_IUPAC_UC: &[u8; 16] = b"ACGTURYSWKMBDHVN";

/// Canonical nucleotide bases + `n` or `N`, no gaps.
pub(crate) const DNA_ACGTN: &[u8; 10] = b"acgtnACGTN";
/// Canonical uppercase nucleotide bases + `N`, no gaps.
pub(crate) const DNA_ACGTN_UC: &[u8; 5] = b"ACGTN";

/// IUPAC Amino Acids in upper case.
pub(crate) const AA_IUPAC_UC: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// IUPAC Amino Acids in upper case, with `X` for ambiguous/unknown residues.
pub(crate) const AA_IUPAC_UC_X: &[u8; 21] = b"ACDEFGHIKLMNPQRSTVWYX";

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use super::*;

    const ALPHABETS: [&[u8]; 6] = [DNA_IUPAC, DNA_IUPAC_UC, DNA_ACGTN, DNA_ACGTN_UC, AA_IUPAC_UC, AA_IUPAC_UC_X];

    const CAP_PAIRS: [(&[u8], &[u8]); 2] = [(DNA_IUPAC, DNA_IUPAC_UC), (DNA_ACGTN, DNA_ACGTN_UC)];

    const SUBSET_PAIRS: [(&[u8], &[u8]); 3] = [
        (DNA_ACGTN_UC, DNA_IUPAC_UC),
        (DNA_ACGTN, DNA_IUPAC),
        (AA_IUPAC_UC, AA_IUPAC_UC_X),
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
