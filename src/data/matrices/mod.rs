//! ## Substitution/Scoring Matrices for Alignment
//!
//! *Zoe* includes several constant [`WeightMatrix`] values for commonly used
//! protein substitution matrices:
//!
//! - PAM matrices: [PAM
//!   matrices](https://en.wikipedia.org/wiki/Point_accepted_mutation) (or Point
//!   Accepted Mutation matrices) were originally introduced by Margaret Dayhoff
//!   in 1978 (1). The PAM$\text{}_{\text{n}}$ matrix is formed via
//!   a Markov assumption by raising a transition matrix to the $n$th power, and
//!   it models the substitution probabilities given that for every 100 amino
//!   acids there are $n$ mutations (possibly at the same location). The
//!   transition matrix was based on 1572 observed mutations in closely related
//!   proteins.
//! - BLOSUM matrices: [BLOSUM matrices](https://en.wikipedia.org/wiki/BLOSUM)
//!   were introduced in 1992 (2). They are first formed by
//!   clustering a database of aligned proteins such that the sequences all have
//!   less than $r$% sequence similarity. Using this data, the log odds of
//!   observing a pair of amino acids is computed, which is then scaled and
//!   rounded to obtain the score.
//!
//! A more in-depth description can be found
//! [here](https://web.archive.org/web/20240527025425/https://cs.rice.edu/~ogilvie/comp571/pam-vs-blosum/).
//!
//! *Zoe* provides compile-time functions for generating PAM matrices based on
//! NCBI's [`pam.tar.gz`](https://ftp.ncbi.nih.gov/blast/matrices/), with
//! ambiguous amino acids handled in a similar way to
//! [Biopython](https://github.com/biopython/biopython/tree/master/Bio/Align/substitution_matrices/data).
//! Additionally, several common PAM matrices are directly defined as well.
//!
//! *Zoe*'s includes files for common BLOSUM matrices, generated using a
//! modified version of NCBI's
//! [`blosum.tar.Z`](https://ftp.ncbi.nih.gov/repository/blocks/unix/blosum/).
//!
//! ## Choice of Matrix
//!
//! PAM and BLOSUM are considered similar, although they are based on different
//! assumptions. PAM matrices assume a Markov property which allows
//! PAM$\text{}_{\text{n}}$ to be derived by raising a transition matrix to a
//! power. The BLOSUM matrices implicitly model this by clustering the data by
//! sequence similarity beforehand. We recommend BLOSUM as a general-use scoring
//! matrix, although testing both and comparing the results may be worthwhile.
//!
//! The following matrices are considered similar (based on entropy and expected
//! score):
//!
//! | PAM    | BLOSUM   |
//! |--------|----------|
//! | PAM100 | BLOSUM90 |
//! | PAM120 | BLOSUM80 |
//! | PAM160 | BLOSUM62 |
//! | PAM200 | BLOSUM50 |
//! | PAM250 | BLOSUM45 |
//!
//! Higher BLOSUM numbers (and lower PAM numbers) are used for closely related
//! proteins. Lower BLOSUM numbers (and higher PAM numbers) are used for
//! distantly related proteins.
//!
//! ## Scale
//!
//! The choice of scale for BLOSUM matrices is relatively standard, chosen based
//! on the entropy of the frequencies of the pairs. However, different sources
//! for BLOSUM80 and BLOSUM100 may switch between 1/2 and 1/3 bit units (a
//! scaling of 2 and 3 respectively). *Zoe* elects to use 1/3 bit units to
//! increase the precision with which the scores are represented.
//!
//! For PAM matrices, *Zoe* uses $\ln(2)/2$ for PAM30 and PAM70 and $\ln(2)/3$
//! for PAM250 (similar to
//! [Biopython](https://github.com/biopython/biopython/tree/master/Bio/Align/substitution_matrices/data)).
//!
//! ## Ambiguous Amino Acids
//!
//! *Zoe* includes the four ambiguous amino acid codes `B`, `Z`, `J`, and `X`.
//! See [IUPAC Standards](crate::data#iupac-standards) for more discussion. The
//! original formulations of the PAM and BLOSUM matrices do not describe how to
//! handle ambiguous amino acids.
//!
//! [`blosum.tar.Z`](https://ftp.ncbi.nih.gov/repository/blocks/unix/blosum/)
//! chooses to perform a weighted average of the log odds scores, but their code
//! has a bug for computing the score of `(X, Z)` and a seemingly incorrect
//! formula for `(B, Z)`. *Zoe* corrects both, which may make some of the scores
//! decrease by 1.
//!
//! Note that *Zoe* also uses a weighted average for computing comparisons
//! against `X`, which may result in different scores depending on the other
//! amino acid. This is different from some sources (such as
//! [here](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62))
//! which use the same fixed penalty.
//!
//! [`pam.tar.Z`](https://ftp.ncbi.nih.gov/blast/matrices/) does not include
//! ambiguous amino acids. To match existing results in
//! [Biopython](https://github.com/biopython/biopython/tree/master/Bio/Align/substitution_matrices/data),
//! *Zoe* adds the joint probabilities for ambiguous amino acids (excluding `X`)
//! before computing the log odds. For `X`, a weighted average of the log odds
//! is performed instead.
//!
//! ## Stop Codon
//!
//! For the stop codon, *Zoe* gives a score of 1 for matching stop codons and
//! the minimum score in the BLOSUM matrix for a stop codon aligned to a
//! residue.
//!
//! ## Module Citations
//!
//! 1. Dayhoff, M., Schwartz, R., & Orcutt, B. (1978). "A model of evolutionary
//!    change in proteins". In M. Dayhoff (Ed.), Atlas of Protein Sequence and
//!    Structure (Vol. 5, pp. 345–352). Washington, D. C.: National Biomedical
//!    Research Foundation.
//!
//! 2. Henikoff, S., & Henikoff, J. G. (1992). "Amino acid substitution matrices
//!    from protein blocks". Proceedings of the National Academy of Sciences of
//!    the United States of America, 89(22), 10915–10919. doi:
//!    <https://doi.org/10.1073/pnas.89.22.10915>

use crate::{
    data::{
        array_types::{elem_max, is_subset},
        mappings::{ByteIndexMap, DNA_PROFILE_MAP},
    },
    math::AnyInt,
};
use std::fmt::Display;

pub(crate) mod aa;
mod parse;

pub use aa::*;
pub use parse::*;

/// Physiochemical distance matrix using the euclidean distance between all
/// amino acid factors.
pub(crate) static PHYSIOCHEMICAL_FACTORS: [[Option<f32>; 256]; 256] = {
    const AA: [u8; 43] = [
        b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R', b'S', b'T', b'V', b'W',
        b'Y', b'X', b'a', b'c', b'd', b'e', b'f', b'g', b'h', b'i', b'k', b'l', b'm', b'n', b'p', b'q', b'r', b's', b't',
        b'v', b'w', b'y', b'x', b'-',
    ];

    /// # Citation
    /// For factor analysis used by the function:
    ///
    /// > Atchley et al. 2008. "Solving the protein sequence metric problem."
    /// > Proc Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Published 2005
    /// > Apr 25.
    const PCF: [[f64; 5]; 43] = [
        [-0.59, -1.3, -0.73, 1.57, -0.15],
        [-1.34, 0.47, -0.86, -1.02, -0.26],
        [1.05, 0.3, -3.66, -0.26, -3.24],
        [1.36, -1.45, 1.48, 0.11, -0.84],
        [-1.01, -0.59, 1.89, -0.4, 0.41],
        [-0.38, 1.65, 1.33, 1.05, 2.06],
        [0.34, -0.42, -1.67, -1.47, -0.08],
        [-1.24, -0.55, 2.13, 0.39, 0.82],
        [1.83, -0.56, 0.53, -0.28, 1.65],
        [-1.02, -0.99, -1.51, 1.27, -0.91],
        [-0.66, -1.52, 2.22, -1.01, 1.21],
        [0.95, 0.83, 1.3, -0.17, 0.93],
        [0.19, 2.08, -1.63, 0.42, -1.39],
        [0.93, -0.18, -3.01, -0.5, -1.85],
        [1.54, -0.06, 1.5, 0.44, 2.9],
        [-0.23, 1.4, -4.76, 0.67, -2.65],
        [-0.03, 0.33, 2.21, 0.91, 1.31],
        [-1.34, -0.28, -0.54, 1.24, -1.26],
        [-0.6, 0.01, 0.67, -2.13, -0.18],
        [0.26, 0.83, 3.1, -0.84, 1.51],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.59, -1.3, -0.73, 1.57, -0.15],
        [-1.34, 0.47, -0.86, -1.02, -0.26],
        [1.05, 0.3, -3.66, -0.26, -3.24],
        [1.36, -1.45, 1.48, 0.11, -0.84],
        [-1.01, -0.59, 1.89, -0.4, 0.41],
        [-0.38, 1.65, 1.33, 1.05, 2.06],
        [0.34, -0.42, -1.67, -1.47, -0.08],
        [-1.24, -0.55, 2.13, 0.39, 0.82],
        [1.83, -0.56, 0.53, -0.28, 1.65],
        [-1.02, -0.99, -1.51, 1.27, -0.91],
        [-0.66, -1.52, 2.22, -1.01, 1.21],
        [0.95, 0.83, 1.3, -0.17, 0.93],
        [0.19, 2.08, -1.63, 0.42, -1.39],
        [0.93, -0.18, -3.01, -0.5, -1.85],
        [1.54, -0.06, 1.5, 0.44, 2.9],
        [-0.23, 1.4, -4.76, 0.67, -2.65],
        [-0.03, 0.33, 2.21, 0.91, 1.31],
        [-1.34, -0.28, -0.54, 1.24, -1.26],
        [-0.6, 0.01, 0.67, -2.13, -0.18],
        [0.26, 0.83, 3.1, -0.84, 1.51],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let mut pcd = [[None; 256]; 256];

    let mut aa1: usize = 0;
    while aa1 < 43 {
        let mut aa2: usize = 0;
        while aa2 < 43 {
            pcd[AA[aa1] as usize][AA[aa2] as usize] = if AA[aa1].eq_ignore_ascii_case(&AA[aa2]) {
                Some(0.0)
            } else {
                let mut d: f64 = 0.0;
                let mut k: usize = 0;
                while k < 5 {
                    d += (PCF[aa1][k] - PCF[aa2][k]) * (PCF[aa1][k] - PCF[aa2][k]);
                    k += 1;
                }
                #[allow(clippy::cast_possible_truncation)]
                Some(crate::math::sqrt_baby(d) as f32)
            };

            aa2 += 1;
        }
        aa1 += 1;
    }

    pcd
};

/// A weight matrix representing the scores for various matches and mismatches
/// when performing sequence alignment.
///
/// Often, a signed weight matrix is needed, in which case `T` should be `i8`.
/// To construct a new signed weight matrix, use [`new`] or [`new_dna_matrix`].
///
/// Unsigned weight matrices are also supported, in which case `T` should be
/// `u8`. To construct a new unsigned weight matrix, first create a signed
/// weight matrix, then call [`to_biased_matrix`]. This will shift all entries
/// in the matrix to become nonnegative, and then store the resulting shift. One
/// can also use [`new_biased_dna_matrix`].
///
/// [`new`]: WeightMatrix::new
/// [`new_dna_matrix`]: WeightMatrix::new_dna_matrix
/// [`new_biased_dna_matrix`]: WeightMatrix::new_biased_dna_matrix
/// [`to_biased_matrix`]: WeightMatrix::to_biased_matrix
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct WeightMatrix<'a, T: AnyInt, const S: usize> {
    pub weights:     [[T; S]; S],
    pub mapping:     &'a ByteIndexMap<S>,
    pub(crate) bias: T,
}

impl<T: AnyInt, const S: usize> WeightMatrix<'_, T, S> {
    /// For a given `ref_residue` and `query_residue`, retrieves the weight
    /// stored in the matrix.
    #[inline]
    #[must_use]
    pub const fn get_weight(&self, ref_residue: u8, query_residue: u8) -> T {
        self.weights[self.mapping.to_index(ref_residue)][self.mapping.to_index(query_residue)]
    }

    /// For a given `ref_residue` and `query_residue`, retrieves a mutable
    /// reference to the weight stored in the matrix.
    #[inline]
    #[must_use]
    const fn get_weight_mut(&mut self, ref_residue: u8, query_residue: u8) -> &mut T {
        &mut self.weights[self.mapping.to_index(ref_residue)][self.mapping.to_index(query_residue)]
    }

    /// Subsets the rows and columns of a weight matrix to only include the
    /// bytes in `subset`, without updating or recomputing the bias.
    const fn get_subset_unadjusted<'a, const L: usize>(&self, subset: &'a ByteIndexMap<L>) -> WeightMatrix<'a, T, L> {
        assert!(
            is_subset(&self.mapping.byte_keys, &subset.byte_keys),
            "The provided keys are not a subset of the original weight matrix!"
        );
        let mut weights = WeightMatrix {
            weights: [[T::ZERO; L]; L],
            mapping: subset,
            bias:    self.bias,
        };
        let mut i = 0;
        while i < L {
            let reference_base = subset.byte_keys[i];
            let mut j = 0;
            while j < L {
                let query_base = subset.byte_keys[j];
                *weights.get_weight_mut(reference_base, query_base) = self.get_weight(reference_base, query_base);
                j += 1;
            }
            i += 1;
        }
        weights
    }
}

impl WeightMatrix<'_, u8, 5> {
    /// Creates a new [`WeightMatrix`] with a fixed `matching` score,
    /// `mismatch` score, and optionally ignoring a base. A pair of bases where
    /// either is the ignored base will always have a score of 0.
    #[inline]
    #[must_use]
    pub const fn new_biased_dna_matrix(matching: i8, mismatch: i8, ignoring: Option<u8>) -> Self {
        WeightMatrix::new(&DNA_PROFILE_MAP, matching, mismatch, ignoring).to_biased_matrix()
    }
}

impl<const S: usize> WeightMatrix<'_, u8, S> {
    /// Recomputes the matrix bias.
    const fn recompute_bias(&mut self) {
        let mut min = u8::MAX;
        let mut i = 0;
        while i < S {
            let mut j = 0;
            while j < S {
                if self.weights[i][j] < min {
                    min = self.weights[i][j];
                }
                j += 1;
            }
            i += 1;
        }

        let new_bias = self.bias.saturating_sub(min);
        let decrement = self.bias - new_bias;

        let mut i = 0;
        while i < S {
            let mut j = 0;
            while j < S {
                self.weights[i][j] -= decrement;
                j += 1;
            }
            i += 1;
        }

        self.bias = new_bias;
    }

    /// Subsets the rows and columns of a weight matrix to only include the
    /// bytes in `subset`. The bias is recomputed.
    ///
    /// This can also be used to reorder the rows/columns in the matrix.
    #[must_use]
    pub const fn get_subset<'a, const L: usize>(&self, subset: &'a ByteIndexMap<L>) -> WeightMatrix<'a, u8, L> {
        let mut out = self.get_subset_unadjusted(subset);
        out.recompute_bias();
        out
    }
}

impl<'a, const S: usize> WeightMatrix<'a, i8, S> {
    /// Creates a new, signed [`WeightMatrix`] with a given alphabet represented
    /// by `mapping`, a fixed `matching` score and `mismatch` score, and an
    /// optionally ignored base. A pair of bases where either is the ignored
    /// base will always have a score of 0.
    ///
    /// If working with DNA, consider using [`new_dna_matrix`]. For more
    /// flexibility, use [`new_custom`].
    ///
    /// [`new_dna_matrix`]: WeightMatrix::new_dna_matrix
    /// [`new_custom`]: WeightMatrix::new_custom
    ///
    /// ## Panics
    ///
    /// Panics if an invalid byte was specified for the `ignoring` field.
    #[must_use]
    pub const fn new(mapping: &'a ByteIndexMap<S>, matching: i8, mismatch: i8, ignoring: Option<u8>) -> Self {
        let mut weights = [[0i8; S]; S];

        let skip_index = match ignoring {
            Some(ignoring) => {
                if mapping.in_byte_keys(ignoring) {
                    Some(mapping.to_index(ignoring))
                } else {
                    panic!("An invalid byte was specified for the ignoring field.")
                }
            }
            None => None,
        };

        let mut i = 0;
        while i < S {
            let mut j = 0;
            while j < S {
                if let Some(k) = skip_index
                    && (k == i || k == j)
                {
                    j += 1;
                    continue;
                }

                if i == j {
                    weights[i][j] = matching;
                } else {
                    weights[i][j] = mismatch;
                }

                j += 1;
            }
            i += 1;
        }

        WeightMatrix {
            weights,
            mapping,
            bias: 0,
        }
    }

    /// Creates a new, signed [`WeightMatrix`] with a given alphabet represented
    /// by `mapping` and a custom weight matrix (where the rows represent the
    /// reference residue, and the columns represent the query residue).
    ///
    /// For simpler weight matrices depending only on matches/mismatches,
    /// consider using [`new`].
    ///
    /// [`new`]: WeightMatrix::new
    #[must_use]
    pub const fn new_custom(mapping: &'a ByteIndexMap<S>, weights: [[i8; S]; S]) -> Self {
        WeightMatrix {
            weights,
            mapping,
            bias: 0,
        }
    }

    /// Creates a new, signed [`WeightMatrix`] from a closure.
    ///
    /// The closure should accept the reference residue as the first argument
    /// and the query residue as the second.
    pub fn new_from_fn<F>(mapping: &'a ByteIndexMap<S>, weight_fn: F) -> Self
    where
        F: Fn(u8, u8) -> i8, {
        let mut out = WeightMatrix::new(mapping, 0, 0, None);
        for ref_residue in mapping.byte_keys {
            for query_residue in mapping.byte_keys {
                let score = weight_fn(ref_residue, query_residue);
                *out.get_weight_mut(ref_residue, query_residue) = score;
            }
        }
        out
    }

    /// Computes the matrix bias, which should be lesser of 0 and the smallest
    /// non-positive number.
    #[must_use]
    const fn get_bias(&self) -> i8 {
        let mut min = 0;
        let mut i = 0;
        while i < S {
            let mut j = 0;
            while j < S {
                if self.weights[i][j] < min {
                    min = self.weights[i][j];
                }
                j += 1;
            }
            i += 1;
        }
        min
    }

    /// Converts the signed [`WeightMatrix`] to an unsigned, biased matrix.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    #[must_use]
    pub const fn to_biased_matrix(&self) -> WeightMatrix<'a, u8, S> {
        let bias = self.get_bias();
        let mut weights = [[0u8; S]; S];
        let mapping = self.mapping;

        let mut i = 0;
        while i < S {
            let mut j = 0;
            while j < S {
                // This quantity must be non-negative.
                weights[i][j] = (self.weights[i][j] as i16 - bias as i16) as u8;
                j += 1;
            }
            i += 1;
        }

        // We can provide the bias as the unsigned version.
        let bias = bias.unsigned_abs();

        WeightMatrix { weights, mapping, bias }
    }

    /// Subsets the rows and columns of a weight matrix to only include the
    /// bytes in `subset`.
    ///
    /// This can also be used to reorder the rows/columns in the matrix.
    #[must_use]
    pub const fn get_subset<'b, const L: usize>(&self, subset: &'b ByteIndexMap<L>) -> WeightMatrix<'b, i8, L> {
        self.get_subset_unadjusted(subset)
    }
}

impl WeightMatrix<'_, i8, 5> {
    /// Creates a new signed [`WeightMatrix`] with a fixed `matching` score,
    /// `mismatch` score, and optionally ignoring a base. A pair of bases where
    /// either is the ignored base will always have a score of 0.
    #[must_use]
    pub const fn new_dna_matrix(matching: i8, mismatch: i8, ignoring: Option<u8>) -> Self {
        WeightMatrix::new(&DNA_PROFILE_MAP, matching, mismatch, ignoring)
    }
}

impl<T: AnyInt, const S: usize> Display for WeightMatrix<'_, T, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut col_widths = self.weights.iter().fold([0; S], |acc, row| {
            let widths = row.map(|val| val.to_string().len());
            elem_max(&acc, &widths)
        });
        for width in col_widths.iter_mut().skip(1) {
            *width += 1;
        }

        let residues = self.mapping.byte_keys();

        write!(f, "   ")?;
        for (residue, width) in residues.iter().zip(&col_widths) {
            write!(f, "{residue:>width$}", residue = *residue as char)?;
        }
        writeln!(f)?;

        for (row, residue) in self.weights.iter().zip(residues) {
            write!(f, "{residue}  ", residue = *residue as char)?;
            for (val, width) in row.iter().zip(&col_widths) {
                write!(f, "{val:>width$}")?;
            }
            writeln!(f)?;
        }
        if self.bias != T::ZERO {
            write!(f, "Bias: {}", self.bias)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    static AA_MATS: [&WeightMatrix<'static, i8, 25>; 22] = [
        &BLOSUM_30,
        &BLOSUM_35,
        &BLOSUM_40,
        &BLOSUM_45,
        &BLOSUM_50,
        &BLOSUM_55,
        &BLOSUM_60,
        &BLOSUM_62,
        &BLOSUM_65,
        &BLOSUM_70,
        &BLOSUM_75,
        &BLOSUM_80,
        &BLOSUM_85,
        &BLOSUM_90,
        &BLOSUM_95,
        &BLOSUM_100,
        &PAM_30,
        &PAM_40,
        &PAM_70,
        &PAM_120,
        &PAM_200,
        &PAM_250,
    ];

    #[test]
    fn create_simple() {
        static RESIDUE_MAP: ByteIndexMap<2> = ByteIndexMap::new(*b"AB", b'A');

        let result1 = WeightMatrix {
            weights: [[1, 0], [0, 1]],
            mapping: &RESIDUE_MAP,
            bias:    0,
        };

        let result2 = WeightMatrix::new(&RESIDUE_MAP, 1, 0, None);

        let result3 = WeightMatrix::new_custom(&RESIDUE_MAP, [[1, 0], [0, 1]]);

        assert_eq!(result1, result2);
        assert_eq!(result1, result3);
    }

    #[test]
    #[should_panic(expected = "An invalid byte was specified for the ignoring field.")]
    fn test_invalid_char() {
        let _ = WeightMatrix::new_dna_matrix(1, 0, Some(b'U'));
    }

    #[allow(non_snake_case)]
    #[test]
    fn create_IRMA_matrix() {
        let result1 = WeightMatrix {
            weights: [
                [2, -5, -5, -5, 0],
                [-5, 2, -5, -5, 0],
                [-5, -5, 2, -5, 0],
                [-5, -5, -5, 2, 0],
                [0, 0, 0, 0, 0],
            ],
            mapping: &DNA_PROFILE_MAP,
            bias:    0,
        };

        let result2 = WeightMatrix::new(&DNA_PROFILE_MAP, 2, -5, Some(b'N'));

        let result3 = WeightMatrix::new_custom(
            &DNA_PROFILE_MAP,
            [
                [2, -5, -5, -5, 0],
                [-5, 2, -5, -5, 0],
                [-5, -5, 2, -5, 0],
                [-5, -5, -5, 2, 0],
                [0, 0, 0, 0, 0],
            ],
        );

        assert_eq!(result1, result2);
        assert_eq!(result1, result3);
    }

    #[test]
    fn test_blosum_matrices() {
        assert!(aa_mat_from_name("BLOSUM30").is_some());
        assert!(aa_mat_from_name("blosum 62").is_some());
        assert!(aa_mat_from_name("BLOSUM_80").is_some());
        assert!(aa_mat_from_name("blosum-100").is_some());
        assert!(aa_mat_from_name("BLOSUM 45").is_some());

        assert!(aa_mat_from_name("BLOSUM29").is_none());
        assert!(aa_mat_from_name("BLOSUM101").is_none());
        assert!(aa_mat_from_name("BLOSUM").is_none());
    }

    #[test]
    fn test_pam_matrices() {
        assert!(aa_mat_from_name("PAM30").is_some());
        assert!(aa_mat_from_name("pam-70").is_some());
        assert!(aa_mat_from_name("PAM_120").is_some());
        assert!(aa_mat_from_name("PAM 200").is_some());

        assert!(aa_mat_from_name("PAM39").is_none());
        assert!(aa_mat_from_name("PAM2500").is_none());
        assert!(aa_mat_from_name("PAM").is_none());
    }

    #[test]
    fn test_invalid_names() {
        assert!(aa_mat_from_name("INVALID30").is_none());
        assert!(aa_mat_from_name("UNKNOWN").is_none());
        assert!(aa_mat_from_name("PAM_").is_none());
        assert!(aa_mat_from_name("").is_none());
    }

    #[test]
    fn test_subset_unbiased() {
        const OTHER_MAP: ByteIndexMap<21> = ByteIndexMap::new_ignoring_case(*b"DIEVRSAWGMCYTXNPFLHKQ", b'X');

        for mat in AA_MATS {
            let other_mat = mat.get_subset(&OTHER_MAP);
            for aa1 in u8::MIN..=u8::MAX {
                for aa2 in u8::MIN..=u8::MAX {
                    let other_aa1 = if OTHER_MAP.byte_keys().contains(&aa1.to_ascii_uppercase()) {
                        aa1
                    } else {
                        b'X'
                    };
                    let other_aa2 = if OTHER_MAP.byte_keys().contains(&aa2.to_ascii_uppercase()) {
                        aa2
                    } else {
                        b'X'
                    };

                    let score = mat.get_weight(other_aa1, other_aa2);
                    let other_score = other_mat.get_weight(aa1, aa2);
                    assert_eq!(score, other_score, "{aa1} vs {aa2}");
                }
            }
        }
    }

    #[test]
    #[allow(clippy::cast_possible_wrap)]
    fn test_subset_biased() {
        const OTHER_MAP: ByteIndexMap<4> = ByteIndexMap::new_ignoring_case(*b"AWGX", b'X');

        for mat in AA_MATS.map(WeightMatrix::to_biased_matrix) {
            let other_mat = mat.get_subset(&OTHER_MAP);
            for aa1 in u8::MIN..=u8::MAX {
                for aa2 in u8::MIN..=u8::MAX {
                    let other_aa1 = if OTHER_MAP.byte_keys().contains(&aa1.to_ascii_uppercase()) {
                        aa1
                    } else {
                        b'X'
                    };
                    let other_aa2 = if OTHER_MAP.byte_keys().contains(&aa2.to_ascii_uppercase()) {
                        aa2
                    } else {
                        b'X'
                    };

                    let score =
                        i8::try_from(mat.get_weight(other_aa1, other_aa2)).unwrap() - i8::try_from(mat.bias).unwrap();
                    let other_score =
                        i8::try_from(other_mat.get_weight(aa1, aa2)).unwrap() - i8::try_from(other_mat.bias).unwrap();
                    assert_eq!(score, other_score, "{aa1} vs {aa2}");
                }
            }
        }
    }

    #[test]
    #[allow(clippy::cast_possible_wrap)]
    fn test_subset_correct_bias() {
        for byte_keys in generate_subsets::<u8, 14, 4>(b"ADFHIKMQTVY*BZ") {
            let mapping =
                ByteIndexMap::new_ignoring_case([byte_keys[0], byte_keys[1], byte_keys[2], byte_keys[3], b'X'], b'X');
            for mat in AA_MATS.map(WeightMatrix::to_biased_matrix) {
                let other_mat = mat.get_subset(&mapping);
                let min_weight = *other_mat.weights.map(|row| *row.iter().min().unwrap()).iter().min().unwrap();
                assert_eq!(min_weight, 0);
            }
        }
    }

    fn generate_subsets<T: Copy + Default, const S: usize, const L: usize>(arr: &[T; S]) -> Vec<[T; L]> {
        fn generate_subsets_inner<T: Copy, const S: usize, const L: usize>(
            arr: &[T; S], start: usize, current: [T; L], current_len: usize, out: &mut Vec<[T; L]>,
        ) {
            if current_len == L {
                out.push(current);
                return;
            }

            for i in start..S {
                let mut new_current = current;
                new_current[current_len] = arr[i];
                generate_subsets_inner::<T, S, L>(arr, i + 1, new_current, current_len + 1, out);
            }
        }

        let mut out = Vec::new();
        generate_subsets_inner::<T, S, L>(arr, 0, [T::default(); L], 0, &mut out);
        out
    }
}
