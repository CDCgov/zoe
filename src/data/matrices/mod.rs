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

pub(crate) mod aa;
mod parse;
mod weight_matrix;
mod weight_row;

pub use aa::*;
pub use parse::*;
pub use weight_matrix::*;

#[cfg(test)]
mod test;

/// Physiochemical distance matrix using the Euclidean distance between all
/// amino acid factors.
pub(crate) static PHYSIOCHEMICAL_FACTORS: [[Option<f32>; 256]; 256] = {
    const AA: [u8; 43] = *b"ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx-";

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
