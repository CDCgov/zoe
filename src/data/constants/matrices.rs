use crate::{
    data::mappings::{ByteIndexMap, DNA_PROFILE_MAP},
    math::AnyInt,
};

// Physiochemical distance matrix using the euclidean distance between all amino acid factors.
pub(crate) static PHYSIOCHEMICAL_FACTORS: [[Option<f32>; 256]; 256] = {
    const AA: [u8; 43] = [
        b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R', b'S', b'T', b'V', b'W',
        b'Y', b'X', b'a', b'c', b'd', b'e', b'f', b'g', b'h', b'i', b'k', b'l', b'm', b'n', b'p', b'q', b'r', b's', b't',
        b'v', b'w', b'y', b'x', b'-',
    ];

    // # Citation
    // For factor analysis used by the function:
    //
    // > Atchley et al. 2008. "Solving the protein sequence metric problem." Proc
    // > Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Published 2005 Apr 25.
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
/// weight matrix, then call [`into_biased_matrix`]. This will shift all entries
/// in the matrix to become nonnegative, and then store the resulting shift. One
/// can also use [`new_biased_dna_matrix`].
///
/// [`new`]: WeightMatrix::new
/// [`new_dna_matrix`]: WeightMatrix::new_dna_matrix
/// [`new_biased_dna_matrix`]: WeightMatrix::new_biased_dna_matrix
/// [`into_biased_matrix`]: WeightMatrix::into_biased_matrix
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub struct WeightMatrix<T: AnyInt, const S: usize> {
    pub weights:     [[T; S]; S],
    pub mapping:     &'static ByteIndexMap<S>,
    pub(crate) bias: T,
}

impl<T: AnyInt, const S: usize> WeightMatrix<T, S> {
    /// For a given `reference_base` and `query_base`, retrieves the weight
    /// stored in the matrix.
    #[inline]
    #[must_use]
    pub const fn get_weight(&self, reference_base: u8, query_base: u8) -> T {
        self.weights[self.mapping.to_index(reference_base)][self.mapping.to_index(query_base)]
    }
}

impl WeightMatrix<u8, 5> {
    /// Creates a new [`WeightMatrix`] with a fixed `matching` score,
    /// `mismatch` score, and optionally ignoring a base. A pair of bases where
    /// either is the ignored base will always have a score of 0.
    #[inline]
    #[must_use]
    pub const fn new_biased_dna_matrix(matching: i8, mismatch: i8, ignoring: Option<u8>) -> Self {
        WeightMatrix::new(&DNA_PROFILE_MAP, matching, mismatch, ignoring).into_biased_matrix()
    }
}

impl<const S: usize> WeightMatrix<i8, S> {
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
    #[must_use]
    pub const fn new(mapping: &'static ByteIndexMap<S>, matching: i8, mismatch: i8, ignoring: Option<u8>) -> Self {
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
    pub const fn new_custom(mapping: &'static ByteIndexMap<S>, weights: [[i8; S]; S]) -> Self {
        WeightMatrix {
            weights,
            mapping,
            bias: 0,
        }
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
    pub const fn into_biased_matrix(self) -> WeightMatrix<u8, S> {
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
}

impl WeightMatrix<i8, 5> {
    /// Creates a new signed [`WeightMatrix`] with a fixed `matching` score,
    /// `mismatch` score, and optionally ignoring a base. A pair of bases where
    /// either is the ignored base will always have a score of 0.
    #[must_use]
    pub const fn new_dna_matrix(matching: i8, mismatch: i8, ignoring: Option<u8>) -> Self {
        WeightMatrix::new(&DNA_PROFILE_MAP, matching, mismatch, ignoring)
    }
}

#[cfg(test)]
mod test {
    use super::*;

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
}
