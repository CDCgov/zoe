use crate::{
    data::{
        array_types::{elem_max, is_subset},
        mappings::{ByteIndexMap, DNA_PROFILE_MAP},
    },
    math::AnyInt,
};
use std::fmt::Display;

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
    /// This field is only nonzero if `T` is unsigned
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
            is_subset(self.mapping.byte_keys(), subset.byte_keys()),
            "The provided keys are not a subset of the original weight matrix!"
        );
        let mut weights = WeightMatrix {
            weights: [[T::ZERO; L]; L],
            mapping: subset,
            bias:    self.bias,
        };
        let mut i = 0;
        while i < L {
            let reference_residue = subset.byte_keys()[i];
            let mut j = 0;
            while j < L {
                let query_residue = subset.byte_keys()[j];
                *weights.get_weight_mut(reference_residue, query_residue) =
                    self.get_weight(reference_residue, query_residue);
                j += 1;
            }
            i += 1;
        }
        weights
    }
}

impl WeightMatrix<'_, u8, 5> {
    /// Creates a new [`WeightMatrix`] with a fixed `matching` score, `mismatch`
    /// score, and optionally ignoring a base.
    ///
    /// A pair of bases where either is the ignored base will always have a
    /// score of 0. Note that this method uses the [`DNA_PROFILE_MAP`] mapping,
    /// which includes `N`.
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
    /// optionally ignored residue.
    ///
    /// A pair of residues where either is the ignored residue will always have
    /// a score of 0.
    ///
    /// If working with DNA, consider using [`new_dna_matrix`]. For more
    /// flexibility, use [`new_custom`].
    ///
    /// ## Panics
    ///
    /// Panics if an invalid byte was specified for the `ignoring` field.
    ///
    /// [`new_dna_matrix`]: WeightMatrix::new_dna_matrix
    /// [`new_custom`]: WeightMatrix::new_custom
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
    /// by `mapping` and a custom weight matrix.
    ///
    /// The rows of `weights` represent the reference residue, and the columns
    /// represent the query residue. For simpler weight matrices depending only
    /// on matches/mismatches, consider using [`new`].
    ///
    /// ## Validity
    ///
    /// If a asymmetric weight matrix is provided, then the alignment profile
    /// must be built from the query sequence, and all alignment functions must
    /// then be provided the reference sequence. Otherwise, the scores for `A →
    /// B` and `B → A` will be erroneously swapped.
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
    ///
    /// ## Validity
    ///
    /// If a asymmetric weight matrix is provided, then the profile must be
    /// built from the query sequence, and all alignment functions must then be
    /// provided the reference sequence. Otherwise, the scores for `A → B` and
    /// `B → A` will be erroneously swapped.
    pub fn new_from_fn<F>(mapping: &'a ByteIndexMap<S>, weight_fn: F) -> Self
    where
        F: Fn(u8, u8) -> i8, {
        let mut out = WeightMatrix::new(mapping, 0, 0, None);
        for ref_residue in mapping.byte_keys() {
            for query_residue in mapping.byte_keys() {
                let score = weight_fn(*ref_residue, *query_residue);
                *out.get_weight_mut(*ref_residue, *query_residue) = score;
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

    /// Generates an unsigned, biased [`WeightMatrix`] from the signed matrix.
    #[must_use]
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
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
    /// `mismatch` score, and optionally ignoring a base.
    ///
    /// A pair of bases where either is the ignored base will always have a
    /// score of 0. Note that this method uses the [`DNA_PROFILE_MAP`] mapping,
    /// which includes `N`.
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
