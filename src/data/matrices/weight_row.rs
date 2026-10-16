//! Defines [`WeightRow`], representing a single row of a [`WeightMatrix`]
//! (allowing alignment algorithms to avoid unnecessary indexing).
//!
//! [`WeightMatrix`]: super::WeightMatrix

use crate::{data::ByteIndexMap, math::AnyInt};

/// A single row of a [`WeightMatrix`], containing the weights when comparing
/// against a fixed reference residue.
///
/// [`WeightMatrix`]: super::WeightMatrix
pub(crate) struct WeightRow<'a, T: AnyInt, const S: usize> {
    /// The row of weights
    pub(crate) weights: [T; S],

    /// The mapping used to index into `weights`
    pub(crate) mapping: &'a ByteIndexMap<S>,
}

impl<T: AnyInt, const S: usize> WeightRow<'_, T, S> {
    /// For a given `query_residue`, retrieves the weight stored in the row.
    #[inline]
    #[must_use]
    pub const fn get_weight(&self, query_residue: u8) -> T {
        self.weights[self.mapping.to_index(query_residue)]
    }
}
