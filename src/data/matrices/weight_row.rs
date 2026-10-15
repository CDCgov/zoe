use crate::{data::ByteIndexMap, math::AnyInt};

pub(crate) struct WeightRow<'a, T: AnyInt, const S: usize> {
    pub(crate) weights: [T; S],
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
