//! A specification struct for generating an arbitrary [`Vec`] whose elements
//! conform to given specifications, and whose length may be bounded or fixed.

use crate::{data::arbitrary::ArbitrarySpecs, iter_utils::ProcessResultsExt};
use arbitrary::{Result, Unstructured};

/// Specifications for generating an arbitrary [`Vec`].
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct VecSpecs<S> {
    /// The specifications for generating each element of the [`Vec`].
    pub element_specs: S,

    /// The minimum length of the [`Vec`].
    ///
    /// This must be less than or equal to `max_len`.
    pub min_len: usize,

    /// The exact length of the [`Vec`] to generate.
    ///
    /// If set, this ignores the `min_len` and `max_len` fields.
    pub len: Option<usize>,

    /// The maximum length of the [`Vec`].
    ///
    /// This must be greater than or equal to `min_len`.
    pub max_len: usize,
}

impl<S> Default for VecSpecs<S>
where
    S: Default,
{
    fn default() -> Self {
        Self {
            element_specs: S::default(),
            min_len:       0,
            len:           None,
            max_len:       usize::MAX,
        }
    }
}

impl<'a, S> ArbitrarySpecs<'a> for VecSpecs<S>
where
    S: ArbitrarySpecs<'a>,
{
    type Output = Vec<S::Output>;

    /// Generates an arbitrary vector conforming to the given specifications.
    ///
    /// ## Errors
    ///
    /// Any errors from the underlying [`arbitrary`] calls are propagated.
    ///
    /// ## Panics
    ///
    /// `min_len` must be less than or equal to `max_len`.
    ///
    /// [`arbitrary`]: arbitrary::Arbitrary::arbitrary
    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let Some(len_range) = self.max_len.checked_sub(self.min_len) else {
            panic!(
                "The min_len field must be less than or equal to the max_len field for VecSpecs (found min_len={min_len} and max_len={max_len}",
                min_len = self.min_len,
                max_len = self.max_len
            );
        };

        let vec = if let Some(len) = self.len {
            std::iter::repeat_with(|| self.element_specs.make_arbitrary(u))
                .take(len)
                .collect::<Result<Vec<_>>>()?
        } else {
            let start = std::iter::repeat_with(|| self.element_specs.make_arbitrary(u)).take(self.min_len);

            let mut out = start.collect::<Result<Vec<_>>>()?;

            let remaining = self.element_specs.make_arbitrary_iter(u).take(len_range);

            remaining.process_results(|iter| {
                out.extend(iter);
            })?;

            out
        };

        Ok(vec)
    }
}
