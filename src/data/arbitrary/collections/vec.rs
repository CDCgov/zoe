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
    pub min_len: usize,

    /// The exact length of the [`Vec`] to generate.
    ///
    /// If set, this ignores the `min_len` and `max_len` fields.
    pub len: Option<usize>,

    /// The maximum length of the [`Vec`].
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
    /// The `len` field must be between `min_len` and `max_len`.
    ///
    /// [`arbitrary`]: arbitrary::Arbitrary::arbitrary
    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let vec = if let Some(len) = self.len {
            std::iter::repeat_with(|| self.element_specs.make_arbitrary(u))
                .take(len)
                .collect::<Result<Vec<_>>>()?
        } else {
            let start = std::iter::repeat_with(|| self.element_specs.make_arbitrary(u)).take(self.min_len);

            let mut out = start.collect::<Result<Vec<_>>>()?;

            let remaining_len = self.max_len - self.min_len;
            let remaining = self.element_specs.make_arbitrary_iter(u).take(remaining_len);

            remaining.process_results(|iter| {
                out.extend(iter);
            })?;

            out
        };

        Ok(vec)
    }
}
