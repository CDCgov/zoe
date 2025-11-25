//! A specification struct for generating arbitrary arrays whose elements
//! conform to given specifications.

use crate::data::arbitrary::ArbitrarySpecs;
use arbitrary::{Result, Unstructured};

/// Specifications for generating an arbitrary array.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct ArraySpecs<S, const N: usize> {
    /// The specifications for generating each element of the array.
    pub element_specs: S,
}

impl<'a, S, const N: usize> ArbitrarySpecs<'a> for ArraySpecs<S, N>
where
    S: ArbitrarySpecs<'a>,
    S::Output: Default,
{
    type Output = [S::Output; N];

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let mut arr = std::array::from_fn(|_| Default::default());
        for val in &mut arr {
            *val = self.element_specs.make_arbitrary(u)?;
        }
        Ok(arr)
    }
}
