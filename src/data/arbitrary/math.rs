//! Specification structs for generating arbitrary numbers with additional
//! guarantees.

use crate::{data::arbitrary::ArbitrarySpecs, math::Float};
use arbitrary::{Arbitrary, Unstructured};
use std::marker::PhantomData;

/// Specifications for generating an arbitrary floating point number.
///
/// ## Parameters
///
/// `T` is the type of the float (`f32` or `f64`).
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct FloatSpecs<T> {
    /// Whether to include `NaN` as a possible output.
    pub include_nan: bool,

    /// Whether to restrict the output to be nonnegative.
    pub nonnegative: bool,

    /// The desired floating point type.
    pub float: PhantomData<T>,
}

impl<T> Default for FloatSpecs<T> {
    fn default() -> Self {
        Self {
            include_nan: true,
            nonnegative: false,
            float:       PhantomData,
        }
    }
}

impl<'a, T> ArbitrarySpecs<'a> for FloatSpecs<T>
where
    T: Float + Arbitrary<'a>,
{
    type Output = T;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        let mut float = T::arbitrary(u)?;

        if !self.include_nan && float.is_nan() {
            float = T::ZERO;
        }

        if self.nonnegative {
            float = float.abs();
        }

        Ok(float)
    }
}
