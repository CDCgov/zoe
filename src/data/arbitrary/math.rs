//! Specification structs for generating arbitrary numbers with additional
//! guarantees.

use crate::{data::arbitrary::ArbitrarySpecs, math::Float};
use arbitrary::{Arbitrary, Result, Unstructured};
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

    /// Whether to include `inf` as a possible output.
    pub include_infinite: bool,

    /// Assuming random bytes, the percent change to generate infinity instead
    /// of a random float.
    ///
    /// The fuzzer is not always good at generating multiple infinite values,
    /// even when the code being fuzzed explicitly has conditions involving
    /// whether the floats are infinite. To resolve this, a 0-100 percent change
    /// of generating infinity can be set.
    ///
    /// Any value above 100 is set to 100. This overrides `include_infinite` and
    /// ignores `max`.
    pub infinite_perc: Option<u8>,

    /// The minimum value to allow for the float.
    ///
    /// Any values not comparable to this value are permitted (e.g., `NaN`). The
    /// values specified here should not be `NaN`.
    pub min: T,

    /// The maximum value to allow for the float.
    ///
    /// Any values not comparable to this value are permitted (e.g., `NaN`). The
    /// values specified here should not be `NaN`.
    pub max: T,

    /// The desired floating point type.
    pub float: PhantomData<T>,
}

impl<T: Float> Default for FloatSpecs<T> {
    fn default() -> Self {
        Self {
            include_nan:      true,
            include_infinite: true,
            infinite_perc:    None,
            min:              T::NEG_INFINITY,
            max:              T::INFINITY,
            float:            PhantomData,
        }
    }
}

impl<'a, T> ArbitrarySpecs<'a> for FloatSpecs<T>
where
    T: Float + Arbitrary<'a>,
{
    type Output = T;

    /// Generates an arbitrary floating point number conforming to the given
    /// specifications.
    ///
    /// ## Errors
    ///
    /// Any errors from the underlying [`arbitrary`] calls are propagated.
    ///
    /// ## Panics
    ///
    /// This panics if `min` is not less than `max`.
    ///
    /// [`arbitrary`]: arbitrary::Arbitrary::arbitrary
    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        assert!(self.min < self.max, "The min field must be less than the max field");

        if let Some(infinite_perc) = self.infinite_perc
            && u.int_in_range(0..=99)? < infinite_perc
        {
            return Ok(T::INFINITY);
        }

        let mut float = T::arbitrary(u)?;

        if !self.include_nan && float.is_nan() {
            float = T::ZERO;
        }

        if !self.include_infinite && float.is_infinite() {
            float = T::ZERO;
        }

        match (self.min > T::NEG_INFINITY, self.max < T::INFINITY) {
            (true, true) => {
                if float < self.min || float > self.max {
                    float = (self.min + float.rem_euclid(self.max - self.min)).min(self.max).max(self.min);
                }
            }
            (true, false) => {
                if float < self.min {
                    float = (float + self.min).max(self.min);
                }
            }
            (false, true) => {
                if float > self.max {
                    float = (float - self.max).min(self.max);
                }
            }
            (false, false) => {}
        }

        Ok(float)
    }
}
