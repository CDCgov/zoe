use crate::math::Float;
use arbitrary::{Arbitrary, Unstructured};

use super::{impl_deref, impl_from};

/// A wrapper around a float such that the [`Arbitrary`] implementation is
/// always not a NaN value.
///
/// ## Parameters
///
/// - `T`: The type of the float (`f32` or `f64`)
#[derive(Debug, Clone, Copy)]
pub struct FloatNotNan<T>(pub T);

impl_deref! {FloatNotNan<T>, T, <T>}
impl_from! {FloatNotNan, f32, f64}

impl<'a, T: Float + Arbitrary<'a>> Arbitrary<'a> for FloatNotNan<T> {
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        let val: T = Arbitrary::arbitrary(u)?;
        let cleaned = if val.is_nan() { T::ZERO } else { val };
        Ok(FloatNotNan(cleaned))
    }
}

/// A wrapper around a float such that the [`Arbitrary`] implementation is
/// always nonnegative (and not a NaN).
///
/// ## Parameters
///
/// - `T`: The type of the float (`f32` or `f64`)
#[derive(Debug, Clone, Copy)]
pub struct NonnegativeFloat<T>(pub T);

impl_deref! {NonnegativeFloat<T>, T, <T>}
impl_from! {NonnegativeFloat, f32, f64}

impl<'a, T: Float + Arbitrary<'a>> Arbitrary<'a> for NonnegativeFloat<T> {
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        Ok(NonnegativeFloat(u.arbitrary::<FloatNotNan<T>>()?.0.abs()))
    }
}
