#![allow(dead_code)]

use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, DivAssign, Mul, Neg, Sub},
    simd::SimdElement,
    str::FromStr,
};

use crate::private::Sealed;

/// Takes square root using the Babylonian Method using 40 iterations.
// Relevant discussion: <https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi>
pub(crate) const fn sqrt_baby(n: f64) -> f64 {
    let mut i: f64 = 0.0;
    while i * i <= n {
        i += 0.1_f64;
    }

    let mut x1 = i;
    let mut x2 = 0.0;

    let mut j: usize = 0;
    while j < 40 {
        x2 = (n / x1).midpoint(x1);
        x1 = x2;
        j += 1;
    }
    x2
}

/// Trait for providing generic functionality over floating point numbers.
pub trait Float:
    Sub<Output = Self>
    + Add<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + DivAssign
    + Sized
    + Default
    + PartialEq
    + PartialOrd
    + Copy
    + Debug
    + Display
    + SimdElement
    + Sized
    + FromStr
    + std::iter::Sum<Self>
    + for<'a> std::iter::Sum<&'a Self>
    + Sealed {
    const MIN: Self;
    const MIN_POSITIVE: Self;
    const MAX: Self;
    const ZERO: Self;
    const ONE: Self;
    const INFINITY: Self;

    /// Generic absolute value for [`Float`]
    #[must_use]
    fn abs(self) -> Self;

    /// Generic minimum of 2 values for [`Float`]
    #[must_use]
    fn min(self, other: Self) -> Self;

    /// Generic log2 for [`Float`]
    #[must_use]
    fn log2(self) -> Self;

    /// Generic `is_nan` calculation for [`Float`]
    #[must_use]
    fn is_nan(self) -> bool;

    /// Generic `is_infinite` for [`Float`]
    #[must_use]
    fn is_infinite(self) -> bool;

    /// Generic natural logarithm for [`Float`]
    #[must_use]
    fn ln(self) -> Self;

    /// Generic exponential for [`Float`]
    #[must_use]
    fn exp(self) -> Self;

    /// Use a primitive cast to convert a usize to the [`Float`]
    fn usize_as_self(a: usize) -> Self;
}

/// Implement [`Float`] for multiple floating point primitive types
macro_rules! impl_float {
    {$($ty:ty),* } => {
        $(
        impl Float for $ty {
            const MIN: Self = <$ty>::MIN;
            const MIN_POSITIVE: Self = <$ty>::MIN_POSITIVE;
            const MAX: Self = <$ty>::MAX;
            const ZERO: Self = 0.0;
            const ONE: Self = 1.0;
            const INFINITY: Self = <$ty>::INFINITY;

            #[inline]
            fn abs(self) -> Self {
                self.abs()
            }

            #[inline]
            fn min(self, other: Self) -> Self {
                self.min(other)
            }

            #[inline]
            fn log2(self) -> Self {
                self.log2()
            }

            #[inline]
            fn is_nan(self) -> bool {
                self.is_nan()
            }

            #[inline]
            fn is_infinite(self) -> bool {
                self.is_infinite()
            }

            #[inline]
            fn ln(self) -> Self {
                self.ln()
            }

            #[inline]
            fn exp(self) -> Self {
                self.exp()
            }

            #[inline]
            #[allow(clippy::cast_precision_loss)]
            fn usize_as_self(a: usize) -> $ty {
                a as $ty
            }
        } )*

     }
}

impl_float!(f32, f64);

/// Utility trait for mapping floats to [`Option`].
pub(crate) trait MapFloat: Float {
    #[inline]
    fn into_option(self) -> Option<Self> {
        if self.is_infinite() || self.is_nan() {
            None
        } else {
            Some(self)
        }
    }
}

impl MapFloat for f32 {}
impl MapFloat for f64 {}
