#![allow(dead_code)]

use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Mul, Neg, Sub},
    simd::SimdElement,
    str::FromStr,
};

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
    + SimdElement
    + Sized
    + FromStr
    + std::iter::Sum<Self>
    + for<'a> std::iter::Sum<&'a Self> {
    const MIN: Self;
    const MAX: Self;
    const ZERO: Self;
    const ONE: Self;
    const INFINITY: Self;

    /// Generic absolute value for [`Float`]
    fn abs(self) -> Self;
    /// Generic minimum of 2 values for [`Float`]
    fn min(self, other: Self) -> Self;
    /// Generic log2 for [`Float`]
    fn log2(self) -> Self;
    /// Generic `is_nan` calculation for [`Float`]
    fn is_nan(self) -> bool;
    /// Generic `is_infinite` for [`Float`]
    fn is_infinite(self) -> bool;
    /// Generic natural logarithm for [`Float`]
    fn ln(self) -> Self;
    /// Generic exponential for [`Float`]
    fn exp(self) -> Self;

    /// Use a primitive cast to convert a usize to the [`Float`]
    fn usize_as_self(a: usize) -> Self;
}

macro_rules! impl_float {
    {$($ty:ty),* } => {
        $(
        impl Float for $ty {
            const MIN: Self = <$ty>::MIN;
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

pub(crate) trait NearlyEqual<T> {
    fn nearly_equal(self, b: Self, eps: T) -> bool;
}

impl<T> NearlyEqual<T> for T
where
    T: Float,
{
    #[inline]
    /// Tests if two floating points are approximately equal within an epsilon.
    /// Port courtesy of <https://floating-point-gui.de/errors/comparison/>
    fn nearly_equal(self, b: Self, eps: T) -> bool {
        let a = self;
        let abs_a = a.abs();
        let abs_b = b.abs();
        let diff = (a - b).abs();
        let zero = Self::ZERO;

        if a == b {
            // shortcut, handles infinities
            true
        } else if a == zero || b == zero || (abs_a + abs_b < Self::MIN) {
            // a or b is zero or both are extremely close to it
            // relative error is less meaningful here
            diff < eps * Self::MIN
        } else {
            // use relative error
            diff / (abs_a + abs_b).min(Self::MAX) < eps
        }
    }
}

impl<T> NearlyEqual<T> for Option<T>
where
    T: Float,
{
    #[inline]
    fn nearly_equal(self, b: Self, eps: T) -> bool {
        match (self, b) {
            (Some(x), Some(y)) => x.nearly_equal(y, eps),
            (None, None) => true,
            _ => false,
        }
    }
}

/// Utlity trait for mapping floats to [`Option`].
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

impl MapFloat for f64 {}
impl MapFloat for f32 {}

#[macro_export]
macro_rules! assert_fp_eq {
    ($a:expr, $b:expr) => {
        assert_fp_eq!($a, $b, 1e-8); // Default epsilon
    };
    ($a:expr, $expected:expr, $epsilon:expr) => {
        let found = $a;
        assert!(
            $crate::math::NearlyEqual::nearly_equal(found, $expected, $epsilon),
            "assertion failed: `(found ≈ expected)`\n left:\t`{:?}`,\n right:\t`{:?}`,\n eps:\t`{}`",
            found,
            $expected,
            $epsilon
        );
    };
}
