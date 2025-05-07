/// A trait unifying the methods for comparing two floating point numbers. This
/// is implemented for [`NearlyEqualStrategy`] (provided methods for comparing
/// floats) and `Fn(T, T) -> bool` closures (custom comparison methods).
///
/// <div class="warning">
///
/// **Warning**
///
/// This is an implementation detail and *should not* be called directly! Use
/// [`is_fp_equal`] or [`assert_fp_equal`] instead.
///
/// </div>
#[doc(hidden)]
pub trait NearlyEqualMethod<T> {
    fn nearly_equal_float(&self, a: T, b: T) -> bool;
}

/// An enum for provided floating point comparison methods.
///
/// <div class="warning">
///
/// **Warning**
///
/// This is an implementation detail and *should not* be called directly! Use
/// [`is_fp_equal`] or [`assert_fp_equal`] instead.
///
/// </div>
#[doc(hidden)]
pub enum NearlyEqualStrategy<T> {
    /// Compare the floating point values using relative error with a tolerance
    /// of `eps`
    Relative { eps: T },
    /// Compare the floating point values using absolute error with a tolerance
    /// of `eps`
    Absolute { eps: T },
}

impl<T: Float> NearlyEqualMethod<T> for NearlyEqualStrategy<T> {
    /// Tests if two floating points are approximately equal within an epsilon.
    /// Relative error port courtesy of
    /// <https://floating-point-gui.de/errors/comparison/>
    fn nearly_equal_float(&self, a: T, b: T) -> bool {
        match self {
            NearlyEqualStrategy::Relative { eps } => {
                let abs_a = a.abs();
                let abs_b = b.abs();
                let diff = (a - b).abs();

                if a == b {
                    // shortcut, handles infinities
                    true
                } else if a == T::ZERO || b == T::ZERO || (abs_a + abs_b < T::MIN_POSITIVE) {
                    // a or b is zero or both are extremely close to it
                    // relative error is less meaningful here
                    diff < *eps * T::MIN_POSITIVE
                } else {
                    // use relative error
                    diff / (abs_a + abs_b).min(T::MAX) < *eps
                }
            }
            NearlyEqualStrategy::Absolute { eps } => a == b || (a - b).abs() < *eps,
        }
    }
}

impl<T: Float, F: Fn(T, T) -> bool> NearlyEqualMethod<T> for F {
    fn nearly_equal_float(&self, a: T, b: T) -> bool {
        self(a, b)
    }
}

/// A trait for enabling equality comparisons involving floating point numbers
/// using a given comparison method. This trait can be implemented for
/// primitives and for more complex types containing floating point numbers. `T`
/// is the underlying floating point type (and the type used for the tolerance).
///
/// <div class="warning">
///
/// **Warning**
///
/// This is an implementation detail and *should not* be called directly! Use
/// [`is_fp_equal`] or [`assert_fp_equal`] instead.
///
/// </div>
#[doc(hidden)]
pub trait NearlyEqual<T> {
    /// Check whether two types are approximately equal using `method`. The
    /// first value in the tuple is whether they are approximately equal.
    ///
    /// The second is optionally a tuple of floats which are responsible for the
    /// inequality. For example, if two unequal `Vec<f32>` are being compared,
    /// then the second value will be `None` if the vecs are unequal due to
    /// being different lengths, or otherwise will be `Some(val1, val2)` where
    /// `val1` and `val2` are the first unequal floats.
    fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, method: &M) -> (bool, Option<(T, T)>);
}

/// Implement [`NearlyEqual`] for floating point primitives.
macro_rules! impl_float_nearly_equal {
    {$($ty:ty),* } => {
        $(
            impl NearlyEqual<$ty> for $ty {
                #[inline]
                #[allow(clippy::float_cmp)]
                fn nearly_equal<M: NearlyEqualMethod<$ty>>(&self, b: &Self, strategy: &M) -> (bool, Option<($ty, $ty)>) {
                    if strategy.nearly_equal_float(*self, *b) {
                        (true, None)
                    } else {
                        (false, Some((*self, *b)))
                    }
                }
            }
        )*
    }
}

impl_float_nearly_equal!(f32, f64);

impl<T, S: NearlyEqual<T>> NearlyEqual<T> for Option<S> {
    #[inline]
    fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, strategy: &M) -> (bool, Option<(T, T)>) {
        match (self, b) {
            (Some(x), Some(y)) => x.nearly_equal(y, strategy),
            (None, None) => (true, None),
            _ => (false, None),
        }
    }
}

impl<T: Copy, S: NearlyEqual<T>, const N: usize> NearlyEqual<T> for [S; N] {
    #[inline]
    fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, strategy: &M) -> (bool, Option<(T, T)>) {
        for (eq, vals) in self.iter().zip(b).map(|(x, y)| x.nearly_equal(y, strategy)) {
            if !eq {
                return (false, vals);
            }
        }
        (true, None)
    }
}

impl<T: Copy, S: NearlyEqual<T>> NearlyEqual<T> for Vec<S> {
    #[inline]
    fn nearly_equal<M: NearlyEqualMethod<T>>(&self, b: &Self, strategy: &M) -> (bool, Option<(T, T)>) {
        if self.len() == b.len() {
            for (eq, vals) in self.iter().zip(b).map(|(x, y)| x.nearly_equal(y, strategy)) {
                if !eq {
                    return (false, vals);
                }
            }
            (true, None)
        } else {
            (false, None)
        }
    }
}

/// Assert that two floating point values are approximately equal.
///
/// ## Comparison Options
///
/// The most basic usage of this macro is:
/// ```
/// # use zoe::assert_fp_eq;
/// assert_fp_eq!(3.0, 1.0 + 2.0);
/// ```
///
/// This uses a default tolerance of $\epsilon=10^{-8}$, and compares the floats
/// using relative error (see
/// <https://floating-point-gui.de/errors/comparison/>). To customize the
/// tolerance, specify a third argument:
/// ```
/// # use zoe::assert_fp_eq;
/// assert_fp_eq!(3.0, 2.99999, 1e-4);
/// ```
///
/// Additionally, the method for comparing the floats can be customized. By
/// default it is relative error (`@relative`), but it can also be changed to
/// absolute error (`@absolute`) or a custom closure (`@custom`).
/// ```
/// # use zoe::{assert_fp_eq, is_fp_eq};
/// assert_fp_eq!(@relative, 3.0, 2.99999, 1e-4);
/// assert_fp_eq!(@absolute, 1e-10, 2e-10, 1e-8);
///
/// let comparison_method = |x: f32, y: f32| {
///     let both_are_small = x.abs() < 1e-8 && y.abs() < 1e-8;
///     both_are_small || is_fp_eq!(@relative, x, y, 1e-8)
/// };
///
/// assert_fp_eq!(@custom, 1e-10, 2e-10, comparison_method);
/// assert_fp_eq!(@custom, 3.0, 1.0+2.0, comparison_method);
/// ```
///
/// ## Allowed Arguments
///
/// This macro allows for `f32` and `f64` to be compared. Additionally, assuming
/// that type `T` can be compared using this macro, so can `Option<T>`, `[T;
/// N]`, and `Vec<T>`. Additionally, the following Zoe structs can be compared:
/// * Profile Hidden Markov Models: [`TransitionParams`], [`EmissionParams`],
///   [`LayerParams`], and [`Phmm`]
///
/// [`TransitionParams`]: crate::alignment::phmm::TransitionParams
/// [`EmissionParams`]: crate::alignment::phmm::EmissionParams
/// [`LayerParams`]: crate::alignment::phmm::LayerParams
/// [`Phmm`]: crate::alignment::phmm::Phmm
#[macro_export]
macro_rules! assert_fp_eq {
    ($(@$method:tt,)? $a:expr, $b:expr) => {
        assert_fp_eq!($(@$method,)? $a, $b, 1e-8); // Default epsilon
    };
    ($(@relative,)? $a:expr, $b:expr, $epsilon:expr) => {
        let (eq, vals) = $crate::math::NearlyEqual::nearly_equal(
            &$a,
            &$b,
            &$crate::math::NearlyEqualStrategy::Relative { eps: $epsilon }
        );
        if !eq {
            if let Some((a, b)) = vals {
                panic!("assertion failed: `(found ≈ expected)`\n left:\t`{:?}`,\n right:\t`{:?}`,\n eps:\t`{}`,\n\n Caused by the comparison of:\n left:\t`{:?}`,\n right:\t`{:?}`", $a, $b, $epsilon, a, b)
            } else {
                panic!("assertion failed: `(found ≈ expected)`\n left:\t`{:?}`,\n right:\t`{:?}`,\n eps:\t`{}`", $a, $b, $epsilon)
            }
        }
    };
    (@absolute, $a:expr, $b:expr, $epsilon:expr) => {
        let (eq, vals) = $crate::math::NearlyEqual::nearly_equal(
            &$a,
            &$b,
            &$crate::math::NearlyEqualStrategy::Absolute { eps: $epsilon }
        );
        if !eq {
            if let Some((a, b)) = vals {
                panic!("assertion failed: `(found ≈ expected)`\n left:\t`{:?}`,\n right:\t`{:?}`,\n eps:\t`{}`,\n\n Caused by the comparison of:\n left:\t`{:?}`,\n right:\t`{:?}`", $a, $b, $epsilon, a, b)
            } else {
                panic!("assertion failed: `(found ≈ expected)`\n left:\t`{:?}`,\n right:\t`{:?}`,\n eps:\t`{}`", $a, $b, $epsilon)
            }
        }
    };
    (@custom, $a:expr, $b:expr, $closure:expr) => {
        let (eq, vals) = $crate::math::NearlyEqual::nearly_equal(&$a, &$b, &$closure);
        if !eq {
            if let Some((a, b)) = vals {
                panic!("assertion failed: `(found ≈ expected)`\n left:\t`{:?}`,\n right:\t`{:?}`,\n\n Caused by the comparison of:\n left:\t`{:?}`,\n right:\t`{:?}`", $a, $b, a, b)
            } else {
                panic!("assertion failed: `(found ≈ expected)`\n left:\t`{:?}`,\n right:\t`{:?}`", $a, $b)
            }
        }
    };
}

/// Check whether two floating point values are approximately equal.
///
/// Similar to [`assert_fp_eq`], but returns a boolean rather than performing an
/// assertion.
#[macro_export]
macro_rules! is_fp_eq {
    ($(@$method:tt,)? $a:expr, $b:expr) => {
        is_fp_eq!($(@$method,)? $a, $b, 1e-8); // Default epsilon
    };
    ($(@relative,)? $a:expr, $b:expr, $epsilon:expr) => {
        $crate::math::NearlyEqual::nearly_equal(
            &$a,
            &$b,
            &$crate::math::NearlyEqualStrategy::Relative { eps: $epsilon }
        ).0
    };
    (@absolute, $a:expr, $b:expr, $epsilon:expr) => {
        $crate::math::NearlyEqual::nearly_equal(
            &$a,
            &$b,
            &$crate::math::NearlyEqualStrategy::Absolute { eps: $epsilon }
        ).0
    };
    (@custom, $a:expr, $b:expr, $closure:expr) => {
        $crate::math::NearlyEqual::nearly_equal(&$a, &$b, &$closure).0
    };
}

use super::Float;
