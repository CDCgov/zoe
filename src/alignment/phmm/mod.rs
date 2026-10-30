//! Structs and algorithms for profile Hidden Markov Models (pHMMs).
//!
//! ## Alignment Modes
//!
//! Traditionally, a pHMM performs global alignment of a query against a model.
//! However, by adding [`modules`] to the beginning and end of the pHMM, it can
//! allow the pHMM to skip bases from the beginning/end of the query, and/or
//! enter the pHMM at any match state.
//!
//! In *Zoe*, [`SamHmmParser`] will by default return a [`GlobalPhmm`]. Then,
//! the methods [`into_local_phmm`], [`into_semilocal_phmm`], and
//! [`into_domain_phmm`] can be called to modify the pHMM for other alignment
//! modes. Some choice of parameters may be required, which are specified with
//! [`LocalConfig`], [`SemiLocalConfig`], and [`DomainConfig`] respectively.
//!
//! A summary of the alignment modes is below:
//!
//! - Global alignment ([`GlobalPhmm`]): aligns a full sequence to the entire
//!   model
//! - Local alignment ([`LocalPhmm`]): aligns a subsequence to a portion of the
//!   model
//! - Semilocal alignment ([`SemiLocalPhmm`]): aligns a full sequence to a
//!   portion of the model
//! - Domain alignment ([`DomainPhmm`]): aligns a subsequence to the full model
//!
//! ## Alignment Algorithms
//!
//! *Zoe* currently supports two algorithms:
//!
//! - *Viterbi*: The viterbi algorithm aligns a sequence $x$ against the model,
//!   returning the best alignment $\pi^\star$. It also returns the joint
//!   probability of that sequence and alignment occurring, $P(x, \pi^\star)$.
//!   Specifically, *Zoe* outputs this probability as $-\ln(P(x, \pi^\star))$.
//!   See [`GlobalPhmm::viterbi`].
//! - *Forward*: The forward algorithm computes the probability of observing a
//!   sequence $x$ from the model, $P(x)$. This considers all alignments the
//!   sequence could form. *Zoe* outputs this probability as $-\ln(P(x))$. See
//!   [`GlobalPhmm::forward`]. This is currently only supported for
//!   [`GlobalPhmm`].
//!
//! ## Log Space Parameters
//!
//! Internally, *Zoe* stores pHMM parameters in log space with $-\ln(\cdot)$.
//! This allows for algorithms to be performed in a numerically stable manner,
//! while also avoiding repeated logarithm computations.
//!
//! When the input model is probabilistic (so that all parameters are in
//! `0..=1`), this will result in stored parameters in `0..=INFINITY`. *Zoe*
//! does not currently provide functionality for parameters outside of this
//! range ([`SamHmmParser`] will return an error). If this is required
//! functionality, please raise an issue (supporting it will require more
//! careful handling of NaNs and negative infinity values).
//!
//! [`SamHmmParser`]: sam_parser::SamHmmParser
//! [`into_local_phmm`]: GlobalPhmm::into_local_phmm
//! [`into_semilocal_phmm`]: GlobalPhmm::into_semilocal_phmm
//! [`into_domain_phmm`]: GlobalPhmm::into_domain_phmm

use crate::math::Float;
use std::{
    fmt::Display,
    ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign},
};

mod forward;
pub mod indexing;
mod models;
pub mod modules;
mod normalize;
pub mod sam_parser;
pub mod state;
pub mod traverse;
mod viterbi;

#[cfg(feature = "rand")]
pub mod sampling;

#[cfg(feature = "alignment-diagnostics")]
pub mod score_from_path;

pub use models::*;
pub use viterbi::*;

/// A trait for numeric types compatible with pHMMs.
///
/// These numeric types are used for performing pHMM calculations in negative
/// log space. Under normal circumstances, the parameters should be
/// non-negative.
pub trait PhmmNumber:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + PartialOrd
    + Default
    + Display {
    /// Infinity, the negative log space score corresponding to probability
    /// zero.
    const INFINITY: Self;
    /// Zero, the negative log space score corresponding to probability one.
    const ZERO: Self;
    /// The natural logarithm of 2.
    const LN_2: Self;

    /// Converts a probability (a floating point value) into negative log space.
    fn from_prob<T: Float>(prob: T) -> Self;

    /// Converts a negative log space score back to a probability.
    fn to_prob<T: Float>(self) -> T;

    /// Returns the negative log space score as a float.
    fn to_float<T: Float>(self) -> T;

    /// Returns the minimum of two negative log space scores.
    #[must_use]
    fn min(self, other: Self) -> Self;

    /// Returns the maximum of two negative log space scores.
    #[must_use]
    fn max(self, other: Self) -> Self;

    /// Multiplies the parameter by a `usize`, used as a more efficient
    /// alternative to repeated addition.
    ///
    /// This must return [`Self::ZERO`] when `times` is 0, even in the presense
    /// of NaNs.
    #[must_use]
    fn mul_usize(self, times: usize) -> Self;

    /// Returns the exponential of the negative log space score, $e^s$.
    #[must_use]
    fn exp(self) -> Self;

    /// Returns `ln(1+s)` for the negative log space score, in a manner more
    /// accurate than if the operations were performed separately.
    #[must_use]
    fn ln_1p(self) -> Self;

    /// Returns whether the parameter corresponds to an "impossible" probability
    /// of 0.
    #[inline]
    #[must_use]
    fn is_impossible(self) -> bool {
        self == Self::INFINITY
    }

    /// Returns whether the parameter corresponds to a "possible" probability
    /// greater than 0.
    #[inline]
    #[must_use]
    fn is_possible(self) -> bool {
        self < Self::INFINITY
    }

    /// Returns the parameter clamped to a probability in `0..=1`.
    #[inline]
    #[must_use]
    fn clamp_valid_prob(self) -> Self {
        // Negative probabilities aren't possible to represent in negative log
        // space, so the only concern is probabilities greater than 1, which is
        // less than 0 in negative log space
        self.max(Self::ZERO)
    }
}

impl PhmmNumber for f32 {
    const INFINITY: Self = f32::INFINITY;
    const ZERO: Self = 0.0;
    const LN_2: Self = std::f32::consts::LN_2;

    #[inline]
    fn from_prob<T: Float>(prob: T) -> Self {
        // Increase precision by converting to f32 last
        let param = (-prob.ln()).cast_as::<f32>();
        if param.is_nan() { Self::INFINITY } else { param }
    }

    #[inline]
    fn to_prob<T: Float>(self) -> T {
        // Increase precision by converting from f32 first
        (-T::cast_from(self)).exp()
    }

    #[inline]
    fn to_float<T: Float>(self) -> T {
        T::cast_from(self)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        self.min(other)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        self.max(other)
    }

    #[inline]
    #[allow(clippy::cast_precision_loss)]
    fn mul_usize(self, times: usize) -> Self {
        if times == 0 { 0.0 } else { times as f32 * self }
    }

    #[inline]
    fn exp(self) -> Self {
        self.exp()
    }

    #[inline]
    fn ln_1p(self) -> Self {
        self.ln_1p()
    }
}

impl PhmmNumber for f64 {
    const INFINITY: Self = f64::INFINITY;
    const ZERO: Self = 0.0;
    const LN_2: Self = std::f64::consts::LN_2;

    #[inline]
    fn from_prob<T: Float>(prob: T) -> Self {
        // Increase precision by converting to f64 first
        let param = -prob.cast_as::<f64>().ln();
        if param.is_nan() { Self::INFINITY } else { param }
    }

    #[inline]
    fn to_prob<T: Float>(self) -> T {
        // Increase precision by converting from f64 last
        T::cast_from((-self).exp())
    }

    #[inline]
    fn to_float<T: Float>(self) -> T {
        T::cast_from(self)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        self.min(other)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        self.max(other)
    }

    #[inline]
    #[allow(clippy::cast_precision_loss)]
    fn mul_usize(self, times: usize) -> Self {
        if times == 0 { 0.0 } else { times as f64 * self }
    }

    #[inline]
    fn exp(self) -> Self {
        self.exp()
    }

    #[inline]
    fn ln_1p(self) -> Self {
        self.ln_1p()
    }
}

/// Similar to [`relative_sum_probs_slice`] but accepting an array for zero-cost
/// handling of edge cases (e.g., `N` as 0, 1, or 2).
///
/// Otherwise, this delegates to [`relative_sum_probs_slice`].
///
/// ## Validity
///
/// Unpredictable results may occur if any value is NaN.
///
/// ## Panics
///
/// `N` must be at least 1.
#[inline]
#[must_use]
fn relative_sum_probs<T: PhmmNumber, const N: usize>(vals: [T; N]) -> (T, T) {
    match N {
        0 => panic!("relative_sum_probs requires N > 0"),
        1 => (T::ZERO, vals[0]),
        2 => relative_add_probs(vals[0], vals[1]),
        _ => relative_sum_probs_slice(&vals),
    }
}

/// Computes the [`PhmmNumber`] corresponding to the sum of the probabilities of
/// the passed [`PhmmNumber`] parameters, but with the minimum [`PhmmNumber`]
/// subtracted off from all parameters to reduce numerical instability.
/// Equivalently, this can be viewed as computing the [`PhmmNumber`]
/// corresponding to the sum of the probabilities, but where each probability
/// has been divided by the maximum probability.
///
/// The minimum [`PhmmNumber`] which was subtracted off is also returned as a
/// second value. To "undo" the normalization, add this parameter back to the
/// normalized sum (the first returned value). This is not done automatically
/// because users of this function can have better numerical stability by
/// rearranging the order of these operations.
///
/// ## Implementation Details
///
/// Let $v_{i}$ be the minimum parameter. This computes $-\ln\left(\sum_{j}
/// e^{-(v_{j}-v_{i})}\right)$. This is equivalent to $-\ln\left(1+\sum_{j\neq
/// i} e^{-(v_{j}-v_{i})}\right)$. By choosing to subtract the minimum $v_{i}$,
/// the exponents are minimized and hence numerical instability due to large
/// $e^{\star}$ values is reduced.
///
/// If $v_{i}=\infty$, then all $v_{j}=\infty$. This makes the above formulation
/// undefined. For consistency with the result as we take the limit as
/// $v_{i}\rightarrow\infty$, we return [`PhmmNumber::ZERO`] for the sum and
/// [`PhmmNumber::INFINITY`] for the minimum.
///
/// Let $N$ be the length of the slice. When using the results, note that the
/// first returned value is bounded between $-\ln(N)$ and $0$, whereas the
/// second returned value is on a similar magnitude as the parameters (it can be
/// anywhere from $0$ to $\infty$). Hence, adding the two directly can lose
/// significant precision from the first returned argument when $v_{i}$ is high.
/// See [`normalize`] for a more stable rearrangement.
///
/// ## Validity
///
/// Unpredictable results may occur if any value is NaN.
///
/// ## Panic
///
/// This will panic on an empty slice. For performance it should only be called
/// on slices of size 3 or more.
#[inline]
#[must_use]
fn relative_sum_probs_slice<T: PhmmNumber>(vals: &[T]) -> (T, T) {
    // Get argmin and min of vals
    let mut argmin = 0;
    let mut min = vals[0];
    for (i, val) in vals.iter().copied().enumerate().skip(1) {
        if val < min {
            argmin = i;
            min = val;
        }
    }

    // First: If any term is infinite probability, then the sum after dividing
    // will be 0. Second: See docs
    if min == -T::INFINITY || min == T::INFINITY {
        return (T::ZERO, min);
    }

    // Compute sum of exponents excluding the minimum
    let mut sum = T::ZERO;
    for val in &vals[..argmin] {
        sum += (min - *val).exp();
    }
    for val in &vals[argmin + 1..] {
        sum += (min - *val).exp();
    }

    // Return final result
    (-sum.ln_1p(), min)
}

/// Similar to [`relative_sum_probs_slice`] but for only two values.
///
/// This method relies on using the [`softplus`] function. Infinite values are
/// handled properly.
///
/// ## Validity
///
/// Unpredictable results may occur if any value is NaN.
fn relative_add_probs<T: PhmmNumber>(a: T, b: T) -> (T, T) {
    if a == b {
        // Short circuit for equal values, so that infinities are handled
        // properly instead of creating NaN
        (-T::LN_2, a)
    } else {
        let difference = a - b;

        // Dispatch on the sign of difference to ensure the exponent is
        // non-positive to prevent numerical instability
        if difference > T::ZERO {
            (-softplus(-difference), b)
        } else {
            (-softplus(difference), a)
        }
    }
}

/// Computes $\ln(1+e^a)$, also known as the softplus operation.
///
/// ## Limitations
///
/// This currently is implemented using `ln_1p` and `exp`. For greater
/// efficiency, it may be worth experimenting with piecewise approximations as
/// in
/// <https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf>,
/// or using a lookup table with interpolation. This also means smaller values
/// are more numerically stable, lest the `exp` operation overflow.
fn softplus<T: PhmmNumber>(a: T) -> T {
    a.exp().ln_1p()
}

/// Computes the [`PhmmNumber`] corresponding to the sum of the probabilities of
/// the passed [`PhmmNumber`] parameters.
///
/// This is equivalent to adding the returned values from
/// [`relative_sum_probs`]. Consider using [`relative_sum_probs`] directly if it
/// can help with numerical stability.
///
/// ## Validity
///
/// Unpredictable results may occur if any value is NaN.
///
/// ## Panics
///
/// `N` must be at least 1.
#[inline]
#[must_use]
fn sum_probs<T: PhmmNumber, const N: usize>(vals: [T; N]) -> T {
    let (rel_sum, min) = relative_sum_probs(vals);
    min + rel_sum
}
