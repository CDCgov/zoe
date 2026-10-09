//! Structs and algorithms for profile Hidden Markov Models (pHMMs).
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

use crate::math::Float;
use std::{
    fmt::Display,
    ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign},
};

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
