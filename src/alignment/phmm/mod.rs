//! Structs and algorithms for profile Hidden Markov Models (pHMMs).
//!
//! <div class="warning note">
//!
//! **Note**
//!
//! You must enable the *dev-phmm* feature in your `Cargo.toml` to use these
//! functions. They are in active development and are not complete.
//!
//! </div>

use crate::{
    alignment::phmm::indexing::PhmmIndex,
    math::{CastAs, CastAsNumeric, CastFrom, CastFromNumeric, Float},
};
use std::{
    fmt::Display,
    ops::{Add, AddAssign, Mul},
};

mod errors;
mod models;
pub mod modules;
pub mod sam_parser;
mod state;
mod traits;
mod viterbi;

#[cfg(not(feature = "alignment-diagnostics"))]
mod editing;
#[cfg(feature = "alignment-diagnostics")]
pub mod editing;

#[cfg(not(feature = "alignment-diagnostics"))]
pub(crate) mod indexing;
#[cfg(feature = "alignment-diagnostics")]
pub mod indexing;

#[cfg(feature = "alignment-diagnostics")]
pub mod score_from_path;

#[cfg(feature = "alignment-diagnostics")]
pub mod visit_params;

pub use errors::*;
pub use models::*;

#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
pub(crate) use state::*;

/// A trait for numeric types compatible with pHMMs.
///
/// These numeric types are used for performing pHMM calculations in negative
/// log space.
pub trait PhmmNumber:
    Copy
    + Add<Output = Self>
    + Mul<Output = Self>
    + AddAssign
    + PartialOrd
    + CastAs
    + CastAsNumeric
    + CastFrom
    + CastFromNumeric
    + Default
    + Display {
    /// Infinity, the negative log space score corresponding to probability zero
    const INFINITY: Self;
    /// Zero, the negative log space score corresponding to probability one
    const ZERO: Self;

    /// Converts a probability (a floating point value) into negative log space
    fn from_prob<T: Float>(prob: T) -> Self;

    /// Converts a negative log space score back to a probability
    fn to_prob<T: Float>(self) -> T;

    /// Returns the negative log space score as a float
    fn to_float<T: Float>(self) -> T;

    /// Computes the minimum of two negative log space scores
    #[must_use]
    fn min(self, other: Self) -> Self;
}

impl PhmmNumber for f32 {
    const INFINITY: Self = f32::INFINITY;
    const ZERO: Self = 0.0;

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
}

impl PhmmNumber for f64 {
    const INFINITY: Self = f64::INFINITY;
    const ZERO: Self = 0.0;

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
}
