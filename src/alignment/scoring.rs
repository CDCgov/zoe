//! Code supporting [`ScoreParams`], a struct holding the parameters for
//! alignment scoring.

use crate::{
    alignment::{ProfileError, validate_profile_args},
    data::WeightMatrix,
    math::AnyInt,
};
use std::borrow::Cow;

/// The parameters for scoring an alignment.
///
/// This includes the parameters for scoring a substitution in [`WeightMatrix`]
/// as well as the [affine gap
/// penalties](crate::alignment::sw#affine-gap-penalties).
///
/// An instance of [`ScoreParams`] can be constructed at run-time or
/// compile-time in many ways:
///
/// - [`new_dna_params`]: Constructs parameters corresponding to a simple DNA
///   alignment with signed parameters. It only supports specifying a match and
///   mismatch score, rather than an arbitrary [`WeightMatrix`]. It also uses
///   the [`DNA_PROFILE_MAP`] alphabet (`ACGTN`).
/// - [`new`]: A more general constructor accepting an arbitrary
///   [`WeightMatrix`] (and hence is compatible with other alphabets). See
///   [here](WeightMatrix) for ways to construct a [`WeightMatrix`].
/// - [`new_borrowed`] and [`new_owned`]: Versions of [`new`] that are
///   const-compatible. The difference involves whether an owned or borrowed
///   [`WeightMatrix`] is stored. Consider using a match statement to handle the
///   returned result.
/// - [`new_biased_dna_params`]: Similar to [`new_dna_params`] but with a biased
///   weight matrix. This is only needed when performing unsigned Striped Smith
///   Waterman alignment.
///
/// ## Parameters
///
/// See [`WeightMatrix`](WeightMatrix#parameters). `'a` is also used to store
/// the lifetime of the matrix, if it is borrowed. `'static` can be used if both
/// the weight matrix and alphabet are available at compile time.
///
/// [`DNA_PROFILE_MAP`]: crate::data::DNA_PROFILE_MAP
/// [`new_dna_params`]: ScoreParams::new_dna_params
/// [`new`]: ScoreParams::new
/// [`new_borrowed`]: ScoreParams::new_borrowed
/// [`new_owned`]: ScoreParams::new_owned
/// [`new_biased_dna_params`]: ScoreParams::new_biased_dna_params
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct ScoreParams<'a, T, const S: usize>
where
    T: AnyInt, {
    pub(crate) matrix:     Cow<'a, WeightMatrix<'a, T, S>>,
    pub(crate) gap_open:   i8,
    pub(crate) gap_extend: i8,
}

impl<'a, T, const S: usize> ScoreParams<'a, T, S>
where
    T: AnyInt,
{
    /// Constructs a new set of alignment scoring parameters from a
    /// [`WeightMatrix`] and non-positive weight scores.
    ///
    /// If the matrix is known at compile-time, consider using [`new_borrowed`]
    /// or [`new_owned`] at compile-time as well. The result can be handled with
    /// a match statement.
    ///
    /// ## Errors
    ///
    /// See [`ProfileError`].
    ///
    /// [`new_borrowed`]: Self::new_borrowed
    /// [`new_owned`]: Self::new_owned
    pub fn new<W>(matrix: W, gap_open: i8, gap_extend: i8) -> Result<Self, ProfileError>
    where
        W: Into<Cow<'a, WeightMatrix<'a, T, S>>>, {
        validate_profile_args(gap_open, gap_extend)?;

        Ok(Self {
            matrix: matrix.into(),
            gap_open,
            gap_extend,
        })
    }

    /// Constructs a new set of alignment scoring parameters from a borrowed
    /// [`WeightMatrix`] and non-positive weight scores.
    ///
    /// Prefer using [`ScoreParams::new`] in run-time code. For use in
    /// compile-time code, this can be called and a `match` statement can handle
    /// the result.
    ///
    /// ## Errors
    ///
    /// See [`ProfileError`].
    pub const fn new_borrowed(
        matrix: &'a WeightMatrix<'a, T, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, ProfileError> {
        if let Err(e) = validate_profile_args(gap_open, gap_extend) {
            return Err(e);
        }

        Ok(Self {
            matrix: Cow::Borrowed(matrix),
            gap_open,
            gap_extend,
        })
    }

    /// Constructs a new set of alignment scoring parameters from an owned
    /// [`WeightMatrix`] and non-positive weight scores.
    ///
    /// Prefer using [`ScoreParams::new`] in run-time code. For use in
    /// compile-time code, this can be called and a `match` statement can handle
    /// the result.
    ///
    /// ## Errors
    ///
    /// See [`ProfileError`].
    pub const fn new_owned(matrix: WeightMatrix<'a, T, S>, gap_open: i8, gap_extend: i8) -> Result<Self, ProfileError> {
        if let Err(e) = validate_profile_args(gap_open, gap_extend) {
            return Err(e);
        }

        Ok(Self {
            matrix: Cow::Owned(matrix),
            gap_open,
            gap_extend,
        })
    }
}

impl ScoreParams<'_, i8, 5> {
    /// Creates a new [`ScoreParams`] using a signed [`WeightMatrix`] with a
    /// fixed `matching` score, `mismatch` score, and optionally ignoring a
    /// base.
    ///
    /// See [`WeightMatrix::new_dna_matrix`] for more details. For use in
    /// compile-time code, this can be called and a `match` statement can handle
    /// the result.
    ///
    /// ## Errors
    ///
    /// See [`ProfileError`].
    pub const fn new_dna_params(
        matching: i8, mismatch: i8, gap_open: i8, gap_extend: i8, ignoring: Option<u8>,
    ) -> Result<Self, ProfileError> {
        if let Err(e) = validate_profile_args(gap_open, gap_extend) {
            return Err(e);
        }

        let matrix = WeightMatrix::new_dna_matrix(matching, mismatch, ignoring);

        Ok(Self {
            matrix: Cow::Owned(matrix),
            gap_open,
            gap_extend,
        })
    }
}

impl ScoreParams<'_, u8, 5> {
    /// Creates a new [`ScoreParams`] using an unsigned (biased)
    /// [`WeightMatrix`] with a fixed `matching` score, `mismatch` score, and
    /// optionally ignoring a base.
    ///
    /// See [`WeightMatrix::new_biased_dna_matrix`] for more details. For use in
    /// compile-time code, this can be called and a `match` statement can handle
    /// the result.
    ///
    /// ## Errors
    ///
    /// See [`ProfileError`].
    pub const fn new_biased_dna_params(
        matching: i8, mismatch: i8, gap_open: i8, gap_extend: i8, ignoring: Option<u8>,
    ) -> Result<Self, ProfileError> {
        if let Err(e) = validate_profile_args(gap_open, gap_extend) {
            return Err(e);
        }

        let matrix = WeightMatrix::new_biased_dna_matrix(matching, mismatch, ignoring);

        Ok(Self {
            matrix: Cow::Owned(matrix),
            gap_open,
            gap_extend,
        })
    }
}

impl<'a, 'b, T, const S: usize> From<&'b ScoreParams<'a, T, S>> for ScoreParams<'b, T, S>
where
    T: AnyInt,
{
    fn from(value: &'b ScoreParams<'a, T, S>) -> Self {
        ScoreParams {
            matrix:     Cow::Borrowed(&value.matrix),
            gap_open:   value.gap_open,
            gap_extend: value.gap_extend,
        }
    }
}

/// A helper function for unwrapping a `Result<ScoreParams, ProfileError>` at
/// compile-time.
///
/// ## Panics
///
/// Panics if `res` is an error.
#[cfg(test)]
#[allow(deprecated)]
pub const fn unwrap_params<T, const S: usize>(res: Result<ScoreParams<T, S>, ProfileError>) -> ScoreParams<T, S>
where
    T: AnyInt, {
    match res {
        Ok(params) => params,
        Err(ProfileError::EmptySequence) => panic!("an empty sequence was provided"),
        Err(ProfileError::GapOpenOutOfRange { .. }) => panic!("gap_open must be in -127..=0"),
        Err(ProfileError::GapExtendOutOfRange { .. }) => panic!("gap_extend must be in -127..=0"),
        Err(ProfileError::BadGapWeights { .. }) => panic!("gap_extend must be greater than or equal to gap_open"),
    }
}
