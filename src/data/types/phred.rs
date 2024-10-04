#![allow(clippy::cast_precision_loss)]

use crate::data::vec_types::CheckSequence;
use std::io::{Error as IOError, ErrorKind};

#[derive(Debug, Default, Clone, PartialEq)]
pub struct QualityScores(pub(crate) Vec<EncodedQS>);

/// Alias for encoded Quality Score data
type EncodedQS = u8;

impl QualityScores {
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        QualityScores(Vec::new())
    }

    #[inline]
    #[must_use]
    /// Creates a [`QualityScores`] struct from a `Vec<u8>` without checking for
    /// values being in the proper range (see [`is_graphic_simd`]).
    ///
    /// # Safety
    /// The [`Vec<u8>`] must be in range `!`..=`~`.
    ///
    /// [`is_graphic_simd`]: crate::data::vec_types::CheckSequence::is_graphic_simd
    pub unsafe fn from_vec_unchecked(v: Vec<EncodedQS>) -> Self {
        QualityScores(v)
    }

    /// Truncates the length of the sequence to the specified `new_length`.
    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
    }

    /// Cuts the 5' end of the [`QualityScores`] just prior to the new starting
    /// index (0-based). Be aware that this clones the internal buffer!
    #[inline]
    pub fn cut_to_start(&mut self, new_start: usize) {
        *self = QualityScores(self.0.drain(new_start..).collect());
    }

    /// Reverse quality scores.
    #[must_use]
    #[inline]
    pub fn reverse(&self) -> Self {
        Self(self.0.iter().rev().copied().collect())
    }

    /// Minimum phred quality score
    #[inline]
    #[must_use]
    pub fn min_q(&self) -> Option<QScoreInt> {
        self.0.iter().min().map(|q| QScoreInt(q - 33))
    }

    /// Maximum phred quality score
    #[inline]
    #[must_use]
    pub fn max_q(&self) -> Option<QScoreInt> {
        self.0.iter().max().map(|q| QScoreInt(q - 33))
    }

    /// Produces a sorted copy of Quality Scores
    #[inline]
    #[must_use]
    pub fn sorted(&self) -> Self {
        let mut s = self.0.clone();
        crate::sort::count_sort_ascii_print(&mut s);
        QualityScores(s)
    }

    /// The median phred quality score.
    #[inline]
    #[must_use]
    pub fn median(&self) -> Option<QScoreFloat> {
        match self.0.len() {
            0 => None,
            1 => Some(self.0[0].into()),
            n => {
                let s = self.sorted();
                if n % 2 == 0 {
                    Some(QScoreFloat((f32::from(s.0[n / 2] + s.0[(n / 2) - 1]) / 2.0) - 33.0))
                } else {
                    Some(s.0[(n - 1) / 2].into())
                }
            }
        }
    }

    /// Efficiently obtain tuple of (min, median, max)
    #[inline]
    #[must_use]
    pub fn min_median_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        match self.0.len() {
            0 => None,
            1 => Some((self.0[0].into(), self.0[0].into(), self.0[0].into())),
            n => {
                let s = self.sorted();
                if n % 2 == 0 {
                    Some((
                        QScoreInt::from(s.0[0]),
                        QScoreFloat((f32::from(s.0[n / 2] + s.0[(n / 2) - 1]) / 2.0) - 33.0),
                        QScoreInt::from(s.0[n - 1]),
                    ))
                } else {
                    Some((
                        QScoreInt::from(s.0[0]),
                        QScoreFloat::from(s.0[(n - 1) / 2]),
                        QScoreInt::from(s.0[n - 1]),
                    ))
                }
            }
        }
    }

    /// The geometric mean error represented as a phred quality score floating
    /// point value; equivalent to the arithmetic mean quality score.
    #[inline]
    #[must_use]
    pub fn geometric_mean(&self) -> Option<QScoreFloat> {
        match self.0.len() {
            0 => None,
            1 => Some(self.0[0].into()),
            n => {
                let sum: usize = self.0.iter().fold(0usize, |sum, x| sum + *x as usize);
                Some(QScoreFloat((sum - (n * 33)) as f32 / n as f32))
            }
        }
    }

    /// Efficiently obtain tuple of min, geometric mean, and max.
    #[inline]
    #[must_use]
    pub fn min_geomean_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        match self.0.len() {
            0 => None,
            1 => Some((self.0[0].into(), self.0[0].into(), self.0[0].into())),
            n => {
                let (min, mean, max) = self.0.iter().fold((self.0[0], 0usize, self.0[0]), |(min, mean, max), x| {
                    let new_min = if *x < min { *x } else { min };
                    let new_max = if *x > min { *x } else { max };
                    (new_min, mean + *x as usize, new_max)
                });
                Some((
                    QScoreInt::from(min),
                    QScoreFloat((mean - (n * 33)) as f32 / n as f32),
                    QScoreInt::from(max),
                ))
            }
        }
    }

    /// The arithetic mean error represented as a phred quality score floating point value.
    #[inline]
    #[must_use]
    pub fn arithmetic_mean(&self) -> Option<QScoreFloat> {
        match self.0.len() {
            0 => None,
            1 => Some(self.0[0].into()),
            n => {
                let sum = self.0.iter().fold(0f32, |sum, x| sum + QScoreInt::encoded_qs_to_error(*x));

                Some(QScoreInt::error_to_q(sum / n as f32))
            }
        }
    }

    /// Efficiently get a tuple of min, arithmetic mean, and maximum.
    #[inline]
    #[must_use]
    pub fn min_average_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        match self.0.len() {
            0 => None,
            1 => Some((self.0[0].into(), self.0[0].into(), self.0[0].into())),
            n => {
                let (min, sum, max) = self.0.iter().fold((self.0[0], 0f32, self.0[0]), |(min, sum, max), x| {
                    let new_min = if *x < min { *x } else { min };
                    let new_max = if *x > min { *x } else { max };
                    (new_min, sum + QScoreInt::encoded_qs_to_error(*x), new_max)
                });

                Some((
                    QScoreInt::from(min),
                    QScoreInt::error_to_q(sum / n as f32),
                    QScoreInt::from(max),
                ))
            }
        }
    }

    /// Length of the [`QualityScores`]
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Is the [`QualityScores`] vector empty?
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Provide [`QualityScores`] as a reference to a raw byte slice ([`EncodedQS`]).
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_slice()
    }

    /// Provide [`QualityScores`] as a mutable reference to a raw byte slice ([`EncodedQS`])
    #[inline]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        self.0.as_mut_slice()
    }

    /// Consumes the [`EncodedQS`] and returns a [`Vec<u8>`].
    #[inline]
    #[must_use]
    pub fn into_vec(self) -> Vec<u8> {
        self.0
    }

    /// Get quality scores for some index or range, returning an optional slice.
    #[inline]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }
}

impl std::fmt::Display for QualityScores {
    // Safety: `QualityScores` is guaranteed to be utf8 because we checked that
    // it is "graphic ascii", which will be valid utf8.
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", unsafe { std::str::from_utf8_unchecked(&self.0) })
    }
}

impl TryFrom<Vec<EncodedQS>> for QualityScores {
    type Error = IOError;
    fn try_from(encoded_scores: Vec<EncodedQS>) -> Result<Self, Self::Error> {
        if encoded_scores.is_graphic_simd::<16>() {
            Ok(QualityScores(encoded_scores))
        } else {
            Err(IOError::new(ErrorKind::InvalidData, "Quality scores contain invalid state!"))
        }
    }
}

impl TryFrom<&[EncodedQS]> for QualityScores {
    type Error = IOError;
    fn try_from(encoded_scores: &[EncodedQS]) -> Result<Self, Self::Error> {
        if encoded_scores.is_graphic_simd::<16>() {
            Ok(QualityScores(encoded_scores.to_vec()))
        } else {
            Err(IOError::new(ErrorKind::InvalidData, "Quality scores contain invalid state!"))
        }
    }
}

impl std::ops::Index<usize> for QualityScores {
    type Output = EncodedQS;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl std::ops::IndexMut<usize> for QualityScores {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl std::ops::Index<std::ops::Range<usize>> for QualityScores {
    type Output = [u8];

    #[inline]
    fn index(&self, index: std::ops::Range<usize>) -> &[u8] {
        &self.0[index]
    }
}

/// Type alias for the quality score when returned as `f32`
///
///  Significant figures for Phred error probabilities is 6 (0.999999).
/// The quality score, being a log of the Phred error times a constant, must retain 6
/// significant figures after the decimal. The maximum value is 60, so 8 digits
/// must be retained at maximum. When converting `usize`/`u64` to `f32`, 8 digits are
/// retained. Therefore, for type conversion and averages, which propogates
/// significant figures, `f32` should be a sufficient choice.
#[repr(transparent)]
pub struct QScoreFloat(pub(crate) f32);

/// Type alias for probabilities as `f32`
type Probability = f32;

impl QScoreFloat {
    /// Obtain the f32 value for the [`QScoreFloat`]
    #[inline]
    #[must_use]
    pub fn as_f32(&self) -> f32 {
        self.0
    }
}

impl From<EncodedQS> for QScoreFloat {
    fn from(encoded: EncodedQS) -> Self {
        QScoreFloat(f32::from(encoded - 33))
    }
}

/// A single phred quality score in numeric (not ASCII encoded) range.
/// The type is functionally a `u8`.
#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct QScoreInt(pub(crate) u8);

impl QScoreInt {
    /// Calculated the expected rate of machine error from an encoded quality
    /// score value (`u8`).
    #[inline]
    #[must_use]
    pub fn encoded_qs_to_error(q: EncodedQS) -> Probability {
        const BASE: f32 = 10.0;
        BASE.powf(-f32::from(q - 33) / BASE)
    }

    /// Calculate the [`QScoreFloat`] from the expected rate of machine error.
    #[inline]
    #[must_use]
    pub fn error_to_q(e: Probability) -> QScoreFloat {
        const BASE: f32 = 10.0;
        QScoreFloat(-BASE * e.log10())
    }

    /// Calculate the expected rate of machine error from the [`QScoreInt`]
    #[inline]
    #[must_use]
    pub fn to_error(&self) -> Probability {
        const BASE: f32 = 10.0;
        BASE.powf(-f32::from(self.0) / BASE)
    }

    /// Obtain the u8 value for the [`QScoreInt`]
    #[inline]
    #[must_use]
    pub fn as_u8(&self) -> u8 {
        self.0
    }
}

impl From<EncodedQS> for QScoreInt {
    fn from(encoded: EncodedQS) -> Self {
        QScoreInt(encoded - 33)
    }
}

impl std::fmt::Display for QScoreInt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
