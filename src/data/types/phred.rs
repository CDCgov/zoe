#![allow(clippy::cast_precision_loss)]

use crate::composition::frequency_table::FrequencyTable;
use crate::data::vec_types::CheckSequence;
use std::io::{Error as IOError, ErrorKind};

#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
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

    /// Minimum phred quality score using a full counting sort
    #[inline]
    #[must_use]
    pub fn min_q(&self) -> Option<QScoreInt> {
        self.0.iter().min().map(|q| QScoreInt(q - 33))
    }

    /// Maximum phred quality score using a full counting sort
    #[inline]
    #[must_use]
    pub fn max_q(&self) -> Option<QScoreInt> {
        self.0.iter().max().map(|q| QScoreInt(q - 33))
    }

    #[inline]
    #[must_use]
    /// Calculates the median phred quality score (not ASCII encoded).
    pub fn median(&self) -> Option<QScoreFloat> {
        // Safety: we have validated that QualityScores are valid
        self.get_median_unchecked().map(|q| QScoreFloat(q - 33.0))
    }

    /// Efficiently obtain tuple of (min, median, max)
    #[inline]
    #[must_use]
    pub fn min_median_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        // Safety: we have validated that QualityScores are valid
        self.get_min_med_max_unchecked()
            .map(|(min, median, max)| (min.into(), median.into(), max.into()))
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
                let sum = self.0.iter().fold(0f32, |sum, x| sum + Self::encoded_qs_to_error(*x));

                Some(QScoreFloat::error_to_q(sum / n as f32))
            }
        }
    }

    /// Efficiently get a tuple of the minimum, arithmetic mean, and maximum.
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
                    (new_min, sum + Self::encoded_qs_to_error(*x), new_max)
                });

                Some((
                    QScoreInt::from(min),
                    QScoreFloat::error_to_q(sum / n as f32),
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

    /// Provide [`QualityScores`] as a reference to a raw byte slice (ASCII encoded quality scores).
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_slice()
    }

    /// Provide [`QualityScores`] as a mutable reference to a raw byte slice (ASCII encoded quality scores).
    #[inline]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        self.0.as_mut_slice()
    }

    /// Consumes [`QualityScores`] and returns a [`Vec<u8>`].
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

    /// Calculated the expected rate of machine error from an encoded quality
    /// score value (`u8`).
    #[inline]
    #[must_use]
    pub fn encoded_qs_to_error(q: EncodedQS) -> Probability {
        const BASE: f32 = 10.0;
        BASE.powf(-f32::from(q - 33) / BASE)
    }
}

impl FrequencyTable for QualityScores {
    const MIN: u8 = b'!';
    const MAX: u8 = b'~';
}

impl AsRef<[u8]> for QualityScores {
    fn as_ref(&self) -> &[u8] {
        &self.0
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
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
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

    /// Calculate the [`QScoreFloat`] from the expected rate of machine error.
    #[inline]
    #[must_use]
    pub fn error_to_q(e: Probability) -> QScoreFloat {
        const BASE: f32 = 10.0;
        QScoreFloat(-BASE * e.log10())
    }
}

impl From<EncodedQS> for QScoreFloat {
    fn from(encoded: EncodedQS) -> Self {
        QScoreFloat(f32::from(encoded - 33))
    }
}

impl std::convert::From<f32> for QScoreFloat {
    fn from(encoded: f32) -> Self {
        QScoreFloat(encoded - 33.0)
    }
}

/// A single phred quality score in numeric (not ASCII encoded) range.
/// The type is functionally a `u8`.
#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct QScoreInt(pub(crate) u8);

impl QScoreInt {
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

#[cfg(test)]
mod bench {
    use std::sync::LazyLock;
    extern crate test;
    use super::*;
    use crate::generate::rand_sequence;
    use test::Bencher;

    const N_A: usize = 1000;
    const N_B: usize = 1059;
    const ALPHA: &[u8] = b"aaBHHFFF@123CA0/f!";
    const BETA: &[u8] = b"aaBHHFFF@123A0/f!";
    const SEED: u64 = 11_072_024;

    static SMALL_EVEN: LazyLock<QualityScores> = LazyLock::new(|| QualityScores(b"aaBHHFFF@123CA0/f!".to_vec()));
    static SMALL_ODD: LazyLock<QualityScores> = LazyLock::new(|| QualityScores(b"aaBHHFFF@123A0/f!".to_vec()));
    static LARGE_EVEN: LazyLock<QualityScores> = LazyLock::new(|| QualityScores(rand_sequence(ALPHA, N_A, SEED)));
    static LARGE_ODD: LazyLock<QualityScores> = LazyLock::new(|| QualityScores(rand_sequence(BETA, N_B, SEED)));

    #[bench]
    fn med_small_even(b: &mut Bencher) {
        b.iter(|| SMALL_EVEN.median());
    }

    #[bench]
    fn med_small_odd(b: &mut Bencher) {
        b.iter(|| SMALL_ODD.median());
    }

    #[bench]
    fn med_large_even(b: &mut Bencher) {
        b.iter(|| LARGE_EVEN.median());
    }

    #[bench]
    fn med_large_odd(b: &mut Bencher) {
        b.iter(|| LARGE_ODD.median());
    }
}

#[cfg(test)]
mod test {
    #![allow(clippy::type_complexity)]
    use super::*;
    use crate::assert_fp_eq;

    #[test]
    fn try_into_quality_scores() {
        assert!(TryInto::<QualityScores>::try_into(b"!!!ABC~~~".as_slice()).is_ok());
        assert!(TryInto::<QualityScores>::try_into(b"\0".as_slice()).is_err());
        assert!(TryInto::<QualityScores>::try_into(b"".as_slice()).is_ok());
    }

    #[test]
    fn median() {
        let median_data = vec![
            (b"BBBBBFFFFBFFGGGGGGGBGGFHGHHGHHHFHHHGHHHFHGGHHHHHHHHHEHHHHHHGHGGCGGHHHHHHGHHFHGHHGGGHGHHHHFHGECGHGEHHHHHFHHHHHHHHHHGHGEHHHHHHGHHGHHHHHHFFFHGHHHHFHGEHHH".to_vec(), Some(39.0)),
            (b"111>1FF1FB1A11BBGFFFGFHFGGAGCGC1EA0AEHFDFGDGAEAFCEEFEEFHHFFFFGFHAFGHDFHH2@1FGF0G1BEGHHCHDGEEF0E//E//<EGEA/FCFGC<<GH2F11F@?G</?EFG<1<?FHBCGGCHHHF<CE/D0".to_vec(),Some(36.0)),
            (b"CCCCCFFFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHHHHHHHHHHHHHGGEGGGHGGGGGGGGGHGHHHHHBHHHGGHHHHGGGHGHHHGHHH0FGHHCGCC############".to_vec(), Some(39.0)),
            (b"&RV".to_vec(),Some(49.0)),
            (b"".to_vec(),None),
            (b"!!!".to_vec(),Some(0.0))
        ];

        for (v, expected) in median_data {
            let qs: QualityScores = v.try_into().unwrap();
            assert_fp_eq!(qs.median().map(|wrapper| wrapper.0), expected);
        }
    }

    #[test]
    fn min_median_max() {
        let all_data = vec![
            (b"BBBBBFFFFBFFGGGGGGGBGGFHGHHGHHHFHHHGHHHFHGGHHHHHHHHHEHHHHHHGHGGCGGHHHHHHGHHFHGHHGGGHGHHHHFHGECGHGEHHHHHFHHHHHHHHHHGHGEHHHHHHGHHGHHHHHHFFFHGHHHHFHGEHHH".to_vec(), Some((QScoreInt(33), QScoreFloat(39.0), QScoreInt(39)))),
            (b"111>1FF1FB1A11BBGFFFGFHFGGAGCGC1EA0AEHFDFGDGAEAFCEEFEEFHHFFFFGFHAFGHDFHH2@1FGF0G1BEGHHCHDGEEF0E//E//<EGEA/FCFGC<<GH2F11F@?G</?EFG<1<?FHBCGGCHHHF<CE/D0".to_vec(), Some((QScoreInt(14), QScoreFloat(36.0), QScoreInt(39)))),
            (b"CCCCCFFFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHHHHHHHHHHHHHGGEGGGHGGGGGGGGGHGHHHHHBHHHGGHHHHGGGHGHHHGHHH0FGHHCGCC############".to_vec(), Some((QScoreInt(2), QScoreFloat(39.0), QScoreInt(39)))),
            (b"&RV".to_vec(),Some((QScoreInt(5), QScoreFloat(49.0), QScoreInt(53)))),
            (b"".to_vec(),None),
            (b"!!!".to_vec(),Some((QScoreInt(0), QScoreFloat(0.0), QScoreInt(0))))
        ];

        //test min, median, max combination calculation
        for (v, expected) in all_data {
            let qs: QualityScores = v.try_into().unwrap();
            assert_eq!(qs.min_median_max(), expected);
        }
    }

    #[test]
    fn median_fuzz_regression() {
        let test_inp: Vec<u8> = [251, 22].iter_mut().map(|x| (*x % 94) + 33).collect();
        let median_val = (f32::from(test_inp[0] + test_inp[1]) / 2.0) - 33.0;
        let new_med = QualityScores(test_inp).median();
        assert_eq!(new_med.map(|wrapper| wrapper.0), Some(median_val));
    }
}
