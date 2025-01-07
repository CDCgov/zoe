use crate::composition::FrequencyTable;

use super::{Len, QScoreFloat, QScoreInt, QualityScores, QualityScoresView, QualityScoresViewMut};

pub trait QualityStats: AsRef<[u8]> + FrequencyTable + Len {
    /// Minimum phred quality score using a full counting sort
    #[inline]
    #[must_use]
    fn min_q(&self) -> Option<QScoreInt> {
        self.as_ref().iter().min().map(|q| QScoreInt(q - 33))
    }

    /// Maximum phred quality score using a full counting sort
    #[inline]
    #[must_use]
    fn max_q(&self) -> Option<QScoreInt> {
        self.as_ref().iter().max().map(|q| QScoreInt(q - 33))
    }

    /// Calculates the median phred quality score (not ASCII encoded).
    #[inline]
    #[must_use]
    fn median(&self) -> Option<QScoreFloat> {
        // Safety: we have validated that QualityScores are valid
        self.get_median_unchecked().map(|q| QScoreFloat(q - 33.0))
    }

    /// Efficiently obtain tuple of (min, median, max)
    #[inline]
    #[must_use]
    fn min_median_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        // Safety: we have validated that QualityScores are valid
        self.get_min_med_max_unchecked()
            .map(|(min, median, max)| (min.into(), median.into(), max.into()))
    }

    /// The geometric mean error represented as a phred quality score floating
    /// point value; equivalent to the arithmetic mean quality score.
    #[inline]
    #[must_use]
    fn geometric_mean(&self) -> Option<QScoreFloat> {
        match self.len() {
            0 => None,
            1 => Some(self.as_ref()[0].into()),
            n => {
                let sum: usize = self.as_ref().iter().fold(0usize, |sum, x| sum + *x as usize);
                Some(QScoreFloat((sum - (n * 33)) as f32 / n as f32))
            }
        }
    }

    /// Efficiently obtain tuple of min, geometric mean, and max.
    #[inline]
    #[must_use]
    fn min_geomean_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        match self.len() {
            0 => None,
            1 => Some((self.as_ref()[0].into(), self.as_ref()[0].into(), self.as_ref()[0].into())),
            n => {
                let (min, mean, max) =
                    self.as_ref()
                        .iter()
                        .fold((self.as_ref()[0], 0usize, self.as_ref()[0]), |(min, mean, max), x| {
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

    /// The arithetic mean error represented as a phred quality score floating
    /// point value.
    #[inline]
    #[must_use]
    fn arithmetic_mean(&self) -> Option<QScoreFloat> {
        match self.len() {
            0 => None,
            1 => Some(self.as_ref()[0].into()),
            n => {
                let sum = self
                    .as_ref()
                    .iter()
                    .fold(0f32, |sum, x| sum + QualityScores::encoded_qs_to_error(*x));

                Some(QScoreFloat::error_to_q(sum / n as f32))
            }
        }
    }

    /// Efficiently get a tuple of the minimum, arithmetic mean, and maximum.
    #[inline]
    #[must_use]
    fn min_average_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        match self.len() {
            0 => None,
            1 => Some((self.as_ref()[0].into(), self.as_ref()[0].into(), self.as_ref()[0].into())),
            n => {
                let (min, sum, max) =
                    self.as_ref()
                        .iter()
                        .fold((self.as_ref()[0], 0f32, self.as_ref()[0]), |(min, sum, max), x| {
                            let new_min = if *x < min { *x } else { min };
                            let new_max = if *x > min { *x } else { max };
                            (new_min, sum + QualityScores::encoded_qs_to_error(*x), new_max)
                        });

                Some((
                    QScoreInt::from(min),
                    QScoreFloat::error_to_q(sum / n as f32),
                    QScoreInt::from(max),
                ))
            }
        }
    }
}

impl QualityStats for QualityScores {}
impl QualityStats for QualityScoresView<'_> {}
impl QualityStats for QualityScoresViewMut<'_> {}
