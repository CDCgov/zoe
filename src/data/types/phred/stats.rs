use crate::composition::{get_median, get_min_med_max, tally_from_unchecked};

use super::{Len, QScoreFloat, QScoreInt, QualityScores, QualityScoresView, QualityScoresViewMut};
use std::convert::Into;

/// Methods for calculating statistics on quality scores. All methods assume
/// that a 33 ASCII offset is used.
pub trait QualityStats: AsRef<[u8]> + Len {
    /// Calculates the minimum phred quality score using a full scan.
    #[inline]
    #[must_use]
    fn min_q(&self) -> Option<QScoreInt> {
        self.as_ref().iter().min().copied().map(QScoreInt::from)
    }

    /// Calculates the maximum phred quality score using a full scan.
    #[inline]
    #[must_use]
    fn max_q(&self) -> Option<QScoreInt> {
        self.as_ref().iter().max().copied().map(QScoreInt::from)
    }

    /// Calculates the median phred quality score (not ASCII encoded).
    #[inline]
    #[must_use]
    fn median(&self) -> Option<QScoreFloat> {
        let mut table = [0; 256];
        // Safety: QualityScores are graphic ASCII
        let num_items = tally_from_unchecked(self, &mut table);
        get_median::<b'!'>(&table, num_items).map(Into::into)
    }

    /// Efficiently calculates a tuple of minimum quality, median quality, and
    /// max quality. See [`min_q`], [`median`], and [`max_q`] for more details.
    ///
    /// [`min_q`]: QualityStats::min_q
    /// [`median`]: QualityStats::median
    /// [`max_q`]: QualityStats::max_q
    #[inline]
    #[must_use]
    fn min_median_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        let mut table = [0; 256];
        // Safety: QualityScores are graphic ASCII
        let num_items = tally_from_unchecked(self, &mut table);
        get_min_med_max::<b'!', b'~'>(&table, num_items).map(|(min, median, max)| (min.into(), median.into(), max.into()))
    }

    /// The geometric mean error represented as a phred quality score floating
    /// point value; equivalent to the arithmetic mean quality score.
    #[inline]
    #[must_use]
    fn geometric_mean(&self) -> Option<QScoreFloat> {
        if self.is_empty() {
            return None;
        }

        let n = self.len();
        let sum = self.as_ref().iter().map(|x| *x as usize).sum::<usize>();
        Some(QScoreFloat((sum - (n * 33)) as f32 / n as f32))
    }

    /// Efficiently calculates a tuple of the minimum quality, the geometric
    /// mean, and the max quality. See [`min_q`], [`geometric_mean`], and
    /// [`max_q`] for more details.
    ///
    /// [`min_q`]: QualityStats::min_q
    /// [`geometric_mean`]: QualityStats::geometric_mean
    /// [`max_q`]: QualityStats::max_q
    #[inline]
    #[must_use]
    fn min_geomean_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        if self.is_empty() {
            return None;
        }

        let n = self.len();
        let (min, sum, max) =
            self.as_ref()
                .iter()
                .copied()
                .fold((self.as_ref()[0], 0usize, self.as_ref()[0]), |(min, sum, max), x| {
                    let min = std::cmp::min(x, min);
                    let max = std::cmp::max(x, max);
                    let sum = sum + x as usize;
                    (min, sum, max)
                });

        let mean = (sum - (n * 33)) as f32 / n as f32;
        Some((QScoreInt::from(min), QScoreFloat(mean), QScoreInt::from(max)))
    }

    /// Calculates the arithetic mean error represented as a phred quality score
    /// floating point value.
    #[inline]
    #[must_use]
    fn arithmetic_mean(&self) -> Option<QScoreFloat> {
        if self.is_empty() {
            return None;
        }

        let sum = self
            .as_ref()
            .iter()
            .copied()
            .map(QualityScores::encoded_qs_to_error)
            .sum::<f32>();

        Some(QScoreFloat::error_to_q(sum / self.len() as f32))
    }

    /// Efficiently calculates a tuple of the minimum quality, arithmetic mean
    /// quality, and maximum quality. See [`min_q`], [`arithmetic_mean`], and
    /// [`max_q`] for more details.
    ///
    /// [`min_q`]: QualityStats::min_q
    /// [`arithmetic_mean`]: QualityStats::arithmetic_mean
    /// [`max_q`]: QualityStats::max_q
    #[inline]
    #[must_use]
    fn min_average_max(&self) -> Option<(QScoreInt, QScoreFloat, QScoreInt)> {
        if self.is_empty() {
            return None;
        }

        let (min, sum, max) =
            self.as_ref()
                .iter()
                .copied()
                .fold((self.as_ref()[0], 0f32, self.as_ref()[0]), |(min, sum, max), x| {
                    let min = std::cmp::min(x, min);
                    let max = std::cmp::max(x, max);
                    let sum = sum + QualityScores::encoded_qs_to_error(x);
                    (min, sum, max)
                });

        Some((
            QScoreInt::from(min),
            QScoreFloat::error_to_q(sum / self.len() as f32),
            QScoreInt::from(max),
        ))
    }
}

impl QualityStats for QualityScores {}
impl QualityStats for QualityScoresView<'_> {}
impl QualityStats for QualityScoresViewMut<'_> {}
