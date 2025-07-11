#![allow(clippy::cast_precision_loss)]

use crate::prelude::*;

mod stats;
mod std_traits;
mod view_traits;

pub use stats::*;

/// [`QualityScores`] is a transparent, new-type wrapper around [`Vec<u8>`] that
/// representing Phred quality scores. It is guaranteed to contain graphic
/// ASCII in range `!`..=`~`.
#[derive(Clone, Eq, PartialEq, Hash, Default)]
#[repr(transparent)]
pub struct QualityScores(pub(crate) Vec<EncodedQS>);

/// The corresponding immutable view type for [`QualityScores`]. See
/// [Views](crate::data#views) for more details. It is guaranteed to contain
/// graphic ASCII in range `!`..=`~`.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[repr(transparent)]
pub struct QualityScoresView<'a>(pub(crate) &'a [u8]);

/// The corresponding mutable view type for [`QualityScores`]. See
/// [Views](crate::data#views) for more details. It is guaranteed to contain
/// graphic ASCII in range `!`..=`~`.
#[derive(Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[repr(transparent)]
pub struct QualityScoresViewMut<'a>(pub(crate) &'a mut [u8]);

/// Alias for encoded quality score data.
type EncodedQS = u8;

impl QualityScores {
    // Conversions and indexing

    /// Creates a new [`QualityScores`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        QualityScores(Vec::new())
    }

    /// Creates a [`QualityScores`] struct from a `Vec<u8>` without checking for
    /// values being in the proper range (see [`is_graphic_simd`]).
    ///
    /// ## Safety
    ///
    /// The [`Vec<u8>`] must be in range `!`..=`~`.
    ///
    /// [`is_graphic_simd`]:
    ///     crate::data::validation::CheckSequence::is_graphic_simd
    #[inline]
    #[must_use]
    pub unsafe fn from_vec_unchecked(v: Vec<EncodedQS>) -> Self {
        QualityScores(v)
    }

    /// Consumes [`QualityScores`] and returns a [`Vec<u8>`].
    #[inline]
    #[must_use]
    pub fn into_vec(self) -> Vec<u8> {
        self.0
    }

    /// Gets the ASCII-encoded quality scores as a byte slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    /// Gets the ASCII encoded quality scores as a mutable byte slice.
    ///
    /// ## Safety
    ///
    /// Any modifications to the bytes must ensure that they remain graphic
    /// ASCII in range `!`..=`~`.
    #[inline]
    #[must_use]
    pub unsafe fn as_mut_bytes(&mut self) -> &mut [u8] {
        &mut self.0
    }

    /// Gets the ASCII-encoded quality score or byte slice at the zero-based
    /// index, returning an [`Option`].
    #[inline]
    #[must_use]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// Creates an iterator over the ASCII encoded quality scores as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }

    // No iter_mut method currently. Note that clippy will complain if we add
    // it.

    // Manipulation

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

    /// Clear the quality scores sequence so that it is empty.
    #[inline]
    pub fn clear(&mut self) {
        self.0.clear();
    }

    // Quality-scores-specific methods

    /// Returns the sequence in reverse.
    #[inline]
    #[must_use]
    pub fn to_reverse(&self) -> QualityScores {
        QualityScores(self.iter().rev().copied().collect())
    }

    /// Reverses the sequence in-place.
    #[inline]
    pub fn make_reverse(&mut self) {
        self.0.reverse();
    }

    // Associated functions

    /// Calculated the expected rate of machine error from an encoded quality
    /// score value (`u8`).
    #[inline]
    #[must_use]
    pub fn encoded_qs_to_error(q: EncodedQS) -> Probability {
        const BASE: f32 = 10.0;
        BASE.powf(-f32::from(q - 33) / BASE)
    }
}

impl<'a> QualityScoresView<'a> {
    // Conversions and indexing

    /// Creates a new [`QualityScoresView`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        QualityScoresView(&[])
    }

    /// Creates a [`QualityScoresView`] struct from a byte slice without
    /// checking for values being in the proper range (see [`is_graphic_simd`]).
    ///
    /// ## Safety
    ///
    /// Bytes in `v` must be in range `!`..=`~`.
    ///
    /// [`is_graphic_simd`]:
    ///     crate::data::validation::CheckSequence::is_graphic_simd
    #[inline]
    #[must_use]
    pub unsafe fn from_bytes_unchecked(v: &'a [u8]) -> Self {
        QualityScoresView(v)
    }

    /// Gets the ASCII encoded quality scores as a byte slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &'a [u8] {
        self.0
    }

    /// Gets the ASCII encoded quality score or byte slice at the zero-based
    /// index, returning an [`Option`].
    #[inline]
    #[must_use]
    pub fn get<I>(&self, index: I) -> Option<&'a I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// Creates an iterator over the ASCII encoded quality scores as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'a, u8> {
        self.0.iter()
    }

    // Quality-scores-specific methods

    /// Returns the sequence in reverse.
    #[inline]
    #[must_use]
    pub fn to_reverse(&self) -> QualityScores {
        QualityScores(self.iter().rev().copied().collect())
    }
}

impl<'a> QualityScoresViewMut<'a> {
    // Conversions and indexing

    /// Creates a new [`QualityScoresViewMut`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        QualityScoresViewMut(&mut [])
    }

    /// Creates a [`QualityScoresViewMut`] struct from a byte slice without
    /// checking for values being in the proper range (see [`is_graphic_simd`]).
    ///
    /// ## Safety
    ///
    /// Bytes in `v` must be in range `!`..=`~`.
    ///
    /// [`is_graphic_simd`]:
    ///     crate::data::validation::CheckSequence::is_graphic_simd
    #[inline]
    #[must_use]
    pub unsafe fn from_bytes_unchecked(v: &'a mut [u8]) -> Self {
        QualityScoresViewMut(v)
    }

    /// Gets the ASCII encoded quality scores as a byte slice.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0
    }

    /// Gets a mutable reference to the underlying raw byte slice (ASCII encoded
    /// quality scores).
    ///
    /// ## Safety
    ///
    /// Any modifications to the bytes must ensure that they remain graphic
    /// ASCII in range `!`..=`~`.
    #[inline]
    #[must_use]
    pub unsafe fn as_mut_bytes(&mut self) -> &mut [u8] {
        self.0
    }

    /// Gets the ASCII encoded quality score or byte slice at the zero-based
    /// index, returning an [`Option`].
    #[inline]
    #[must_use]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// Creates an iterator over the ASCII encoded quality scores as `&u8`.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }

    // Quality-scores-specific methods

    /// Returns the sequence in reverse.
    #[inline]
    #[must_use]
    pub fn to_reverse(&self) -> QualityScores {
        QualityScores(self.iter().rev().copied().collect())
    }

    /// Reverses the sequence in-place.
    #[inline]
    pub fn make_reverse(&mut self) {
        self.0.reverse();
    }
}

/// Wrapper type for the quality score when returned as `f32`.
///
/// ## Significant Figures
///
/// Significant figures for Phred error probabilities is 6 (0.999999). The use
/// of [`powf`] has an unspecified precision in Rust, but empirical tests reveal
/// that an f32 achieve the 6 significant figures for quality scores in the
/// range `~`..=`J`, which covers the vast majority of use-cases. For higher
/// qualities, sometimes only 5 significant figures are preserved. For type
/// conversion and averages, which propogate significant figures, `f32` should
/// therefore be a sufficient choice.
///
/// [`powf`]: f32::powf
#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct QScoreFloat(pub(crate) f32);

/// Type alias for probabilities as `f32`.
type Probability = f32;

impl QScoreFloat {
    /// Obtains the `f32` value for the [`QScoreFloat`].
    #[inline]
    #[must_use]
    pub fn as_f32(&self) -> f32 {
        self.0
    }

    /// Calculates the [`QScoreFloat`] from the expected rate of machine error.
    #[inline]
    #[must_use]
    pub fn error_to_q(e: Probability) -> QScoreFloat {
        const BASE: f32 = 10.0;
        QScoreFloat(-BASE * e.log10())
    }
}

impl From<EncodedQS> for QScoreFloat {
    /// Converts an ASCII-encoded quality score to a [`QScoreFloat`]. This
    /// assumes that a 33 ASCII offset is used.
    #[inline]
    fn from(encoded: EncodedQS) -> Self {
        QScoreFloat(f32::from(encoded - 33))
    }
}

impl std::convert::From<f32> for QScoreFloat {
    /// Converts a `f32` representing an average/mean ASCII-encoded quality
    /// score to a [`QScoreFloat`]. This assumes that a 33 ASCII offset is used.
    #[inline]
    fn from(encoded: f32) -> Self {
        QScoreFloat(encoded - 33.0)
    }
}

/// A single phred quality score in numeric (not ASCII encoded) range. The type
/// is functionally a `u8`.
#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct QScoreInt(pub(crate) u8);

impl QScoreInt {
    /// Calculates the expected rate of machine error from the [`QScoreInt`].
    #[inline]
    #[must_use]
    pub fn to_error(&self) -> Probability {
        const BASE: f32 = 10.0;
        BASE.powf(-f32::from(self.0) / BASE)
    }

    /// Obtains the `u8` value for the [`QScoreInt`].
    #[inline]
    #[must_use]
    pub fn as_u8(&self) -> u8 {
        self.0
    }
}

impl From<EncodedQS> for QScoreInt {
    /// Converts an ASCII-encoded quality score to a [`QScoreInt`]. This assumes
    /// that a 33 ASCII offset is used.
    #[inline]
    fn from(encoded: EncodedQS) -> Self {
        QScoreInt(encoded - 33)
    }
}

impl std::fmt::Display for QScoreInt {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut buff = itoa::Buffer::new();
        f.write_str(buff.format(self.0))
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
