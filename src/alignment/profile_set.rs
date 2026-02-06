use crate::{
    alignment::{Alignment, MaybeAligned, ProfileError, ScoreAndRanges, SeqSrc, StripedProfile, validate_profile_args},
    data::matrices::WeightMatrix,
};
use std::{cell::OnceCell, sync::OnceLock};

/// A trait supporting sets of striped alignment profiles.
///
/// [`ProfileSets`] offer an abstraction around [`StripedProfile`], providing
/// convenience methods for automatically increasing the integer width and
/// rerunning the alignment when overflow occurs.
///
/// ## Type Parameters
///
/// - `M` - The number of SIMD lanes for `i8` profiles.
/// - `N` - The number of SIMD lanes for `i16` profiles.
/// - `O` - The number of SIMD lanes for `i32` profiles.
/// - `S` - The size of the alphabet (usually 5 for DNA including `N`).
pub trait ProfileSets<'a, const M: usize, const N: usize, const O: usize, const S: usize>: Sized {
    /// Gets or initializes [`StripedProfile`] with elements of `i8` and `M`
    /// SIMD lanes and returns a reference to the field.
    fn get_i8(&self) -> &StripedProfile<'a, i8, M, S>;

    /// Gets or initializes [`StripedProfile`] with elements of `i16` and `N`
    /// SIMD lanes and returns a reference to the field.
    fn get_i16(&self) -> &StripedProfile<'a, i16, N, S>;

    /// Gets or initializes [`StripedProfile`] with elements of `i32` and `O`
    /// SIMD lanes and returns a reference to the field.
    fn get_i32(&self) -> &StripedProfile<'a, i32, O, S>;

    /// Retrieves the sequence from which the profiles are built.
    fn sequence(&self) -> &[u8];

    /// Retrieves the gap open score being used.
    fn gap_open(&self) -> i8;

    /// Retrieves the gap extend score being used.
    fn gap_extend(&self) -> i8;

    /// Retrieves the weight matrix being used.
    fn matrix(&self) -> &WeightMatrix<'a, i8, S>;

    /// Lazily execute [`StripedProfile::sw_score`] starting with
    /// the `i8` profile.
    ///
    /// Lazily initializes the profiles and works its way up to the `i32`
    /// profile. Execution stops when the score returned no longer overflows the
    /// profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #    alignment::{LocalProfiles, ProfileSets, SeqSrc, sw::sw_simd_score},
    /// #    data::matrices::WeightMatrix
    /// # };
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.sw_score_from_i8(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    fn sw_score_from_i8<T>(&self, seq: &T) -> MaybeAligned<u32>
    where
        T: AsRef<[u8]> + ?Sized, {
        self.get_i8()
            .sw_score(seq)
            .or_else_overflowed(|| self.get_i16().sw_score(seq))
            .or_else_overflowed(|| self.get_i32().sw_score(seq))
    }

    /// Lazily execute [`StripedProfile::sw_score`] starting with
    /// the `i16` profile, skipping the `i8` profile.
    ///
    /// Lazily initializes the profiles and works its way up to the `i32`
    /// profile. Execution stops when the score returned no longer overflows the
    /// profile's integer range.
    ///
    /// See [`LocalProfiles::sw_score_from_i8`] for an example.
    #[inline]
    #[must_use]
    fn sw_score_from_i16<T>(&self, seq: &T) -> MaybeAligned<u32>
    where
        T: AsRef<[u8]> + ?Sized, {
        self.get_i16()
            .sw_score(seq)
            .or_else_overflowed(|| self.get_i32().sw_score(seq))
    }

    /// Execute [`StripedProfile::sw_score`] with the `i32` profile,
    /// skipping the `i8` and `i16` profiles.
    ///
    /// See [`LocalProfiles::sw_score_from_i8`] for an example.
    #[inline]
    #[must_use]
    fn sw_score_from_i32<T>(&self, seq: &T) -> MaybeAligned<u32>
    where
        T: AsRef<[u8]> + ?Sized, {
        self.get_i32().sw_score(seq)
    }

    /// Lazily execute [`StripedProfile::sw_align`] starting
    /// with the `i8` profile.
    ///
    /// Lazily initializes the profiles and works its way up to the `i32`
    /// profile. Execution stops when the alignment returned no longer overflows
    /// the profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     alignment::{LocalProfiles, ProfileSets, SeqSrc, sw::sw_simd_score},
    /// #     data::matrices::WeightMatrix
    /// # };
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let alignment = profile.sw_align_from_i8(SeqSrc::Reference(reference)).unwrap();
    /// ```
    #[inline]
    #[must_use]
    fn sw_align_from_i8<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<Alignment<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i8()
            .sw_align(seq)
            .or_else_overflowed(|| self.get_i16().sw_align(seq))
            .or_else_overflowed(|| self.get_i32().sw_align(seq))
    }

    /// Lazily execute [`StripedProfile::sw_align`] starting
    /// with the `i16` profile, skipping the `i8` profile.
    ///
    /// Lazily initializes the profiles and works its way up to the `i32`
    /// profile. Execution stops when the alignment returned no longer overflows
    /// the profile's integer range.
    ///
    /// See [`LocalProfiles::sw_align_from_i8`] for an example.
    #[inline]
    #[must_use]
    fn sw_align_from_i16<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<Alignment<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i16()
            .sw_align(seq)
            .or_else_overflowed(|| self.get_i32().sw_align(seq))
    }

    /// Execute [`StripedProfile::sw_align`] with the `i32`
    /// profile, skipping the `i8` and `i16` profiles.
    ///
    /// See [`LocalProfiles::sw_align_from_i8`] for an example.
    #[inline]
    #[must_use]
    fn sw_align_from_i32<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<Alignment<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i32().sw_align(seq)
    }

    // TODO: we will add dispatching instead if the method needs to be hybrid based on size considerations
    /// Lazily executes a 3-pass version of Smith-Waterman local alignment,
    /// starting with an `i8` profile.
    ///
    /// For more details on the algorithm, see [`sw_align_3pass`].
    ///
    /// The method automatically increases the integer width until `i32` upon
    /// each overflow. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #    alignment::{LocalProfiles, ProfileSets, SeqSrc, sw::sw_simd_score},
    /// #    data::matrices::WeightMatrix
    /// # };
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.sw_align_from_i8_3pass(SeqSrc::Reference(reference)).unwrap().score;
    /// assert_eq!(score, 26);
    /// ```
    ///
    /// [`sw_align_3pass`]: crate::alignment::sw::sw_align_3pass
    #[inline]
    #[must_use]
    fn sw_align_from_i8_3pass<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<Alignment<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i8()
            .sw_align_3pass(seq, self.sequence(), self.matrix(), self.gap_open(), self.gap_extend())
            .or_else_overflowed(|| {
                self.get_i16()
                    .sw_align_3pass(seq, self.sequence(), self.matrix(), self.gap_open(), self.gap_extend())
            })
            .or_else_overflowed(|| {
                self.get_i32()
                    .sw_align_3pass(seq, self.sequence(), self.matrix(), self.gap_open(), self.gap_extend())
            })
    }

    // TODO: we will add dispatching instead if the method needs to be hybrid based on size considerations
    /// Lazily executes a 3-pass version of Smith-Waterman local alignment,
    /// starting with an `i16` profile.
    ///
    /// For more details on the algorithm, see [`sw_align_3pass`]. See
    /// [`sw_align_from_i8_3pass`] for an example.
    ///
    /// The method automatically increases the integer width until `i32` upon
    /// each overflow. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    ///
    /// [`sw_align_from_i8_3pass`]:
    ///     ProfileSets::sw_align_from_i8_3pass
    /// [`sw_align_3pass`]: crate::alignment::sw::sw_align_3pass
    #[inline]
    #[must_use]
    fn sw_align_from_i16_3pass<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<Alignment<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i16()
            .sw_align_3pass(seq, self.sequence(), self.matrix(), self.gap_open(), self.gap_extend())
            .or_else_overflowed(|| {
                self.get_i32()
                    .sw_align_3pass(seq, self.sequence(), self.matrix(), self.gap_open(), self.gap_extend())
            })
    }

    // TODO: we will add dispatching instead if the method needs to be hybrid based on size considerations
    /// Lazily executes a 3-pass version of Smith-Waterman local alignment using
    /// an `i32` profile.
    ///
    /// For more details on the algorithm, see [`sw_align_3pass`]. See
    /// [`sw_align_from_i8_3pass`] for an example.
    ///
    /// If the score overflows the range allowed by an `i32`, then
    /// [`MaybeAligned::Overflowed`] is returned, since this is the highest
    /// profile supported by *Zoe*'s profile sets.
    ///
    /// [`sw_align_from_i8_3pass`]:
    ///     ProfileSets::sw_align_from_i8_3pass
    /// [`sw_align_3pass`]: crate::alignment::sw::sw_align_3pass
    #[inline]
    #[must_use]
    fn sw_align_from_i32_3pass<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<Alignment<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i32()
            .sw_align_3pass(seq, self.sequence(), self.matrix(), self.gap_open(), self.gap_extend())
    }

    /// Lazily execute [`StripedProfile::sw_score_ranges`] starting
    /// with the `i8` profile.
    ///
    /// Lazily initializes the profiles and works its way up to the `i32`
    /// profile. Execution stops when the alignment returned no longer overflows
    /// the profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     alignment::{LocalProfiles, ProfileSets, SeqSrc, sw::sw_simd_score},
    /// #     data::matrices::WeightMatrix
    /// # };
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let scores_and_ranges = profile.sw_score_ranges_from_i8(SeqSrc::Reference(reference)).unwrap();
    /// assert_eq!(scores_and_ranges.score, 26);
    /// assert_eq!(scores_and_ranges.query_range, 0..15);
    /// assert_eq!(scores_and_ranges.ref_range, 14..31);
    /// ```
    #[inline]
    #[must_use]
    fn sw_score_ranges_from_i8<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<ScoreAndRanges<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i8()
            .sw_score_ranges(seq)
            .or_else_overflowed(|| self.get_i16().sw_score_ranges(seq))
            .or_else_overflowed(|| self.get_i32().sw_score_ranges(seq))
    }

    /// Lazily execute [`StripedProfile::sw_score_ranges`] starting
    /// with the `i16` profile, skipping the `i8` profile.
    ///
    /// Lazily initializes the profiles and works its way up to the `i32`
    /// profile. Execution stops when the alignment returned no longer overflows
    /// the profile's integer range.
    ///
    /// See [`LocalProfiles::sw_score_ranges_from_i8`] for an
    /// example.
    #[inline]
    #[must_use]
    fn sw_score_ranges_from_i16<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<ScoreAndRanges<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i16()
            .sw_score_ranges(seq)
            .or_else_overflowed(|| self.get_i32().sw_score_ranges(seq))
    }

    /// Execute [`StripedProfile::sw_score_ranges`] with the `i32`
    /// profile, skipping the `i8` and `i16` profiles.
    ///
    /// See [`LocalProfiles::sw_score_ranges_from_i8`] for an
    /// example.
    #[inline]
    #[must_use]
    fn sw_score_ranges_from_i32<T>(&self, seq: SeqSrc<&T>) -> MaybeAligned<ScoreAndRanges<u32>>
    where
        T: AsRef<[u8]> + ?Sized, {
        let seq = seq.map(AsRef::as_ref);

        self.get_i32().sw_score_ranges(seq)
    }
}

/// A lazily-evaluated set of striped alignment profiles for local
/// (thread-specific) use.
///
/// Provides convenience methods around [`StripedProfile`] for automatically
/// increasing the integer width and rerunning the alignment when overflow
/// occurs. This implements signed versions of the algorithms.
///
/// If it is necessary to share between multiple threads, consider using
/// [`SharedProfiles`].
///
/// ## Type Parameters
///
/// - `M` - The number of SIMD lanes for `i8` profiles.
/// - `N` - The number of SIMD lanes for `i16` profiles.
/// - `O` - The number of SIMD lanes for `i32` profiles.
/// - `S` - The size of the alphabet (usually 5 for DNA including `N`).
#[derive(Debug, Clone)]
pub struct LocalProfiles<'a, const M: usize, const N: usize, const O: usize, const S: usize> {
    pub(crate) seq:         &'a [u8],
    pub(crate) matrix:      &'a WeightMatrix<'a, i8, S>,
    pub(crate) gap_open:    i8,
    pub(crate) gap_extend:  i8,
    pub(crate) profile_i8:  OnceCell<StripedProfile<'a, i8, M, S>>,
    pub(crate) profile_i16: OnceCell<StripedProfile<'a, i16, N, S>>,
    pub(crate) profile_i32: OnceCell<StripedProfile<'a, i32, O, S>>,
}

impl<'a, const M: usize, const N: usize, const O: usize, const S: usize> LocalProfiles<'a, M, N, O, S> {
    /// Creates an empty [`LocalProfiles`].
    ///
    /// Usually you instead want to use [`new_with_w128`], [`new_with_w256`], or
    /// [`new_with_w512`], based on the SIMD register width of your target.
    /// [`new_with_w256`] is a good default.
    ///
    /// ## Errors
    ///
    /// - [`ProfileError::EmptySequence`] if `seq` is empty
    /// - [`ProfileError::GapOpenOutOfRange`] if `gap_open` is not between -127
    ///   and 0, inclusive
    /// - [`ProfileError::GapExtendOutOfRange`] if `gap_extend` is not between
    ///   -127 and 0, inclusive
    /// - [`ProfileError::BadGapWeights`] if `gap_extend` is less than
    ///   `gap_open`
    ///
    /// [`new_with_w128`]: LocalProfiles::new_with_w128
    /// [`new_with_w256`]: LocalProfiles::new_with_w256
    /// [`new_with_w512`]: LocalProfiles::new_with_w512
    #[inline]
    pub fn new<T>(
        seq: &'a T, matrix: &'a WeightMatrix<'a, i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, ProfileError>
    where
        T: AsRef<[u8]> + ?Sized, {
        validate_profile_args(seq.as_ref(), gap_open, gap_extend)?;

        Ok(LocalProfiles {
            seq: seq.as_ref(),
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceCell::new(),
            profile_i16: OnceCell::new(),
            profile_i32: OnceCell::new(),
        })
    }
}

impl<'a, const S: usize> LocalProfiles<'a, 16, 8, 4, S> {
    /// Creates an empty [`LocalProfiles`] optimized for 128-bit SIMD width.
    ///
    /// This sets `M=16`, `N=8`, and `O=4` for `i8`, `i16`, and `i32` profiles
    /// respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`LocalProfiles::new`].
    #[inline]
    pub fn new_with_w128<T>(
        seq: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, ProfileError>
    where
        T: AsRef<[u8]> + ?Sized, {
        Self::new(seq, matrix, gap_open, gap_extend)
    }
}

impl<'a, const S: usize> LocalProfiles<'a, 32, 16, 8, S> {
    /// Creates an empty [`LocalProfiles`] optimized for 256-bit SIMD width.
    ///
    /// This sets `M=32`, `N=16`, and `O=8` for `i8`, `i16`, and `i32` profiles
    /// respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`LocalProfiles::new`].
    #[inline]
    pub fn new_with_w256<T>(
        seq: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, ProfileError>
    where
        T: AsRef<[u8]> + ?Sized, {
        Self::new(seq, matrix, gap_open, gap_extend)
    }
}

impl<'a, const S: usize> LocalProfiles<'a, 64, 32, 16, S> {
    /// Creates an empty [`LocalProfiles`] optimized for 512-bit SIMD width.
    ///
    /// This sets `M=64`, `N=32`, and `O=16` for `i8`, `i16`, and `i32` profiles
    /// respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`LocalProfiles::new`].
    #[inline]
    pub fn new_with_w512<T>(
        seq: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, ProfileError>
    where
        T: AsRef<[u8]> + ?Sized, {
        Self::new(seq, matrix, gap_open, gap_extend)
    }
}

impl<'a, const M: usize, const N: usize, const O: usize, const S: usize> ProfileSets<'a, M, N, O, S>
    for LocalProfiles<'a, M, N, O, S>
{
    #[inline]
    fn get_i8(&self) -> &StripedProfile<'a, i8, M, S> {
        // Validity: We already validated profile
        self.profile_i8
            .get_or_init(|| StripedProfile::new_unchecked(self.seq, self.matrix, self.gap_open, self.gap_extend))
    }

    #[inline]
    fn get_i16(&self) -> &StripedProfile<'a, i16, N, S> {
        // Validity: We already validated profile
        self.profile_i16
            .get_or_init(|| StripedProfile::new_unchecked(self.seq, self.matrix, self.gap_open, self.gap_extend))
    }

    #[inline]
    fn get_i32(&self) -> &StripedProfile<'a, i32, O, S> {
        // Validity: We already validated profile
        self.profile_i32
            .get_or_init(|| StripedProfile::new_unchecked(self.seq, self.matrix, self.gap_open, self.gap_extend))
    }

    #[inline]
    fn sequence(&self) -> &[u8] {
        self.seq
    }

    #[inline]
    fn gap_open(&self) -> i8 {
        self.gap_open
    }

    #[inline]
    fn gap_extend(&self) -> i8 {
        self.gap_extend
    }

    #[inline]
    fn matrix(&self) -> &WeightMatrix<'a, i8, S> {
        self.matrix
    }
}

/// A lazily-evaluated set of striped alignment profiles which can be shared
/// across threads.
///
/// This is an abstraction around [`StripedProfile`], providing convenience
/// methods for automatically increasing the integer width and rerunning the
/// alignment when overflow occurs. This only supports the unsigned version of
/// the algorithm.
///
/// When sharing between threads is not needed, consider using [`LocalProfiles`]
/// instead.
///
/// ## Type Parameters
///
/// - `M` - The number of SIMD lanes for `i8` profiles
/// - `N` - The number of SIMD lanes for `i16` profiles
/// - `O` - The number of SIMD lanes for `i32` profiles
/// - `S` - The size of the alphabet (usually 5 for DNA including *N*)
#[derive(Debug, Clone)]
pub struct SharedProfiles<'a, const M: usize, const N: usize, const O: usize, const S: usize> {
    pub(crate) seq:         &'a [u8],
    pub(crate) matrix:      &'a WeightMatrix<'a, i8, S>,
    pub(crate) gap_open:    i8,
    pub(crate) gap_extend:  i8,
    pub(crate) profile_i8:  OnceLock<StripedProfile<'a, i8, M, S>>,
    pub(crate) profile_i16: OnceLock<StripedProfile<'a, i16, N, S>>,
    pub(crate) profile_i32: OnceLock<StripedProfile<'a, i32, O, S>>,
}

impl<'a, const M: usize, const N: usize, const O: usize, const S: usize> SharedProfiles<'a, M, N, O, S> {
    /// Creates an empty [`SharedProfiles`].
    ///
    /// Usually you instead want to use [`new_with_w128`], [`new_with_w256`], or
    /// [`new_with_w512`], based on the SIMD register width of your target.
    /// [`new_with_w256`] is a good default.
    ///
    /// ## Errors
    ///
    /// - [`ProfileError::EmptySequence`] if `seq` is empty
    /// - [`ProfileError::GapOpenOutOfRange`] if `gap_open` is not between -127
    ///   and 0, inclusive
    /// - [`ProfileError::GapExtendOutOfRange`] if `gap_extend` is not between
    ///   -127 and 0, inclusive
    /// - [`ProfileError::BadGapWeights`] if `gap_extend` is less than
    ///   `gap_open`
    ///
    /// [`new_with_w128`]: SharedProfiles::new_with_w128
    /// [`new_with_w256`]: SharedProfiles::new_with_w256
    /// [`new_with_w512`]: SharedProfiles::new_with_w512
    #[inline]
    pub fn new<T>(
        seq: &'a T, matrix: &'a WeightMatrix<'a, i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, ProfileError>
    where
        T: AsRef<[u8]> + ?Sized, {
        validate_profile_args(seq.as_ref(), gap_open, gap_extend)?;

        Ok(SharedProfiles {
            seq: seq.as_ref(),
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceLock::new(),
            profile_i16: OnceLock::new(),
            profile_i32: OnceLock::new(),
        })
    }
}

impl<'a, const S: usize> SharedProfiles<'a, 16, 8, 4, S> {
    /// Creates an empty [`SharedProfiles`] optimized for 128-bit SIMD width.
    ///
    /// This sets `M=16`, `N=8`, and `O=4` for `i8`, `i16`, and `i32` profiles
    /// respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`SharedProfiles::new`].
    #[inline]
    pub fn new_with_w128<T>(
        seq: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<SharedProfiles<'a, 16, 8, 4, S>, ProfileError>
    where
        T: AsRef<[u8]> + ?Sized, {
        Self::new(seq, matrix, gap_open, gap_extend)
    }
}

impl<'a, const S: usize> SharedProfiles<'a, 32, 16, 8, S> {
    /// Creates an empty [`SharedProfiles`] optimized for 256-bit SIMD width.
    ///
    /// This sets `M=32`, `N=16`, and `O=8` for `i8`, `i16`, and `i32` profiles
    /// respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`SharedProfiles::new`].
    #[inline]
    pub fn new_with_w256<T>(
        seq: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<SharedProfiles<'a, 32, 16, 8, S>, ProfileError>
    where
        T: AsRef<[u8]> + ?Sized, {
        Self::new(seq, matrix, gap_open, gap_extend)
    }
}

impl<'a, const S: usize> SharedProfiles<'a, 64, 32, 16, S> {
    /// Creates an empty [`SharedProfiles`] optimized for 512-bit SIMD width.
    ///
    /// This sets `M=64`, `N=32`, and `O=16` for `i8`, `i16`, and `i32` profiles
    /// respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`SharedProfiles::new`].
    #[inline]
    pub fn new_with_w512<T>(
        seq: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<SharedProfiles<'a, 64, 32, 16, S>, ProfileError>
    where
        T: AsRef<[u8]> + ?Sized, {
        Self::new(seq, matrix, gap_open, gap_extend)
    }
}

impl<'a, const M: usize, const N: usize, const O: usize, const S: usize> ProfileSets<'a, M, N, O, S>
    for SharedProfiles<'a, M, N, O, S>
{
    #[inline]
    fn get_i8(&self) -> &StripedProfile<'a, i8, M, S> {
        // Validity: We already validated profile
        self.profile_i8
            .get_or_init(|| StripedProfile::new_unchecked(self.seq, self.matrix, self.gap_open, self.gap_extend))
    }

    #[inline]
    fn get_i16(&self) -> &StripedProfile<'a, i16, N, S> {
        // Validity: We already validated profile
        self.profile_i16
            .get_or_init(|| StripedProfile::new_unchecked(self.seq, self.matrix, self.gap_open, self.gap_extend))
    }

    #[inline]
    fn get_i32(&self) -> &StripedProfile<'a, i32, O, S> {
        // Validity: We already validated profile
        self.profile_i32
            .get_or_init(|| StripedProfile::new_unchecked(self.seq, self.matrix, self.gap_open, self.gap_extend))
    }

    #[inline]
    fn sequence(&self) -> &[u8] {
        self.seq
    }

    #[inline]
    fn gap_open(&self) -> i8 {
        self.gap_open
    }

    #[inline]
    fn gap_extend(&self) -> i8 {
        self.gap_extend
    }

    #[inline]
    fn matrix(&self) -> &WeightMatrix<'a, i8, S> {
        self.matrix
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::alignment::sw::test_data::{GAP_EXTEND, GAP_OPEN, WEIGHTS};

    #[allow(clippy::cast_possible_wrap)]
    #[test]
    fn sw_simd_profile_set() {
        let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
        let profiles = LocalProfiles::new_with_w256(&v, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        let profile1 = profiles.get_i8();
        let profile2 = StripedProfile::<i8, 32, 5>::new(v, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        assert_eq!(profile1, &profile2);
    }
}
