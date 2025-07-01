use super::{StripedProfile, validate_profile_args};
use crate::data::{err::QueryProfileError, matrices::WeightMatrix};
use std::{
    cell::OnceCell,
    simd::{LaneCount, SupportedLaneCount},
    sync::OnceLock,
};

/// A lazily-evaluated set of striped alignment profiles for local
/// (thread-specific) use. This is an abstraction around [`StripedProfile`],
/// providing convenience methods for automatically increasing the integer width
/// and rerunning the alignment when overflow occurs. This only supports the
/// unsigned version of the algorithm.
///
/// If it is necessary to share between multiple threads, consider using
/// [`SharedProfiles`].
///
/// ## Type Parameters
///
/// * `M` - The number of SIMD lanes for i8 profiles
/// * `N` - The number of SIMD lanes for i16 profiles
/// * `O` - The number of SIMD lanes for i32 profiles
/// * `P` - The number of SIMD lanes for i64 profiles
/// * `S` - The size of the alphabet (usually 5 for DNA including *N*)
#[derive(Debug, Clone)]
pub struct LocalProfiles<'a, const M: usize, const N: usize, const O: usize, const P: usize, const S: usize>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
    LaneCount<P>: SupportedLaneCount, {
    pub(crate) query:       &'a [u8],
    pub(crate) matrix:      &'a WeightMatrix<i8, S>,
    pub(crate) gap_open:    i8,
    pub(crate) gap_extend:  i8,
    pub(crate) profile_i8:  OnceCell<StripedProfile<i8, M, S>>,
    pub(crate) profile_i16: OnceCell<StripedProfile<i16, N, S>>,
    pub(crate) profile_i32: OnceCell<StripedProfile<i32, O, S>>,
    pub(crate) profile_i64: OnceCell<StripedProfile<i64, P, S>>,
}

impl<'a, const M: usize, const N: usize, const O: usize, const P: usize, const S: usize> LocalProfiles<'a, M, N, O, P, S>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
    LaneCount<P>: SupportedLaneCount,
{
    /// Creates an empty [`LocalProfiles`]. Usually you want to use
    /// [`LocalProfiles::new_with_w256`] or [`LocalProfiles::new_with_w512`]
    /// instead.
    ///
    /// ## Errors
    ///
    /// The following errors are possible:
    /// * [`QueryProfileError::EmptyQuery`] if `query` is empty
    /// * [`QueryProfileError::GapOpenOutOfRange`] if `gap_open` is not between
    ///   -127 and 0, inclusive
    /// * [`QueryProfileError::GapExtendOutOfRange`] if `gap_extend` is not
    ///   between -127 and 0, inclusive
    /// * [`QueryProfileError::BadGapWeights`] if `gap_extend` is less than
    ///   `gap_open`
    #[inline]
    pub fn new<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, QueryProfileError> {
        validate_profile_args(query.as_ref(), gap_open, gap_extend)?;

        Ok(LocalProfiles {
            query: query.as_ref(),
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceCell::new(),
            profile_i16: OnceCell::new(),
            profile_i32: OnceCell::new(),
            profile_i64: OnceCell::new(),
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i8` and `M` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i8(&self) -> &StripedProfile<i8, M, S> {
        // We already validated profile
        self.profile_i8.get_or_init(|| {
            StripedProfile::<i8, M, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i16` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i16(&self) -> &StripedProfile<i16, N, S> {
        // We already validated profile
        self.profile_i16.get_or_init(|| {
            StripedProfile::<i16, N, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i32` and `O` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i32(&self) -> &StripedProfile<i32, O, S> {
        // We already validated profile
        self.profile_i32.get_or_init(|| {
            StripedProfile::<i32, O, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i64` and `P` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i64(&self) -> &StripedProfile<i64, P, S> {
        // We already validated profile
        self.profile_i64.get_or_init(|| {
            StripedProfile::<i64, P, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `i8` profile. Lazily initializes the profiles and works its way up
    /// to the `i64` profile. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.smith_waterman_score_from_i8(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score_from_i8<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_i8()
            .smith_waterman_score(query)
            .or_else(|| self.get_i16().smith_waterman_score(query))
            .or_else(|| self.get_i32().smith_waterman_score(query))
            .or_else(|| self.get_i64().smith_waterman_score(query))
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `i16` profile, skipping the `i8` profile. Lazily initializes the
    /// profiles and works its way up to the `i64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.smith_waterman_score_from_i16(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score_from_i16<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_i16()
            .smith_waterman_score(query)
            .or_else(|| self.get_i32().smith_waterman_score(query))
            .or_else(|| self.get_i64().smith_waterman_score(query))
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `i32` profile, skipping the `i8` and `i16` profiles. Lazily
    /// initializes the profiles and works its way up to the `i64` profile.
    /// Execution stops when the score returned no longer overflows the
    /// profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = LocalProfiles::new_with_w256(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.smith_waterman_score_from_i32(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score_from_i32<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_i32()
            .smith_waterman_score(query)
            .or_else(|| self.get_i64().smith_waterman_score(query))
    }
}

impl<'a, const S: usize> LocalProfiles<'a, 16, 8, 4, 2, S> {
    /// Creates an empty [`LocalProfiles`] optimized for 128-bit SIMD width.
    /// Sets M=16, N=8, O=4, P=2 for i8, i16, i32, i64 profiles respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`LocalProfiles::new`]
    #[inline]
    pub fn new_with_w128<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<LocalProfiles<'a, 16, 8, 4, 2, S>, QueryProfileError> {
        validate_profile_args(query.as_ref(), gap_open, gap_extend)?;

        Ok(LocalProfiles {
            query: query.as_ref(),
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceCell::new(),
            profile_i16: OnceCell::new(),
            profile_i32: OnceCell::new(),
            profile_i64: OnceCell::new(),
        })
    }
}

impl<'a, const S: usize> LocalProfiles<'a, 32, 16, 8, 4, S> {
    /// Creates an empty [`LocalProfiles`] optimized for 256-bit SIMD width.
    /// Sets M=32, N=16, O=8, P=4 for i8, i16, i32, i64 profiles respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`LocalProfiles::new`]
    #[inline]
    pub fn new_with_w256<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<LocalProfiles<'a, 32, 16, 8, 4, S>, QueryProfileError> {
        validate_profile_args(query.as_ref(), gap_open, gap_extend)?;

        Ok(LocalProfiles {
            query: query.as_ref(),
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceCell::new(),
            profile_i16: OnceCell::new(),
            profile_i32: OnceCell::new(),
            profile_i64: OnceCell::new(),
        })
    }
}

impl<'a, const S: usize> LocalProfiles<'a, 64, 32, 16, 8, S> {
    /// Creates an empty [`LocalProfiles`] optimized for 512-bit SIMD width.
    /// Sets M=64, N=32, O=16, P=8 for i8, i16, i32, i64 profiles respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`LocalProfiles::new`]
    #[inline]
    pub fn new_with_w512<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<LocalProfiles<'a, 64, 32, 16, 8, S>, QueryProfileError> {
        validate_profile_args(query.as_ref(), gap_open, gap_extend)?;

        Ok(LocalProfiles {
            query: query.as_ref(),
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceCell::new(),
            profile_i16: OnceCell::new(),
            profile_i32: OnceCell::new(),
            profile_i64: OnceCell::new(),
        })
    }
}

/// A lazily-evaluated set of striped alignment profiles which can be shared
/// across threads. This is an abstraction around [`StripedProfile`], providing
/// convenience methods for automatically increasing the integer width and
/// rerunning the alignment when overflow occurs. This only supports the
/// unsigned version of the algorithm.
///
/// When sharing between threads is not needed, consider using [`LocalProfiles`]
/// instead.
///
/// ## Type Parameters
///
/// * `M` - The number of SIMD lanes for i8 profiles
/// * `N` - The number of SIMD lanes for i16 profiles
/// * `O` - The number of SIMD lanes for i32 profiles
/// * `P` - The number of SIMD lanes for i64 profiles
/// * `S` - The size of the alphabet (usually 5 for DNA including *N*)
#[derive(Debug, Clone)]
pub struct SharedProfiles<'a, const M: usize, const N: usize, const O: usize, const P: usize, const S: usize>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
    LaneCount<P>: SupportedLaneCount, {
    pub(crate) query:       Box<[u8]>,
    pub(crate) matrix:      &'a WeightMatrix<i8, S>,
    pub(crate) gap_open:    i8,
    pub(crate) gap_extend:  i8,
    pub(crate) profile_i8:  OnceLock<StripedProfile<i8, M, S>>,
    pub(crate) profile_i16: OnceLock<StripedProfile<i16, N, S>>,
    pub(crate) profile_i32: OnceLock<StripedProfile<i32, O, S>>,
    pub(crate) profile_i64: OnceLock<StripedProfile<i64, P, S>>,
}

impl<'a, const M: usize, const N: usize, const O: usize, const P: usize, const S: usize> SharedProfiles<'a, M, N, O, P, S>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
    LaneCount<P>: SupportedLaneCount,
{
    /// Creates an empty [`SharedProfiles`]. Usually you want to use
    /// [`SharedProfiles::new_with_w256`] or [`SharedProfiles::new_with_w512`]
    /// instead.
    ///
    /// ## Errors
    ///
    /// The following errors are possible:
    /// * [`QueryProfileError::EmptyQuery`] if `query` is empty
    /// * [`QueryProfileError::GapOpenOutOfRange`] if `gap_open` is not between
    ///   -127 and 0, inclusive
    /// * [`QueryProfileError::GapExtendOutOfRange`] if `gap_extend` is not
    ///   between -127 and 0, inclusive
    /// * [`QueryProfileError::BadGapWeights`] if `gap_extend` is less than
    ///   `gap_open`
    #[inline]
    pub fn new(
        query: Box<[u8]>, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, QueryProfileError> {
        validate_profile_args(&query, gap_open, gap_extend)?;

        Ok(SharedProfiles {
            query,
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceLock::new(),
            profile_i16: OnceLock::new(),
            profile_i32: OnceLock::new(),
            profile_i64: OnceLock::new(),
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i8` and `M` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i8(&self) -> &StripedProfile<i8, M, S> {
        // We already validated profile
        self.profile_i8.get_or_init(|| {
            StripedProfile::<i8, M, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i16` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i16(&self) -> &StripedProfile<i16, N, S> {
        self.profile_i16.get_or_init(|| {
            // We already validated profile
            StripedProfile::<i16, N, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i32` and `O` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i32(&self) -> &StripedProfile<i32, O, S> {
        self.profile_i32.get_or_init(|| {
            // We already validated profile
            StripedProfile::<i32, O, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i64` and `P` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i64(&self) -> &StripedProfile<i64, P, S> {
        self.profile_i64.get_or_init(|| {
            // We already validated profile
            StripedProfile::<i64, P, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `i8` profile. Lazily initializes the profiles and works its way up
    /// to the `i64` profile. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = SharedProfiles::new_with_w256(query.into(), &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.smith_waterman_score_from_i8(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score_from_i8<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_i8()
            .smith_waterman_score(query)
            .or_else(|| self.get_i16().smith_waterman_score(query))
            .or_else(|| self.get_i32().smith_waterman_score(query))
            .or_else(|| self.get_i64().smith_waterman_score(query))
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `i16` profile, skipping the `i8` profile. Lazily initializes the
    /// profiles and works its way up to the `i64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = SharedProfiles::new_with_w256(query.into(), &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.smith_waterman_score_from_i16(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score_from_i16<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_i16()
            .smith_waterman_score(query)
            .or_else(|| self.get_i32().smith_waterman_score(query))
            .or_else(|| self.get_i64().smith_waterman_score(query))
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `i32` profile, skipping the `i8` and `i16` profiles. Lazily
    /// initializes the profiles and works its way up to the `i64` profile.
    /// Execution stops when the score returned no longer overflows the
    /// profile's integer range.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let profile = SharedProfiles::new_with_w256(query.into(), &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let score = profile.smith_waterman_score_from_i32(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[inline]
    #[must_use]
    pub fn smith_waterman_score_from_i32<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_i32()
            .smith_waterman_score(query)
            .or_else(|| self.get_i64().smith_waterman_score(query))
    }
}

impl<'a, const S: usize> SharedProfiles<'a, 16, 8, 4, 2, S> {
    /// Creates an empty [`SharedProfiles`] optimized for 128-bit SIMD width.
    /// Sets M=16, N=8, O=4, P=2 for i8, i16, i32, i64 profiles respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`SharedProfiles::new`]
    #[inline]
    pub fn new_with_w128(
        query: Box<[u8]>, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<SharedProfiles<'a, 16, 8, 4, 2, S>, QueryProfileError> {
        validate_profile_args(&query, gap_open, gap_extend)?;

        Ok(SharedProfiles {
            query,
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceLock::new(),
            profile_i16: OnceLock::new(),
            profile_i32: OnceLock::new(),
            profile_i64: OnceLock::new(),
        })
    }
}

impl<'a, const S: usize> SharedProfiles<'a, 32, 16, 8, 4, S> {
    /// Creates an empty [`SharedProfiles`] optimized for 256-bit SIMD width.
    /// Sets M=32, N=16, O=8, P=4 for i8, i16, i32, i64 profiles respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`SharedProfiles::new`]
    #[inline]
    pub fn new_with_w256(
        query: Box<[u8]>, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<SharedProfiles<'a, 32, 16, 8, 4, S>, QueryProfileError> {
        validate_profile_args(&query, gap_open, gap_extend)?;

        Ok(SharedProfiles {
            query,
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceLock::new(),
            profile_i16: OnceLock::new(),
            profile_i32: OnceLock::new(),
            profile_i64: OnceLock::new(),
        })
    }
}

impl<'a, const S: usize> SharedProfiles<'a, 64, 32, 16, 8, S> {
    /// Creates an empty [`SharedProfiles`] optimized for 512-bit SIMD width.
    /// Sets M=64, N=32, O=16, P=8 for i8, i16, i32, i64 profiles respectively.
    ///
    /// ## Errors
    ///
    /// Same as [`SharedProfiles::new`]
    #[inline]
    pub fn new_with_w512(
        query: Box<[u8]>, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<SharedProfiles<'a, 64, 32, 16, 8, S>, QueryProfileError> {
        validate_profile_args(&query, gap_open, gap_extend)?;

        Ok(SharedProfiles {
            query,
            matrix,
            gap_open,
            gap_extend,
            profile_i8: OnceLock::new(),
            profile_i16: OnceLock::new(),
            profile_i32: OnceLock::new(),
            profile_i64: OnceLock::new(),
        })
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
