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
/// and rerunning the alignment when overflow occurs. The number of SIMD lanes
/// `N` must be specified.
///
/// If it is necessary to share between multiple threads, consider using
/// [`SharedProfiles`].
#[derive(Debug, Clone)]
pub struct LocalProfiles<'a, const N: usize, const S: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) query:       &'a [u8],
    pub(crate) matrix:      &'a WeightMatrix<i8, S>,
    pub(crate) gap_open:    i8,
    pub(crate) gap_extend:  i8,
    pub(crate) profile_i8:  OnceCell<StripedProfile<i8, N, S>>,
    pub(crate) profile_i16: OnceCell<StripedProfile<i16, N, S>>,
    pub(crate) profile_i32: OnceCell<StripedProfile<i32, N, S>>,
    pub(crate) profile_i64: OnceCell<StripedProfile<i64, N, S>>,
}

impl<'a, const N: usize, const S: usize> LocalProfiles<'a, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates an empty [`LocalProfiles`]. Usually you want to use
    /// [`LocalProfiles::new_with_i8`] or [`LocalProfiles::new_with_i16`] instead.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
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

    /// Creates an empty [`LocalProfiles`] and eagerly initializes the `i8`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_i8<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, QueryProfileError> {
        let p = LocalProfiles::new(query, matrix, gap_open, gap_extend)?;
        p.get_i8();
        Ok(p)
    }

    /// Creates an empty [`LocalProfiles`] and eagerly initializes the `i16`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_i16<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, QueryProfileError> {
        let p = LocalProfiles::new(query, matrix, gap_open, gap_extend)?;
        p.get_i16();
        Ok(p)
    }

    /// Get or initialize [`StripedProfile`] with elements of `i8` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i8(&self) -> &StripedProfile<i8, N, S> {
        // We already validated profile
        self.profile_i8.get_or_init(|| {
            StripedProfile::<i8, N, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
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

    /// Get or initialize [`StripedProfile`] with elements of `i32` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i32(&self) -> &StripedProfile<i32, N, S> {
        // We already validated profile
        self.profile_i32.get_or_init(|| {
            StripedProfile::<i32, N, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i64` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i64(&self) -> &StripedProfile<i64, N, S> {
        // We already validated profile
        self.profile_i64.get_or_init(|| {
            StripedProfile::<i64, N, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `i8` profile. Lazily initializes the profiles and works its way up
    /// to the `i64` profile. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = LocalProfiles::<32, 5>::new_with_i8(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_i8(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
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
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = LocalProfiles::<32, 5>::new_with_i16(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_i16(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
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
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = LocalProfiles::<32, 5>::new(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_i32(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_i32<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_i32()
            .smith_waterman_score(query)
            .or_else(|| self.get_i64().smith_waterman_score(query))
    }
}

/// A lazily-evaluated set of striped alignment profiles which can be shared
/// across threads. This is an abstraction around [`StripedProfile`], providing
/// convenience methods for automatically increasing the integer width and
/// rerunning the alignment when overflow occurs. The number of SIMD lanes `N`
/// must be specified.
///
/// When sharing between threads is not needed, consider using [`LocalProfiles`]
/// instead.
#[derive(Debug, Clone)]
pub struct SharedProfiles<'a, const N: usize, const S: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) query:       Box<[u8]>,
    pub(crate) matrix:      &'a WeightMatrix<i8, S>,
    pub(crate) gap_open:    i8,
    pub(crate) gap_extend:  i8,
    pub(crate) profile_i8:  OnceLock<StripedProfile<i8, N, S>>,
    pub(crate) profile_i16: OnceLock<StripedProfile<i16, N, S>>,
    pub(crate) profile_i32: OnceLock<StripedProfile<i32, N, S>>,
    pub(crate) profile_i64: OnceLock<StripedProfile<i64, N, S>>,
}

impl<'a, const N: usize, const S: usize> SharedProfiles<'a, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates an empty [`SharedProfiles`]. Usually you want to use
    /// [`SharedProfiles::new_with_i8`] or [`SharedProfiles::new_with_i16`]
    /// instead.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
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

    /// Creates an empty [`SharedProfiles`] and eagerly initializes the `i8`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_i8(
        query: Box<[u8]>, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, QueryProfileError> {
        let p = SharedProfiles::new(query, matrix, gap_open, gap_extend)?;
        p.get_i8();
        Ok(p)
    }

    /// Creates an empty [`SharedProfiles`] and eagerly initializes the `i16`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_i16(
        query: Box<[u8]>, matrix: &'a WeightMatrix<i8, S>, gap_open: i8, gap_extend: i8,
    ) -> Result<Self, QueryProfileError> {
        let p = SharedProfiles::new(query, matrix, gap_open, gap_extend)?;
        p.get_i16();
        Ok(p)
    }

    /// Get or initialize [`StripedProfile`] with elements of `i8` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i8(&self) -> &StripedProfile<i8, N, S> {
        // We already validated profile
        self.profile_i8.get_or_init(|| {
            StripedProfile::<i8, N, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
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

    /// Get or initialize [`StripedProfile`] with elements of `i32` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i32(&self) -> &StripedProfile<i32, N, S> {
        self.profile_i32.get_or_init(|| {
            // We already validated profile
            StripedProfile::<i32, N, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `i64` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_i64(&self) -> &StripedProfile<i64, N, S> {
        self.profile_i64.get_or_init(|| {
            // We already validated profile
            StripedProfile::<i64, N, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `i8` profile. Lazily initializes the profiles and works its way up
    /// to the `i64` profile. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = SharedProfiles::<32, 5>::new(query.into(), &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_i8(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
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
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = SharedProfiles::<32, 5>::new(query.into(), &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_i16(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
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
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::WeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = SharedProfiles::<32, 5>::new(query.into(), &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_i32(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_i32<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_i32()
            .smith_waterman_score(query)
            .or_else(|| self.get_i64().smith_waterman_score(query))
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
        let profiles = LocalProfiles::<32, 5>::new(&v, &WEIGHTS, GAP_OPEN as i8, GAP_EXTEND as i8).unwrap();
        let profile1 = profiles.get_i8();
        let profile2 = StripedProfile::<i8, 32, 5>::new(v, &WEIGHTS, GAP_OPEN as i8, GAP_EXTEND as i8).unwrap();
        assert_eq!(profile1, &profile2);
    }
}
