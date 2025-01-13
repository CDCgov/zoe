use super::{StripedProfile, validate_profile_args};
use crate::data::{err::QueryProfileError, matrices::BiasedWeightMatrix};
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
    pub(crate) matrix:      &'a BiasedWeightMatrix<S>,
    pub(crate) gap_open:    u8,
    pub(crate) gap_extend:  u8,
    pub(crate) profile_u8:  OnceCell<StripedProfile<u8, N, S>>,
    pub(crate) profile_u16: OnceCell<StripedProfile<u16, N, S>>,
    pub(crate) profile_u32: OnceCell<StripedProfile<u32, N, S>>,
    pub(crate) profile_u64: OnceCell<StripedProfile<u64, N, S>>,
}

impl<'a, const N: usize, const S: usize> LocalProfiles<'a, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates an empty [`LocalProfiles`]. Usually you want to use
    /// [`LocalProfiles::new_with_u8`] or [`LocalProfiles::new_with_u16`] instead.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, QueryProfileError> {
        validate_profile_args(query.as_ref(), gap_open, gap_extend)?;

        Ok(LocalProfiles {
            query: query.as_ref(),
            matrix,
            gap_open,
            gap_extend,
            profile_u8: OnceCell::new(),
            profile_u16: OnceCell::new(),
            profile_u32: OnceCell::new(),
            profile_u64: OnceCell::new(),
        })
    }

    /// Creates an empty [`LocalProfiles`] and eagerly initializes the `u8`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_u8<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, QueryProfileError> {
        let p = LocalProfiles::new(query, matrix, gap_open, gap_extend)?;
        p.get_u8();
        Ok(p)
    }

    /// Creates an empty [`LocalProfiles`] and eagerly initializes the `u16`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_u16<T: AsRef<[u8]> + ?Sized>(
        query: &'a T, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, QueryProfileError> {
        let p = LocalProfiles::new(query, matrix, gap_open, gap_extend)?;
        p.get_u16();
        Ok(p)
    }

    /// Get or initialize [`StripedProfile`] with elements of `u8` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_u8(&self) -> &StripedProfile<u8, N, S> {
        // We already validated profile
        self.profile_u8.get_or_init(|| {
            StripedProfile::<u8, N, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u16` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_u16(&self) -> &StripedProfile<u16, N, S> {
        // We already validated profile
        self.profile_u16.get_or_init(|| {
            StripedProfile::<u16, N, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u32` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_u32(&self) -> &StripedProfile<u32, N, S> {
        // We already validated profile
        self.profile_u32.get_or_init(|| {
            StripedProfile::<u32, N, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u64` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_u64(&self) -> &StripedProfile<u64, N, S> {
        // We already validated profile
        self.profile_u64.get_or_init(|| {
            StripedProfile::<u64, N, S>::new_unchecked(self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `u8` profile. Lazily initializes the profiles and works its way up
    /// to the `u64` profile. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::BiasedWeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: BiasedWeightMatrix<5> = BiasedWeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = LocalProfiles::<32, 5>::new_with_u8(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_u8(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u8<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_u8()
            .smith_waterman_score(query)
            .or_else(|| self.get_u16().smith_waterman_score(query))
            .or_else(|| self.get_u32().smith_waterman_score(query))
            .or_else(|| self.get_u64().smith_waterman_score(query))
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `u16` profile, skipping the `u8` profile. Lazily initializes the
    /// profiles and works its way up to the `u64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::BiasedWeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: BiasedWeightMatrix<5> = BiasedWeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = LocalProfiles::<32, 5>::new_with_u16(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_u16(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u16<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_u16()
            .smith_waterman_score(query)
            .or_else(|| self.get_u32().smith_waterman_score(query))
            .or_else(|| self.get_u64().smith_waterman_score(query))
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `u32` profile, skipping the `u8` and `u16` profiles. Lazily
    /// initializes the profiles and works its way up to the `u64` profile.
    /// Execution stops when the score returned no longer overflows the
    /// profile's integer range.
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{LocalProfiles, sw::sw_simd_score}, data::BiasedWeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: BiasedWeightMatrix<5> = BiasedWeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = LocalProfiles::<32, 5>::new(query, &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_u32(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u32<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_u32()
            .smith_waterman_score(query)
            .or_else(|| self.get_u64().smith_waterman_score(query))
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
    pub(crate) matrix:      &'a BiasedWeightMatrix<S>,
    pub(crate) gap_open:    u8,
    pub(crate) gap_extend:  u8,
    pub(crate) profile_u8:  OnceLock<StripedProfile<u8, N, S>>,
    pub(crate) profile_u16: OnceLock<StripedProfile<u16, N, S>>,
    pub(crate) profile_u32: OnceLock<StripedProfile<u32, N, S>>,
    pub(crate) profile_u64: OnceLock<StripedProfile<u64, N, S>>,
}

impl<'a, const N: usize, const S: usize> SharedProfiles<'a, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates an empty [`SharedProfiles`]. Usually you want to use
    /// [`SharedProfiles::new_with_u8`] or [`SharedProfiles::new_with_u16`]
    /// instead.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new(
        query: Box<[u8]>, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, QueryProfileError> {
        validate_profile_args(&query, gap_open, gap_extend)?;

        Ok(SharedProfiles {
            query,
            matrix,
            gap_open,
            gap_extend,
            profile_u8: OnceLock::new(),
            profile_u16: OnceLock::new(),
            profile_u32: OnceLock::new(),
            profile_u64: OnceLock::new(),
        })
    }

    /// Creates an empty [`SharedProfiles`] and eagerly initializes the `u8`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_u8(
        query: Box<[u8]>, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, QueryProfileError> {
        let p = SharedProfiles::new(query, matrix, gap_open, gap_extend)?;
        p.get_u8();
        Ok(p)
    }

    /// Creates an empty [`SharedProfiles`] and eagerly initializes the `u16`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`QueryProfileError::EmptyQuery`] if `query` is empty or
    /// [`QueryProfileError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_u16(
        query: Box<[u8]>, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, QueryProfileError> {
        let p = SharedProfiles::new(query, matrix, gap_open, gap_extend)?;
        p.get_u16();
        Ok(p)
    }

    /// Get or initialize [`StripedProfile`] with elements of `u8` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_u8(&self) -> &StripedProfile<u8, N, S> {
        // We already validated profile
        self.profile_u8.get_or_init(|| {
            StripedProfile::<u8, N, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u16` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_u16(&self) -> &StripedProfile<u16, N, S> {
        self.profile_u16.get_or_init(|| {
            // We already validated profile
            StripedProfile::<u16, N, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u32` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_u32(&self) -> &StripedProfile<u32, N, S> {
        self.profile_u32.get_or_init(|| {
            // We already validated profile
            StripedProfile::<u32, N, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u64` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    pub fn get_u64(&self) -> &StripedProfile<u64, N, S> {
        self.profile_u64.get_or_init(|| {
            // We already validated profile
            StripedProfile::<u64, N, S>::new_unchecked(&self.query, self.matrix, self.gap_open, self.gap_extend)
        })
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `u8` profile. Lazily initializes the profiles and works its way up
    /// to the `u64` profile. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::BiasedWeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: BiasedWeightMatrix<5> = BiasedWeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = SharedProfiles::<32, 5>::new(query.into(), &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_u8(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u8<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_u8()
            .smith_waterman_score(query)
            .or_else(|| self.get_u16().smith_waterman_score(query))
            .or_else(|| self.get_u32().smith_waterman_score(query))
            .or_else(|| self.get_u64().smith_waterman_score(query))
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `u16` profile, skipping the `u8` profile. Lazily initializes the
    /// profiles and works its way up to the `u64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::BiasedWeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: BiasedWeightMatrix<5> = BiasedWeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = SharedProfiles::<32, 5>::new(query.into(), &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_u16(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u16<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_u16()
            .smith_waterman_score(query)
            .or_else(|| self.get_u32().smith_waterman_score(query))
            .or_else(|| self.get_u64().smith_waterman_score(query))
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `u32` profile, skipping the `u8` and `u16` profiles. Lazily
    /// initializes the profiles and works its way up to the `u64` profile.
    /// Execution stops when the score returned no longer overflows the
    /// profile's integer range.
    ///
    /// ### Example
    ///
    /// ```
    /// # use zoe::{alignment::{SharedProfiles, sw::sw_simd_score}, data::BiasedWeightMatrix};
    /// let reference: &[u8] = b"ATGCATCGATCGATCGATCGATCGATCGATGC";
    /// let query: &[u8] = b"CGTTCGCCATAAAGGGGG";
    ///
    /// const WEIGHTS: BiasedWeightMatrix<5> = BiasedWeightMatrix::new_biased_dna_matrix(4, -2, Some(b'N'));
    ///
    /// let profile = SharedProfiles::<32, 5>::new(query.into(), &WEIGHTS, 3, 1).unwrap();
    /// let score = profile.smith_waterman_score_from_u32(reference).unwrap();
    /// assert_eq!(score, 26);
    /// ```
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u32<T: AsRef<[u8]> + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_u32()
            .smith_waterman_score(query)
            .or_else(|| self.get_u64().smith_waterman_score(query))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::alignment::sw::test_data::{BIASED_WEIGHTS, GAP_EXTEND, GAP_OPEN};

    #[test]
    fn sw_simd_profile_set() {
        let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
        let profiles = LocalProfiles::<32, 5>::new(&v, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        let profile1 = profiles.get_u8();
        let profile2 = StripedProfile::<u8, 32, 5>::new(v, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        assert_eq!(profile1, &profile2);
    }
}
