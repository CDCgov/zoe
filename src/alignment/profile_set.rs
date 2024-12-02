use super::{validate_profile_args, StripedProfile};
use crate::data::{err::AlignmentError, matrices::BiasedWeightMatrix, types::nucleotides::MaybeNucleic};
use std::{
    cell::OnceCell,
    simd::{LaneCount, SupportedLaneCount},
    sync::OnceLock,
};

/// Creates a lazy set of striped profiles for local (thread-specific) use.
/// The number of SIMD lanes `N` must be specified.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalProfile<'a, const N: usize, const S: usize>
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

impl<'a, const N: usize, const S: usize> LocalProfile<'a, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates an empty [`LocalProfile`]. Usually you want to use
    /// [`LocalProfile::new_with_u8`] or [`LocalProfile::new_with_u16`] instead.
    ///
    /// # Errors
    ///
    /// Will return [`AlignmentError::EmptyQuery`] if `query` is empty or
    /// [`AlignmentError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        query: &'a T, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, AlignmentError> {
        validate_profile_args(query.as_ref(), gap_open, gap_extend)?;

        Ok(LocalProfile {
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

    /// Creates an empty [`LocalProfile`] and eagerly initializes the
    /// `u8` profile.
    ///
    /// # Errors
    ///
    /// Will return [`AlignmentError::EmptyQuery`] if `query` is empty or
    /// [`AlignmentError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_u8<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        query: &'a T, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, AlignmentError> {
        let p = LocalProfile::new(query, matrix, gap_open, gap_extend)?;
        p.get_u8();
        Ok(p)
    }

    /// Creates an empty [`LocalProfile`] and eagerly initializes the
    /// `u16` profile.
    ///
    /// # Errors
    ///
    /// Will return [`AlignmentError::EmptyQuery`] if `query` is empty or
    /// [`AlignmentError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_u16<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        query: &'a T, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, AlignmentError> {
        let p = LocalProfile::new(query, matrix, gap_open, gap_extend)?;
        p.get_u16();
        Ok(p)
    }

    /// Get or initialize [`StripedProfile`] with elements of `u8` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub fn get_u8(&self) -> &StripedProfile<u8, N, S> {
        // Unwrap will not panic since we already validated profile
        self.profile_u8.get_or_init(|| {
            StripedProfile::<u8, N, S>::new(self.query, self.matrix, self.gap_open, self.gap_extend).unwrap()
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u16` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub fn get_u16(&self) -> &StripedProfile<u16, N, S> {
        // Unwrap will not panic since we already validated profile
        self.profile_u16.get_or_init(|| {
            StripedProfile::<u16, N, S>::new(self.query, self.matrix, u16::from(self.gap_open), u16::from(self.gap_extend))
                .unwrap()
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u32` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub fn get_u32(&self) -> &StripedProfile<u32, N, S> {
        // Unwrap will not panic since we already validated profile
        self.profile_u32.get_or_init(|| {
            StripedProfile::<u32, N, S>::new(self.query, self.matrix, u32::from(self.gap_open), u32::from(self.gap_extend))
                .unwrap()
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u64` and `N` SIMD
    /// lanes and returns a reference to the field.
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub fn get_u64(&self) -> &StripedProfile<u64, N, S> {
        // Unwrap will not panic since we already validated profile
        self.profile_u64.get_or_init(|| {
            StripedProfile::<u64, N, S>::new(self.query, self.matrix, u64::from(self.gap_open), u64::from(self.gap_extend))
                .unwrap()
        })
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `u8` profile. Lazily initializes the profiles and works its way up
    /// to the `u64` profile. Execution stops when the score returned no longer
    /// overflows the profile's integer range.
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u8<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(&self, query: &T) -> Option<u64> {
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
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u16<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(&self, query: &T) -> Option<u64> {
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
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u32<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(&self, query: &T) -> Option<u64> {
        let query = query.as_ref();
        self.get_u32()
            .smith_waterman_score(query)
            .or_else(|| self.get_u64().smith_waterman_score(query))
    }
}

/// Creates a thread-safe, lazy set of striped profiles that can be shared
/// across threads. The number of SIMD lanes `N` must be specified.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SharedProfile<'a, const N: usize, const S: usize>
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

impl<'a, const N: usize, const S: usize> SharedProfile<'a, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates an empty [`SharedProfile`]. Usually you want to use
    /// [`SharedProfile::new_with_u8`] or [`SharedProfile::new_with_u16`]
    /// instead.
    ///
    /// # Errors
    ///
    /// Will return [`AlignmentError::EmptyQuery`] if `query` is empty or
    /// [`AlignmentError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new(
        query: Box<[u8]>, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, AlignmentError> {
        validate_profile_args(&query, gap_open, gap_extend)?;

        Ok(SharedProfile {
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

    /// Creates an empty [`SharedProfile`] and eagerly initializes the `u8`
    /// profile.
    ///
    /// # Errors
    ///
    /// Will return [`AlignmentError::EmptyQuery`] if `query` is empty or
    /// [`AlignmentError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_u8(
        query: Box<[u8]>, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, AlignmentError> {
        let p = SharedProfile::new(query, matrix, gap_open, gap_extend)?;
        p.get_u8();
        Ok(p)
    }

    /// Creates an empty [`SharedProfile`] and eagerly initializes the
    /// `u16` profile.
    ///
    /// # Errors
    ///
    /// Will return [`AlignmentError::EmptyQuery`] if `query` is empty or
    /// [`AlignmentError::BadGapWeights`] if `gap_extend` is greater than
    /// `gap_open`.
    #[inline]
    pub fn new_with_u16(
        query: Box<[u8]>, matrix: &'a BiasedWeightMatrix<S>, gap_open: u8, gap_extend: u8,
    ) -> Result<Self, AlignmentError> {
        let p = SharedProfile::new(query, matrix, gap_open, gap_extend)?;
        p.get_u16();
        Ok(p)
    }

    /// Get or initialize [`StripedProfile`] with elements of `u8` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub fn get_u8(&self) -> &StripedProfile<u8, N, S> {
        // Unwrap will not panic since we already validated profile
        self.profile_u8.get_or_init(|| {
            StripedProfile::<u8, N, S>::new(&self.query, self.matrix, self.gap_open, self.gap_extend).unwrap()
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u16` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub fn get_u16(&self) -> &StripedProfile<u16, N, S> {
        self.profile_u16.get_or_init(|| {
            // Unwrap will not panic since we already validated profile
            StripedProfile::<u16, N, S>::new(&self.query, self.matrix, u16::from(self.gap_open), u16::from(self.gap_extend))
                .unwrap()
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u32` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub fn get_u32(&self) -> &StripedProfile<u32, N, S> {
        self.profile_u32.get_or_init(|| {
            // Unwrap will not panic since we already validated profile
            StripedProfile::<u32, N, S>::new(&self.query, self.matrix, u32::from(self.gap_open), u32::from(self.gap_extend))
                .unwrap()
        })
    }

    /// Get or initialize [`StripedProfile`] with elements of `u64` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    #[allow(clippy::missing_panics_doc)]
    pub fn get_u64(&self) -> &StripedProfile<u64, N, S> {
        self.profile_u64.get_or_init(|| {
            // Unwrap will not panic since we already validated profile
            StripedProfile::<u64, N, S>::new(&self.query, self.matrix, u64::from(self.gap_open), u64::from(self.gap_extend))
                .unwrap()
        })
    }

    /// Lazily execute [`StripedProfile::smith_waterman_score`] starting with
    /// the `u8` profile, skipping the `u8` profile. Lazily initializes the
    /// profiles and works its way up to the `u64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u8<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(&self, query: &T) -> Option<u64> {
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
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u16<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(&self, query: &T) -> Option<u64> {
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
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u32<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(&self, query: &T) -> Option<u64> {
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
        let profiles = LocalProfile::<32, 5>::new(&v, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        let profile1 = profiles.get_u8();
        let profile2 = StripedProfile::<u8, 32, 5>::new(v, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        assert_eq!(profile1, &profile2);
    }
}
