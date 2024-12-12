use super::StripedDNAProfile;
use crate::data::{matrices::BiasedWeightMatrix, types::nucleotides::MaybeNucleic};
use std::{
    cell::OnceCell,
    simd::{LaneCount, SupportedLaneCount},
    sync::OnceLock,
};

/// Creates a lazy set of striped DNA profiles for local (thread-specific) use.
/// The number of SIMD lanes `N` must be specified.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalDNAProfile<'a, const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) profile_seq: &'a [u8],
    pub(crate) matrix:      &'a BiasedWeightMatrix<5>,
    pub(crate) profile_u8:  OnceCell<StripedDNAProfile<u8, N>>,
    pub(crate) profile_u16: OnceCell<StripedDNAProfile<u16, N>>,
    pub(crate) profile_u32: OnceCell<StripedDNAProfile<u32, N>>,
    pub(crate) profile_u64: OnceCell<StripedDNAProfile<u64, N>>,
}

impl<'a, const N: usize> Default for LocalDNAProfile<'a, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn default() -> Self {
        Self {
            profile_seq: &[],
            matrix:      &BiasedWeightMatrix {
                index:   &[0; 5],
                mapping: [[0; 5]; 5],
                bias:    0,
            },
            profile_u8:  OnceCell::new(),
            profile_u16: OnceCell::new(),
            profile_u32: OnceCell::new(),
            profile_u64: OnceCell::new(),
        }
    }
}

impl<'a, const N: usize> LocalDNAProfile<'a, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates an empty [`LocalDNAProfile`]. Usually you want to use
    /// [`LocalDNAProfile::new_with_u8`] or
    /// [`LocalDNAProfile::new_with_u16`] instead.
    #[inline]
    #[must_use]
    pub fn new<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(seq: &'a T, matrix: &'a BiasedWeightMatrix<5>) -> Self {
        LocalDNAProfile {
            profile_seq: seq.as_ref(),
            matrix,
            ..Default::default()
        }
    }

    /// Creates an empty [`LocalDNAProfile`] and eagerly initializes the
    /// `u8` profile.
    #[inline]
    #[must_use]
    pub fn new_with_u8<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(seq: &'a T, matrix: &'a BiasedWeightMatrix<5>) -> Self {
        let p = LocalDNAProfile {
            profile_seq: seq.as_ref(),
            matrix,
            ..Default::default()
        };

        p.get_u8();
        p
    }

    /// Creates an empty [`LocalDNAProfile`] and eagerly initializes the
    /// `u16` profile.
    #[inline]
    #[must_use]
    pub fn new_with_u16<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(seq: &'a T, matrix: &'a BiasedWeightMatrix<5>) -> Self {
        let p = LocalDNAProfile {
            profile_seq: seq.as_ref(),
            matrix,
            ..Default::default()
        };
        p.get_u16();
        p
    }

    /// Get or initialize [`StripedDNAProfile`] with elements of `u8` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    pub fn get_u8(&self) -> &StripedDNAProfile<u8, N> {
        self.profile_u8
            .get_or_init(|| StripedDNAProfile::<u8, N>::new(self.profile_seq, self.matrix))
    }

    /// Get or initialize [`StripedDNAProfile`] with elements of `u16` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    pub fn get_u16(&self) -> &StripedDNAProfile<u16, N> {
        self.profile_u16
            .get_or_init(|| StripedDNAProfile::<u16, N>::new(self.profile_seq, self.matrix))
    }

    /// Get or initialize [`StripedDNAProfile`] with elements of `u32` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    pub fn get_u32(&self) -> &StripedDNAProfile<u32, N> {
        self.profile_u32
            .get_or_init(|| StripedDNAProfile::<u32, N>::new(self.profile_seq, self.matrix))
    }

    /// Get or initialize [`StripedDNAProfile`] with elements of `u64` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    pub fn get_u64(&self) -> &StripedDNAProfile<u64, N> {
        self.profile_u64
            .get_or_init(|| StripedDNAProfile::<u64, N>::new(self.profile_seq, self.matrix))
    }

    /// Lazily execute [`StripedDNAProfile::smith_waterman_score`] starting with
    /// the `u8` profile, skipping the `u8` profile. Lazily initializes the
    /// profiles and works its way up to the `u64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u8<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        &self, query: &T, gap_open: u8, gap_extend: u8,
    ) -> Option<u64> {
        let query = query.as_ref();
        self.get_u8()
            .smith_waterman_score(query, gap_open, gap_extend)
            .or_else(|| self.get_u16().smith_waterman_score(query, gap_open, gap_extend))
            .or_else(|| self.get_u32().smith_waterman_score(query, gap_open, gap_extend))
            .or_else(|| self.get_u64().smith_waterman_score(query, gap_open, gap_extend))
    }

    /// Lazily execute [`StripedDNAProfile::smith_waterman_score`] starting with
    /// the `u16` profile, skipping the `u8` profile. Lazily initializes the
    /// profiles and works its way up to the `u64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u16<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        &self, query: &T, gap_open: u8, gap_extend: u8,
    ) -> Option<u64> {
        let query = query.as_ref();
        self.get_u16()
            .smith_waterman_score(query, gap_open, gap_extend)
            .or_else(|| self.get_u32().smith_waterman_score(query, gap_open, gap_extend))
            .or_else(|| self.get_u64().smith_waterman_score(query, gap_open, gap_extend))
    }

    /// Lazily execute [`StripedDNAProfile::smith_waterman_score`] starting with
    /// the `u32` profile, skipping the `u8` and `u16` profiles. Lazily
    /// initializes the profiles and works its way up to the `u64` profile.
    /// Execution stops when the score returned no longer overflows the
    /// profile's integer range.
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u32<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        &self, query: &T, gap_open: u8, gap_extend: u8,
    ) -> Option<u64> {
        let query = query.as_ref();
        self.get_u32()
            .smith_waterman_score(query, gap_open, gap_extend)
            .or_else(|| self.get_u64().smith_waterman_score(query, gap_open, gap_extend))
    }
}

/// Creates a thread-safe, lazy set of striped DNA profiles that can be shared
/// across threads. The number of SIMD lanes `N` must be specified.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SharedDNAProfile<'a, const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    pub(crate) seq:         Box<[u8]>,
    pub(crate) matrix:      &'a BiasedWeightMatrix<5>,
    pub(crate) profile_u8:  OnceLock<StripedDNAProfile<u8, N>>,
    pub(crate) profile_u16: OnceLock<StripedDNAProfile<u16, N>>,
    pub(crate) profile_u32: OnceLock<StripedDNAProfile<u32, N>>,
    pub(crate) profile_u64: OnceLock<StripedDNAProfile<u64, N>>,
}

impl<'a, const N: usize> Default for SharedDNAProfile<'a, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn default() -> Self {
        Self {
            seq:         Box::new([]),
            matrix:      &BiasedWeightMatrix {
                index:   &[0; 5],
                mapping: [[0; 5]; 5],
                bias:    0,
            },
            profile_u8:  OnceLock::new(),
            profile_u16: OnceLock::new(),
            profile_u32: OnceLock::new(),
            profile_u64: OnceLock::new(),
        }
    }
}

impl<'a, const N: usize> SharedDNAProfile<'a, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Creates an empty [`SharedDNAProfile`]. Usually you want to use
    /// [`SharedDNAProfile::new_with_u8`] or
    /// [`SharedDNAProfile::new_with_u16`] instead.
    #[inline]
    #[must_use]
    pub fn new(seq: Box<[u8]>, matrix: &'a BiasedWeightMatrix<5>) -> Self {
        SharedDNAProfile {
            seq,
            matrix,
            ..Default::default()
        }
    }

    /// Creates an empty [`SharedDNAProfile`] and eagerly initializes the
    /// `u8` profile.
    #[inline]
    #[must_use]
    pub fn new_with_u8(seq: Box<[u8]>, matrix: &'a BiasedWeightMatrix<5>) -> Self {
        let p = SharedDNAProfile {
            seq,
            matrix,
            ..Default::default()
        };

        p.get_u8();
        p
    }

    /// Creates an empty [`SharedDNAProfile`] and eagerly initializes the
    /// `u16` profile.
    #[inline]
    #[must_use]
    pub fn new_with_u16(seq: Box<[u8]>, matrix: &'a BiasedWeightMatrix<5>) -> Self {
        let p = SharedDNAProfile {
            seq,
            matrix,
            ..Default::default()
        };
        p.get_u16();
        p
    }

    /// Get or initialize [`StripedDNAProfile`] with elements of `u8` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    pub fn get_u8(&self) -> &StripedDNAProfile<u8, N> {
        self.profile_u8
            .get_or_init(|| StripedDNAProfile::<u8, N>::new(&self.seq, self.matrix))
    }

    /// Get or initialize [`StripedDNAProfile`] with elements of `u16` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    pub fn get_u16(&self) -> &StripedDNAProfile<u16, N> {
        self.profile_u16
            .get_or_init(|| StripedDNAProfile::<u16, N>::new(&self.seq, self.matrix))
    }

    /// Get or initialize [`StripedDNAProfile`] with elements of `u32` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    pub fn get_u32(&self) -> &StripedDNAProfile<u32, N> {
        self.profile_u32
            .get_or_init(|| StripedDNAProfile::<u32, N>::new(&self.seq, self.matrix))
    }

    /// Get or initialize [`StripedDNAProfile`] with elements of `u64` and `N`
    /// SIMD lanes and returns a reference to the field.
    #[inline]
    pub fn get_u64(&self) -> &StripedDNAProfile<u64, N> {
        self.profile_u64
            .get_or_init(|| StripedDNAProfile::<u64, N>::new(&self.seq, self.matrix))
    }

    /// Lazily execute [`StripedDNAProfile::smith_waterman_score`] starting with
    /// the `u8` profile, skipping the `u8` profile. Lazily initializes the
    /// profiles and works its way up to the `u64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u8<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        &self, query: &T, gap_open: u8, gap_extend: u8,
    ) -> Option<u64> {
        let query = query.as_ref();
        self.get_u8()
            .smith_waterman_score(query, gap_open, gap_extend)
            .or_else(|| self.get_u16().smith_waterman_score(query, gap_open, gap_extend))
            .or_else(|| self.get_u32().smith_waterman_score(query, gap_open, gap_extend))
            .or_else(|| self.get_u64().smith_waterman_score(query, gap_open, gap_extend))
    }

    /// Lazily execute [`StripedDNAProfile::smith_waterman_score`] starting with
    /// the `u16` profile, skipping the `u8` profile. Lazily initializes the
    /// profiles and works its way up to the `u64` profile. Execution stops when
    /// the score returned no longer overflows the profile's integer range.
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u16<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        &self, query: &T, gap_open: u8, gap_extend: u8,
    ) -> Option<u64> {
        let query = query.as_ref();
        self.get_u16()
            .smith_waterman_score(query, gap_open, gap_extend)
            .or_else(|| self.get_u32().smith_waterman_score(query, gap_open, gap_extend))
            .or_else(|| self.get_u64().smith_waterman_score(query, gap_open, gap_extend))
    }

    /// Lazily execute [`StripedDNAProfile::smith_waterman_score`] starting with
    /// the `u32` profile, skipping the `u8` and `u16` profiles. Lazily
    /// initializes the profiles and works its way up to the `u64` profile.
    /// Execution stops when the score returned no longer overflows the
    /// profile's integer range.
    #[must_use]
    #[inline]
    pub fn smith_waterman_score_from_u32<T: AsRef<[u8]> + MaybeNucleic + ?Sized>(
        &self, query: &T, gap_open: u8, gap_extend: u8,
    ) -> Option<u64> {
        let query = query.as_ref();
        self.get_u32()
            .smith_waterman_score(query, gap_open, gap_extend)
            .or_else(|| self.get_u64().smith_waterman_score(query, gap_open, gap_extend))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::data::constants::matrices::SimpleWeightMatrix;

    #[test]
    fn sw_simd_profile_set() {
        let v: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/CY137594.txt"));
        let matrix = SimpleWeightMatrix::<5>::new_dna_matrix(2, -5, Some(b'N')).into_biased_matrix();
        let profiles = LocalDNAProfile::<32>::new(&v, &matrix);
        let profile1 = profiles.get_u8();
        let profile2 = StripedDNAProfile::<u8, 32>::new(v, &matrix);
        assert_eq!(profile1, &profile2);
    }
}
