use crate::{
    data::{
        constants::alphas::DNA_IUPAC_WITH_GAPS,
        types::nucleotides::{Nucleotides, NucleotidesReadable},
    },
    math::Uint,
};
use std::ops::{Add, AddAssign};
use std::simd::{SimdElement, prelude::*};

/// Nucleotide count statistics for A, C, G, T/U, N, - (or gaps), other valid
/// IUPAC codes, and invalid codes.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct NucleotideCounts<T: Uint> {
    inner: [T; 8],
}

impl<T> NucleotideCounts<T>
where
    T: Uint,
{
    /// Creates a new [`NucleotideCounts`] object with counts initialized to 0.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        NucleotideCounts { inner: [T::ZERO; 8] }
    }

    /// Retrives the counts as an array of size 8.
    #[inline]
    #[must_use]
    pub fn into_inner(self) -> [T; 8] {
        self.inner
    }

    /// Gets the count of the base A.
    #[inline]
    #[must_use]
    pub fn a(&self) -> T {
        self.inner[0]
    }

    /// Gets the count of the base C.
    #[inline]
    #[must_use]
    pub fn c(&self) -> T {
        self.inner[1]
    }

    /// Gets the count of the base G.
    #[inline]
    #[must_use]
    pub fn g(&self) -> T {
        self.inner[2]
    }

    /// Gets the count of the bases T and U.
    #[inline]
    #[must_use]
    pub fn t(&self) -> T {
        self.inner[3]
    }

    /// Gets the count of the base N.
    #[inline]
    #[must_use]
    pub fn n(&self) -> T {
        self.inner[4]
    }

    /// Gets the number of gaps.
    #[inline]
    #[must_use]
    pub fn gap(&self) -> T {
        self.inner[5]
    }

    /// Gets the number of other characters which are valid IUPAC.
    #[inline]
    #[must_use]
    pub fn other(&self) -> T {
        self.inner[6]
    }

    /// Gets the count of characters which are not valid IUPAC.
    #[inline]
    #[must_use]
    pub fn invalid(&self) -> T {
        self.inner[7]
    }

    /// Gets the total count of the bases A and T/U.
    #[inline]
    #[must_use]
    pub fn total_at(&self) -> T
    where
        T: Add<Output = T>, {
        self.a() + self.t()
    }

    /// Gets the total count of the bases G and C.
    #[inline]
    #[must_use]
    pub fn total_gc(&self) -> T
    where
        T: Add<Output = T>, {
        self.g() + self.c()
    }

    /// Gets the total count of the bases A, C, G, and T/U.
    #[inline]
    #[must_use]
    pub fn total_acgt(&self) -> T
    where
        T: Add<Output = T>, {
        self.inner[0..4].iter().fold(T::ZERO, |acc, &v| acc + v)
    }

    /// Gets the total count of the bases A, C, G, T/U, and N
    #[inline]
    #[must_use]
    pub fn total_acgtn(&self) -> T
    where
        T: Add<Output = T>, {
        self.inner[0..5].iter().fold(T::ZERO, |acc, &v| acc + v)
    }

    /// Gets the total count of the bases A, C, G, T/U, and N, as well as gaps.
    #[inline]
    #[must_use]
    pub fn total_acgtn_gap(&self) -> T
    where
        T: Add<Output = T>, {
        self.inner[0..6].iter().fold(T::ZERO, |acc, &v| acc + v)
    }

    /// Gets the total count of valid IUPAC codes, including regular gaps.
    #[inline]
    #[must_use]
    pub fn total_valid(&self) -> T
    where
        T: Add<Output = T>, {
        self.inner[0..7].iter().fold(T::ZERO, |acc, &v| acc + v)
    }

    /// The most frequent allele for A, C, G, T, and N. In case of ties, the
    /// order `A > C > G > T > N` is used.
    #[inline]
    #[must_use]
    pub fn plurality_acgtn(&self) -> u8 {
        let mut max_count = self.inner[0];
        let mut max_idx = 0;
        for i in 1..5 {
            if self.inner[i] > max_count {
                max_count = self.inner[i];
                max_idx = i;
            }
        }
        b"ACGTN"[max_idx]
    }
}

impl<T> Default for NucleotideCounts<T>
where
    T: Uint,
    T: Add<Output = T>,
    T: Iterator<Item = T>,
{
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

macro_rules! impl_reduce {
    {  $($ty:ty),* } => {
        $(
            impl NucleotideCounts<$ty> {
                /// Gets the total count of all characters processed.
                #[inline]
                #[must_use]
                pub fn total_any(&self) -> $ty {
                    Simd::from_array(self.inner).reduce_sum()
                }
            }
        )*
    }
}
impl_reduce!(u8, u16, u32, u64, usize);

impl<T> From<u8> for NucleotideCounts<T>
where
    T: Uint,
{
    #[inline]
    fn from(b: u8) -> Self {
        let mut nc = NucleotideCounts { inner: [T::ZERO; 8] };
        nc.inner[base_to_inner_index(b)] = T::ONE;
        nc
    }
}

impl<T> AddAssign for NucleotideCounts<T>
where
    T: Uint + SimdElement,
    Simd<T, 8>: Add<Output = Simd<T, 8>>,
{
    #[inline]
    fn add_assign(&mut self, other: Self) {
        let added = Simd::from_array(self.inner) + Simd::from_array(other.inner);
        *self = Self { inner: added.to_array() };
    }
}

impl<T> Add for NucleotideCounts<T>
where
    T: Uint + SimdElement,
    Simd<T, 8>: Add<Output = Simd<T, 8>>,
{
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        let added: Simd<T, 8> = Simd::from_array(self.inner) + Simd::from_array(other.inner);
        Self { inner: added.to_array() }
    }
}

impl<T> Add<u8> for NucleotideCounts<T>
where
    T: Uint + AddAssign,
{
    type Output = Self;

    #[inline]
    fn add(mut self, other: u8) -> Self {
        self.inner[base_to_inner_index(other)] += T::ONE;
        self
    }
}

impl<T> AddAssign<u8> for NucleotideCounts<T>
where
    T: Uint + AddAssign,
{
    #[inline]
    fn add_assign(&mut self, other: u8) {
        self.inner[base_to_inner_index(other)] += T::ONE;
    }
}

/// Extension trait to allow a per-site vector of [`NucleotideCounts`] to be
/// generated from an iterator of sequences, where the sequences are assumed to
/// be from the same multiple sequence alignment.
pub trait AlignmentComposition {
    /// Gets a per-site vector of [`NucleotideCounts`] from an iterator of
    /// sequences, where the sequences are assumed to be from the same multiple
    /// sequence alignment.
    fn per_site_counts<T>(&mut self) -> Option<Vec<NucleotideCounts<T>>>
    where
        Self: Sized + Iterator,
        Self::Item: IntoIterator<Item = u8> + Sized,
        T: Uint, {
        if let Some(first) = self.next() {
            let mut counts: Vec<NucleotideCounts<T>> = first.into_iter().map(NucleotideCounts::from).collect();

            for v in self {
                for (c, b) in counts.iter_mut().zip(v) {
                    *c += b;
                }
            }

            Some(counts)
        } else {
            None
        }
    }
}
impl<I: Iterator> AlignmentComposition for I {}

/// Converts a base into an index for [`NucleotideCounts`].
const fn base_to_inner_index(b: u8) -> usize {
    const TO_INNER: [u8; 256] = {
        const OTHER: u8 = 6;
        const INVALID: u8 = 7;

        // Default to invalid states
        let mut v = [INVALID; 256];

        // Valid IUPAC is OTHER
        let mut i = 0;
        while i < DNA_IUPAC_WITH_GAPS.len() {
            v[DNA_IUPAC_WITH_GAPS[i] as usize] = OTHER;
            i += 1;
        }

        // Map specifically to our array
        let from = b"acgtunACGTUN-";
        let dest = b"0123340123345";
        i = 0;
        while i < from.len() {
            v[from[i] as usize] = dest[i] - b'0';
            i += 1;
        }

        v
    };

    TO_INNER[b as usize] as usize
}

/// Extension trait for creating a [`NucleotideCounts`] object from a sequence.
pub trait ToBaseCounts: NucleotidesReadable {
    /// Creates a [`NucleotideCounts`] object from a sequence using the
    /// specified [`Uint`].
    #[inline]
    #[must_use]
    fn to_base_counts<T: Uint>(&self) -> NucleotideCounts<T> {
        self.nucleotide_bytes()
            .iter()
            .fold(NucleotideCounts::new(), |acc, &b| acc + b)
    }
}

impl<T: NucleotidesReadable> ToBaseCounts for T {}

/// Extension trait to allow for a consensus sequence to be generated from a
/// vector of [`NucleotideCounts`]. See [`plurality_acgtn`] for more details.
///
/// [`plurality_acgtn`]: CreateConsensus::plurality_acgtn
pub trait CreateConsensus {
    /// Given a vector of [`NucleotideCounts`] representing the counts in each
    /// position of a multiple sequence alignment, generate a consensus
    /// sequence. See [`NucleotideCounts::plurality_acgtn`] for more details.
    #[must_use]
    fn plurality_acgtn(&self) -> Nucleotides;
}

impl<T> CreateConsensus for Vec<NucleotideCounts<T>>
where
    T: Uint + Ord + Add<Output = T>,
{
    #[inline]
    fn plurality_acgtn(&self) -> Nucleotides {
        self.iter().map(NucleotideCounts::plurality_acgtn).collect()
    }
}

#[cfg(test)]
mod test;
