use crate::data::types::nucleotides::Nucleotides;
use crate::data::{constants::alphas::NUCLEIC_IUPAC, types::Uint};
use std::ops::{Add, AddAssign};
use std::simd::prelude::*;

/// Nucleotide count statistics for A, C, G, T/U, N, -, other, invalid.
#[derive(Debug, Clone, Copy)]
pub struct NucleotideCounts<T: Uint> {
    inner: [T; 8],
}

impl<T> NucleotideCounts<T>
where
    T: Uint,
{
    #[must_use]
    pub fn new() -> Self {
        NucleotideCounts { inner: [T::zero(); 8] }
    }

    #[must_use]
    pub fn into_inner(self) -> [T; 8] {
        self.inner
    }

    #[inline]
    pub fn a(&self) -> T {
        self.inner[0]
    }
    #[inline]
    pub fn c(&self) -> T {
        self.inner[1]
    }
    #[inline]
    pub fn g(&self) -> T {
        self.inner[2]
    }
    #[inline]
    pub fn t(&self) -> T {
        self.inner[3]
    }
    #[inline]
    pub fn n(&self) -> T {
        self.inner[4]
    }
    #[inline]
    pub fn gap(&self) -> T {
        self.inner[5]
    }
    #[inline]
    pub fn other(&self) -> T {
        self.inner[6]
    }
    #[inline]
    pub fn invalid(&self) -> T {
        self.inner[7]
    }

    /// Total count for A, T
    #[inline]
    pub fn total_at(&self) -> T
    where
        T: Add<Output = T>, {
        self.a() + self.t()
    }

    /// Total count for A, T
    #[inline]
    pub fn total_gc(&self) -> T
    where
        T: Add<Output = T>, {
        self.g() + self.c()
    }

    /// Total count for A, C, G, T
    #[inline]
    pub fn total_actg(&self) -> T
    where
        T: Add<Output = T>, {
        self.inner[0..4].iter().fold(T::zero(), |acc, &v| acc + v)
    }

    /// Total count for A, C, G, T, N
    #[inline]
    pub fn total_actgn(&self) -> T
    where
        T: Add<Output = T>, {
        self.inner[0..5].iter().fold(T::zero(), |acc, &v| acc + v)
    }

    /// Total count for A, C, G, T, N
    #[inline]
    pub fn total_acgtn_gap(&self) -> T
    where
        T: Add<Output = T>, {
        self.inner[0..5].iter().fold(T::zero(), |acc, &v| acc + v)
    }

    /// Total valid IUPAC codes, including regular gaps.
    #[inline]
    pub fn total_valid(&self) -> T
    where
        T: Add<Output = T>, {
        self.inner[0..7].iter().fold(T::zero(), |acc, &v| acc + v)
    }

    /// The most frequent allele for A
    #[inline]
    pub fn plurality_acgtn(&self) -> u8 {
        const TO_BASE: &[u8; 6] = b"ACGTN-";
        const N: usize = length_needed(b'N');

        let mut counts = [BaseCount {
            base:  0,
            count: T::zero(),
            index: 0,
        }; N];
        for i in 0..N {
            counts[i] = BaseCount {
                base:  TO_BASE[i],
                count: self.inner[i],
                index: i,
            };
        }

        counts.sort_by(|a, b| b.cmp(a));
        counts[0].base
    }
}

impl<T> Default for NucleotideCounts<T>
where
    T: Uint,
    T: Add<Output = T>,
    T: Iterator<Item = T>,
{
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
struct BaseCount<T> {
    base:  u8,
    count: T,
    index: usize,
}

impl<T> PartialOrd for BaseCount<T>
where
    T: PartialOrd,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.count.partial_cmp(&other.count)
    }
}

impl<T> Ord for BaseCount<T>
where
    T: Ord + Eq,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.count.cmp(&other.count).then_with(|| self.index.cmp(&other.index))
    }
}

macro_rules! impl_reduce {
    {  $($ty:ty),* } => {
        $(
            impl NucleotideCounts<$ty> {
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
        let mut nc = NucleotideCounts { inner: [T::zero(); 8] };
        nc.inner[base_to_inner_index(b)] = T::one();
        nc
    }
}

impl<T> AddAssign for NucleotideCounts<T>
where
    T: Uint,
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
    T: Uint,
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
        self.inner[base_to_inner_index(other)] += T::one();
        self
    }
}

impl<T> AddAssign<u8> for NucleotideCounts<T>
where
    T: Uint + AddAssign,
{
    #[inline]
    fn add_assign(&mut self, other: u8) {
        self.inner[base_to_inner_index(other)] += T::one();
    }
}

pub trait AlignmentComposition {
    fn per_site_counts<T>(&mut self) -> Option<Vec<NucleotideCounts<T>>>
    where
        Self: Sized + Iterator,
        Self::Item: IntoIterator<Item = u8> + Sized,
        T: Uint,
        NucleotideCounts<T>: AddAssign<u8>, {
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

const fn base_to_inner_index(b: u8) -> usize {
    const TO_INNER: [u8; 256] = {
        const OTHER: u8 = 6;
        const INVALID: u8 = 7;

        // Default to invalid states
        let mut v = [INVALID; 256];

        // Valid IUPAC is OTHER
        let mut i = 0;
        while i < NUCLEIC_IUPAC.len() {
            v[NUCLEIC_IUPAC[i] as usize] = OTHER;
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

const fn length_needed(b: u8) -> usize {
    base_to_inner_index(b) + 1
}

pub trait CreateConsensus {
    fn plurality_acgtn(&self) -> Nucleotides;
}

impl<T> CreateConsensus for Vec<NucleotideCounts<T>>
where
    T: Uint + Ord + Add<Output = T>,
{
    fn plurality_acgtn(&self) -> Nucleotides {
        self.iter().map(NucleotideCounts::plurality_acgtn).collect()
    }
}
