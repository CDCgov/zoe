//! A module providing implementations of
//! [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
//! for many of *Zoe*'s types, as well as other tools for generating arbitrary
//! testing data.
//!
//! *Zoe*'s implementations of [`Arbitrary`] generate the data with the least
//! guarantees, excluding:
//!
//! - Safety guarantees, such as graphic ASCII for [`QualityScores`]
//! - Some basic guarantees that would otherwise prevent any utility, such as
//!   [`ThreeBitKmerSet`] containing the encoded k-mers of the proper size
//!
//! To generate data with stronger assumptions, *Zoe* uses the following
//! framework:
//!
//! 1. In your fuzz test file, define an input struct containing the arbitrary
//!    types you would like to generate. This could be a wrapper around a single
//!    *Zoe* type, or it could include multiple types.
//! 2. Begin implementing [`Arbitrary`] for the input struct. For any types not
//!    requiring additional guarantees, generate these using
//!    [`Arbitrary::arbitrary`].
//! 3. For any types requiring more guarantees, find the corresponding
//!    specifications struct in this module. For example, to specify the
//!    guarantees for a CIGAR string, use [`CigarSpecs`].
//! 4. Within the implementation, initialize the specifications struct. Any
//!    guarantees that should be upheld will need to be specified, and then
//!    `..Default::default()` can be used to initialize the rest of the fields
//!    (the default for the specification structs is no guarantees).
//! 5. Call [`ArbitrarySpecs::make_arbitrary`] on the specifications struct to
//!    generate the arbitrary data.
//!
//! [`QualityScores`]: crate::prelude::QualityScores
//! [`ThreeBitKmerSet`]: crate::kmer::encoders::three_bit::ThreeBitKmerSet

use arbitrary::{Arbitrary, Result, Unstructured};
use std::{collections::HashMap, marker::PhantomData, ops::RangeInclusive};

mod alignment;
mod byte;
mod collections;
mod kmer;
mod math;
mod records;
mod string;
mod types;

pub use alignment::*;
pub use byte::*;
pub use collections::*;
pub use math::*;
pub use records::*;
pub use string::*;
pub use types::*;

/// Implemented on structs which provide specifications for how to generate
/// arbitrary data.
///
/// See [`arbitrary`] for more details.
///
/// [`arbitrary`]: crate::data::arbitrary
pub trait ArbitrarySpecs<'a>: Sized {
    /// The type of the arbitrary data that this specification generates.
    type Output;

    /// Generates arbitrary data conforming to the given specifications.
    ///
    /// ## Errors
    ///
    /// Any errors from the underlying [`arbitrary`] calls are propagated.
    ///
    /// [`arbitrary`]: Arbitrary::arbitrary
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output>;

    /// Generates an iterator of arbitrary data, each item conforming to the
    /// given specifications.
    #[inline]
    fn make_arbitrary_iter<'b>(&self, u: &'b mut Unstructured<'a>) -> MakeArbitraryIter<'a, 'b, '_, Self> {
        MakeArbitraryIter { specs: self, u }
    }
}

/// A zero-sized type representing the specs for generating `T` without any
/// constraints (i.e., using the standard [`Arbitrary`] impl).
///
/// This is provided for compatibility with generic interfaces, like
/// [`VecSpecs`] which expects the element specs to implement
/// [`ArbitrarySpecs`].
pub struct NoConstraintSpecs<T> {
    inner: PhantomData<T>,
}

impl<T> Default for NoConstraintSpecs<T> {
    fn default() -> Self {
        Self { inner: PhantomData }
    }
}

impl<'a, T> ArbitrarySpecs<'a> for NoConstraintSpecs<T>
where
    T: Arbitrary<'a>,
{
    type Output = T;

    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        T::arbitrary(u)
    }
}

impl<'a, S> ArbitrarySpecs<'a> for Option<S>
where
    S: ArbitrarySpecs<'a, Output: Arbitrary<'a>>,
{
    type Output = S::Output;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        if let Some(specs) = self {
            specs.make_arbitrary(u)
        } else {
            S::Output::arbitrary(u)
        }
    }
}

/// Utility iterator produced by [`ArbitrarySpecs::make_arbitrary_iter`].
///
/// Each call to [`next`] arbitrarily determines whether to continue the
/// iterator or end it. The iterator is not fused.
///
/// [`next`]: Iterator::next
#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct MakeArbitraryIter<'a, 'b, 'c, Specs> {
    specs: &'c Specs,
    u:     &'b mut Unstructured<'a>,
}

impl<'a, Specs: ArbitrarySpecs<'a>> Iterator for MakeArbitraryIter<'a, '_, '_, Specs> {
    type Item = Result<Specs::Output>;

    #[inline]
    fn next(&mut self) -> Option<Result<Specs::Output>> {
        let keep_going = self.u.arbitrary().unwrap_or(false);
        if keep_going {
            Some(self.specs.make_arbitrary(self.u))
        } else {
            None
        }
    }
}

/// An extension trait for [`Unstructured`] providing additional methods.
pub trait UnstructuredExt {
    /// Chooses `n` distinct values from an inclusive range.
    ///
    /// ## Errors
    ///
    /// Propagates errors from the underlying arbitrary integer generation.
    ///
    /// ## Panics
    ///
    /// Panics when the range contains fewer than `n` distinct values.
    fn distinct_in_range(&mut self, range: RangeInclusive<usize>, n: usize) -> Result<Vec<usize>>;

    /// Divides `sum` into `num_summands` positive integers.
    ///
    /// ## Errors
    ///
    /// Propagates errors from [`Self::distinct_in_range`].
    ///
    /// ## Panics
    ///
    /// Panics if `num_summands` is not in `1..=sum`, unless both `sum` and
    /// `num_summands` are 0.
    fn integer_composition(&mut self, sum: usize, num_summands: usize) -> Result<Vec<usize>> {
        if sum == 0 && num_summands == 0 {
            return Ok(Vec::new());
        }
        assert!(
            num_summands > 0,
            "Cannot write {sum} (which is non-zero) as the sum of 0 positive summands"
        );
        assert!(
            num_summands <= sum,
            "Cannot write {sum} as the sum of {num_summands} positive summands"
        );

        let mut cuts = self.distinct_in_range(1..=(sum - 1), num_summands - 1)?;
        cuts.sort_unstable();

        let mut last = 0;
        let mut summands = Vec::with_capacity(num_summands);

        for cut in cuts {
            summands.push(cut - last);
            last = cut;
        }
        summands.push(sum - last);

        Ok(summands)
    }

    /// Randomly interleaves two exact-size iterators while preserving the
    /// ordering within each iterator.
    ///
    /// ## Errors
    ///
    /// Propagates errors from the underlying arbitrary byte generation.
    fn interleave<T>(
        &mut self, left: impl ExactSizeIterator<Item = T>, right: impl ExactSizeIterator<Item = T>,
    ) -> Result<Vec<T>>;
}

impl UnstructuredExt for Unstructured<'_> {
    fn distinct_in_range(&mut self, range: RangeInclusive<usize>, n: usize) -> Result<Vec<usize>> {
        let Some(n_m1) = n.checked_sub(1) else {
            return Ok(Vec::new());
        };

        assert!(range.start() <= range.end(), "Cannot choose values from an empty range");
        let mut remaining_m1 = *range.end() - *range.start();

        assert!(
            remaining_m1 >= n_m1,
            "There are not {n} distinct values in the provided range"
        );

        let mut swaps = HashMap::with_capacity(n);
        let mut values = Vec::with_capacity(n);

        for draw in 0..n {
            let index = self.int_in_range(0..=remaining_m1)?;
            let picked = swaps.get(&index).copied().unwrap_or(index);

            if index != remaining_m1 {
                let last_value = swaps.get(&remaining_m1).copied().unwrap_or(remaining_m1);
                swaps.insert(index, last_value);
            }
            swaps.remove(&remaining_m1);

            values.push(*range.start() + picked);
            if draw < n_m1 {
                remaining_m1 -= 1;
            }
        }

        Ok(values)
    }

    fn interleave<T>(
        &mut self, left: impl ExactSizeIterator<Item = T>, right: impl ExactSizeIterator<Item = T>,
    ) -> Result<Vec<T>> {
        let (mut left, mut right) = (left.into_iter(), right.into_iter());

        let total = left.len() + right.len();
        let mut bits = vec![0; total.saturating_sub(1).div_ceil(8)];
        // fill_buffer pads with zeroes if the input runs out.
        self.fill_buffer(&mut bits)?;

        let mut out = Vec::with_capacity(total);
        let mut choice = 0;

        while left.len() > 0 && right.len() > 0 {
            let take_left = (bits[choice / 8] >> (choice % 8)) & 1 == 0;
            choice += 1;

            if take_left {
                out.push(left.next().unwrap());
            } else {
                out.push(right.next().unwrap());
            }
        }

        out.extend(left);
        out.extend(right);
        Ok(out)
    }
}
