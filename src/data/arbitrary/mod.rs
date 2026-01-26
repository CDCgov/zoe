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
//! <div class="warning note">
//!
//! **Note**
//!
//! You must enable the *fuzzing* feature in your `Cargo.toml` to use this API.
//!
//! </div>
//!
//! [`QualityScores`]: crate::prelude::QualityScores
//! [`ThreeBitKmerSet`]: crate::kmer::encoders::three_bit::ThreeBitKmerSet

use arbitrary::{Arbitrary, Result, Unstructured};

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
