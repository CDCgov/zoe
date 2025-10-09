use crate::{
    alignment::{Alignment, AlignmentStates, CheckedCigar},
    data::arbitrary::{ensure_no_match_len_overflow, ensure_no_query_len_overflow},
};
use arbitrary::{Arbitrary, Unstructured};
use std::{marker::PhantomData, ops::Range};

impl<'a, T: Arbitrary<'a>> Arbitrary<'a> for Alignment<T> {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        Ok(Alignment {
            score:       T::arbitrary(u)?,
            ref_range:   Range::<usize>::arbitrary(u)?,
            query_range: Range::<usize>::arbitrary(u)?,
            states:      AlignmentStates::arbitrary(u)?,
            ref_len:     usize::arbitrary(u)?,
            query_len:   usize::arbitrary(u)?,
        })
    }
}

/// A wrapper type for [`Alignment`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
/// ## Parameters:
///
/// - `T`: The score type for the [`Alignment`]
/// - `S`: [`AlignmentStates`] or arbitrary wrapper for generating the alignment
///   states
/// - `R_RANGE`: If true, ensures the length of the reference range is
///   compatible with the generated [`AlignmentStates`]
/// - `Q_RANGE`: If true, ensures the length of the query range is compatible
///   with the generated [`AlignmentStates`]
/// - `R_LEN`: If true, ensures the length of the reference is at least the end
///   of the reference range
/// - `Q_LEN`: If true, ensures the length of the query is at least the end of
///   the query range
/// - `CLIP`: If true, prepends `query_range.start` soft clipping and appends
///   `query_len-query_range.end` soft clipping
///
/// [`AlignmentStatesArbitrary`]:
///     crate::data::arbitrary::types::AlignmentStatesArbitrary
#[derive(Debug)]
pub struct AlignmentArbitrary<
    T,
    S,
    const R_RANGE: bool,
    const Q_RANGE: bool,
    const R_LEN: bool,
    const Q_LEN: bool,
    const CLIP: bool,
>(pub Alignment<T>, pub PhantomData<S>);

impl<'a, T, S, const R_RANGE: bool, const Q_RANGE: bool, const R_LEN: bool, const Q_LEN: bool, const CLIP: bool>
    Arbitrary<'a> for AlignmentArbitrary<T, S, R_RANGE, Q_RANGE, R_LEN, Q_LEN, CLIP>
where
    T: Arbitrary<'a>,
    S: Arbitrary<'a> + Into<AlignmentStates>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        let mut states = S::arbitrary(u)?.into();

        if R_RANGE && states.num_ref_consumed_checked().is_none() {
            ensure_no_match_len_overflow(&mut states);
        }

        if Q_RANGE && states.num_query_consumed_checked().is_none() {
            ensure_no_query_len_overflow(&mut states);
        }

        let ref_range = if R_RANGE {
            let len = states.num_ref_consumed_checked().unwrap();
            let start = u.int_in_range(0..=(usize::MAX - len))?;
            let end = start + len;
            start..end
        } else {
            Range::<usize>::arbitrary(u)?
        };

        let query_range = if Q_RANGE {
            let len = states.num_query_consumed_checked().unwrap();
            let start = u.int_in_range(0..=(usize::MAX - len))?;
            let end = start + len;
            start..end
        } else {
            Range::<usize>::arbitrary(u)?
        };

        let ref_len = if R_LEN {
            u.int_in_range(ref_range.end..=usize::MAX)?
        } else {
            usize::arbitrary(u)?
        };

        let query_len = if Q_LEN {
            u.int_in_range(query_range.end..=usize::MAX)?
        } else {
            usize::arbitrary(u)?
        };

        if CLIP {
            states.prepend_soft_clip(query_range.start);
            states.soft_clip(query_len - query_range.end);
        }

        Ok(AlignmentArbitrary(
            Alignment {
                score: T::arbitrary(u)?,
                ref_range,
                query_range,
                states,
                ref_len,
                query_len,
            },
            PhantomData,
        ))
    }
}
