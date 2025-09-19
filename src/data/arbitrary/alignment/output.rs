use crate::{
    alignment::{Alignment, AlignmentStates, CheckedCigar, StatesSequence},
    data::arbitrary::{clamp_match_len, clamp_query_len},
    prelude::Len,
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
///   the query range plus the amount of soft clipping at the end of the
///   [`AlignmentStates`]
/// - `CLIP`: If true, ensures the full query length is accounted for with soft
///   clipping. If false, some or none of it may be marked with soft clipping
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

        if R_RANGE {
            clamp_match_len(&mut states, usize::MAX);
        }

        if Q_RANGE {
            clamp_query_len(&mut states, usize::MAX);
        }

        // Get the amount of clipping at the start and end
        let mut ciglets = states.as_slice();
        let clip_at_start = ciglets.next_if_op(|op| op == b'S').map_or(0, |c| c.inc);
        let clip_at_end = ciglets.next_back_if_op(|op| op == b'S').map_or(0, |c| c.inc);

        let ref_range = if R_RANGE {
            let len = states.num_ref_consumed_checked().unwrap();
            let start = u.int_in_range(0..=(usize::MAX - len))?;
            let end = start + len;
            start..end
        } else {
            Range::<usize>::arbitrary(u)?
        };

        let query_range = if Q_RANGE {
            // This length includes soft clipping
            let len = states.num_query_consumed_checked().unwrap();
            // The query_range should not include soft clipped portions, but
            // must be compatible so that `query_range.start - clip_at_start`
            // does not underflow and `query_range.end + clip_at_end` does not
            // overflow
            let len_with_clip_at_end = len - clip_at_start;
            let start = u.int_in_range(clip_at_start..=usize::MAX - len_with_clip_at_end)?;
            let len_no_clipping = len_with_clip_at_end - clip_at_end;
            let end = start + len_no_clipping;
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
            u.int_in_range(query_range.end + clip_at_end..=usize::MAX)?
        } else {
            usize::arbitrary(u)?
        };

        if CLIP {
            states.prepend_soft_clip(query_range.start - clip_at_start);
            states.soft_clip(query_len - (query_range.end + clip_at_end));
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

/// A combined arbitrary implementation for [`Alignment`] along with a
/// corresponding query and reference sequence. This ensures that the
/// `query_len`, `query_range`, `ref_len`, and `ref_range` fields are compatible
/// with the generated sequences and with the [`AlignmentStates`].
///
/// ## Parameters:
///
/// - `Q`: The type for the query sequence
/// - `R`: The type for the reference sequence
/// - `T`: The score type for the [`Alignment`]
/// - `S`: [`AlignmentStates`] or arbitrary wrapper for generating the alignment
///   states
/// - `CLIP`: If true, ensures the full query length is accounted for with soft
///   clipping. If false, some or none of it may be marked with soft clipping
///
/// [`AlignmentStatesArbitrary`]:
///     crate::data::arbitrary::types::AlignmentStatesArbitrary
#[derive(Debug)]
pub struct AlignmentAndSeqs<Q, R, T, S, const CLIP: bool> {
    pub alignment: Alignment<T>,
    pub query:     Q,
    pub reference: R,
    pub phantom:   PhantomData<S>,
}

impl<'a, Q, R, T, S, const CLIP: bool> Arbitrary<'a> for AlignmentAndSeqs<Q, R, T, S, CLIP>
where
    T: Arbitrary<'a>,
    Q: Arbitrary<'a> + Len,
    R: Arbitrary<'a> + Len,
    S: Arbitrary<'a> + Into<AlignmentStates>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        let query = Q::arbitrary(u)?;
        let reference = R::arbitrary(u)?;

        let mut states = S::arbitrary(u)?.into();

        clamp_match_len(&mut states, reference.len());
        clamp_query_len(&mut states, query.len());

        // Get the amount of clipping at the start and end
        let mut ciglets = states.as_slice();
        let clip_at_start = ciglets.next_if_op(|op| op == b'S').map_or(0, |c| c.inc);
        let clip_at_end = ciglets.next_back_if_op(|op| op == b'S').map_or(0, |c| c.inc);

        let ref_range = {
            let len = states.num_ref_consumed_checked().unwrap();
            let start = u.int_in_range(0..=(reference.len() - len))?;
            let end = start + len;
            start..end
        };

        let query_range = {
            let len = states.num_query_consumed_checked().unwrap();
            let len_with_clip_at_end = len - clip_at_start;
            let start = u.int_in_range(clip_at_start..=query.len() - len_with_clip_at_end)?;
            let len_no_clipping = len_with_clip_at_end - clip_at_end;
            let end = start + len_no_clipping;
            start..end
        };

        if CLIP {
            states.prepend_soft_clip(query_range.start - clip_at_start);
            states.soft_clip(query.len() - (query_range.end + clip_at_end));
        }

        Ok(AlignmentAndSeqs {
            alignment: Alignment {
                score: T::arbitrary(u)?,
                ref_range,
                query_range,
                states,
                ref_len: reference.len(),
                query_len: query.len(),
            },
            query,
            reference,
            phantom: PhantomData,
        })
    }
}
