//! Arbitrary implementations and specification structs for the output of
//! alignment algorithms.
//!
//! To also generate the corresponding sequences, see [`AlignmentAndSeqsSpecs`]
//! and [`AlignmentAndQuerySpecs`].

use crate::{
    alignment::{Alignment, AlignmentStates, CheckedCigar, StatesSequence, StatesSequenceMut},
    data::arbitrary::{AlignmentStatesSpecs, ArbitrarySpecs, ByteSpecs, ClampAlignment, VecSpecs},
    prelude::Len,
};
use arbitrary::{Arbitrary, Result, Unstructured};
use std::{marker::PhantomData, ops::Range};

impl<'a, T: Arbitrary<'a>> Arbitrary<'a> for Alignment<T> {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
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

/// Specifications for generating an arbitrary [`Alignment`].
///
/// ## Parameters
///
/// `Score` is the numeric type to use for the `score` field.
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct AlignmentSpecs<Score> {
    /// The specifications for generating the [`AlignmentStates`].
    pub states_specs: AlignmentStatesSpecs,

    /// Ensures that the `ref_range` field is compatible with the number of
    /// residues consumed from the reference in the `states` field.
    pub compatible_ref_range_and_states: bool,

    /// Ensures that the `ref_len` field is at least as large as the number of
    /// residues consumed by `states` from the reference.
    pub compatible_ref_len_and_states: bool,

    /// Ensures that the `ref_len` field is at least the end of `ref_range`.
    pub compatible_ref_len_and_range: bool,

    /// Ensures that the `query_range` field has a length equal to the number of
    /// residues consumed by `states` from the query, excluding clipping at the
    /// start or end. Furthermore, the start of the range must be at least the
    /// amount of soft clipping at the start, and adding the soft clipping at
    /// the end to the end of the query range must not overflow a `usize`.
    pub compatible_query_range_and_states: bool,

    /// Ensures that the `query_len` field is at least as large as the number of
    /// residues consumed by `states` from the query, including soft clipping at
    /// the start and end.
    pub compatible_query_len_and_states: bool,

    /// If set, ensures that the `query_len` field is at least the end of
    /// `query_range`.
    pub compatible_query_len_and_range: bool,

    /// Whether to annotate the beginning and end of the alignment with soft
    /// clipping so that the full query length is accounted for.
    ///
    /// To ensure that this works properly, you may also want to set
    /// `compatible_query_len_and_range`.
    pub annotate_clipping: bool,

    /// The maximum value that the start of `ref_range` can be.
    ///
    /// For example, if a 1-based position is calculated from the start of
    /// `ref_range`, then capping it at `usize::MAX - 1` would prevent overflow.
    pub max_ref_range_start: usize,

    /// The maximum value that the end of `ref_range` can be.
    pub max_ref_range_end: usize,

    /// The maximum value that `ref_len` can be.
    pub max_ref_len: usize,

    /// The exact value that `ref_len` should be.
    ///
    /// If this is set, then `max_ref_len` is ignored.
    pub ref_len: Option<usize>,

    /// The maximum value that the start of `query_range` can be.
    pub max_query_range_start: usize,

    /// The maximum value that the end of `query_range` can be.
    pub max_query_range_end: usize,

    /// The maximum value that `query_len` can be.
    pub max_query_len: usize,

    /// The exact value that `query_len` should be.
    ///
    /// If this is set, the `max_query_len` is ignored.
    pub query_len: Option<usize>,

    /// The type for the `score` field.
    pub score: PhantomData<Score>,
}

impl<Score> Default for AlignmentSpecs<Score> {
    /// Generates the default specifications for [`AlignmentSpecs`], which
    /// enforces no guarantees.
    #[inline]
    fn default() -> Self {
        Self {
            states_specs: AlignmentStatesSpecs::default(),
            compatible_ref_range_and_states: false,
            compatible_ref_len_and_states: false,
            compatible_ref_len_and_range: false,
            compatible_query_range_and_states: false,
            compatible_query_len_and_states: false,
            compatible_query_len_and_range: false,
            annotate_clipping: false,
            max_ref_range_start: usize::MAX,
            max_ref_range_end: usize::MAX,
            max_ref_len: usize::MAX,
            ref_len: None,
            max_query_range_start: usize::MAX,
            max_query_range_end: usize::MAX,
            max_query_len: usize::MAX,
            query_len: None,
            score: PhantomData,
        }
    }
}

impl<'a, Score> ArbitrarySpecs<'a> for AlignmentSpecs<Score>
where
    Score: Arbitrary<'a>,
{
    type Output = Alignment<Score>;

    #[inline]
    #[allow(clippy::too_many_lines)]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        // Generate the arbitrary alignment
        let mut states = self.states_specs.make_arbitrary(u)?;

        // Copy the maximum bounds, which we will then update as needed
        let mut max_ref_len = self.max_ref_len;
        let mut max_ref_range_end = self.max_ref_range_end;
        let mut max_ref_range_start = self.max_ref_range_start;
        let mut max_query_len = self.max_query_len;
        let mut max_query_range_end = self.max_query_range_end;
        let mut max_query_range_start = self.max_query_range_start;

        // If exact lengths are set for ref_len or query_len, set the maximum
        // bounds equal to those
        if let Some(ref_len) = self.ref_len {
            max_ref_len = ref_len;
        }
        if let Some(query_len) = self.query_len {
            max_query_len = query_len;
        }

        // Make bounds more restrictive based on compatibilities
        if self.compatible_ref_len_and_range {
            max_ref_range_end = max_ref_range_end.min(max_ref_len);
        }
        max_ref_range_start = max_ref_range_start.min(max_ref_range_end);

        if self.compatible_query_len_and_range {
            max_query_range_end = max_query_range_end.min(max_query_len);
        }
        max_query_range_start = max_query_range_start.min(max_query_range_end);

        // Clamp the states if needed
        if self.compatible_ref_range_and_states {
            states.clamp_match_len(max_ref_range_end);
        }
        if self.compatible_ref_len_and_states {
            states.clamp_match_len(max_ref_len);
        }
        if self.compatible_query_range_and_states {
            // Restrict so that the soft clipping at the back isn't counted, but
            // the clipping at the front is
            states.clamp_query_len_exclude_tail(max_query_range_end);
            // When clipping at the front and end are included, we are not
            // allowed to overflow usize::MAX
            states.clamp_query_len(usize::MAX);
        }
        if self.compatible_query_len_and_states {
            states.clamp_query_len(max_query_len);
        }

        // Generate the reference range: generate the starting value, then
        // compute/generate the end
        let ref_range = if self.compatible_ref_range_and_states {
            // Validity: the number of bases consumed in the reference was
            // clamped above when compatible_ref_range_and_states is set
            let len = states.num_ref_consumed_checked().unwrap();
            // Validity: subtraction will not overflow because we clamped len to
            // at most max_ref_range_end
            max_ref_range_start = max_ref_range_start.min(max_ref_range_end - len);

            let start = u.int_in_range(0..=max_ref_range_start)?;
            let end = start + len;
            start..end
        } else {
            let start = u.int_in_range(0..=max_ref_range_start)?;
            let end = u.int_in_range(start..=max_ref_range_end)?;
            start..end
        };

        // Get the soft clipping ciglets at the start and end. The options are
        // mutable to allow the inner references to be mutated without requiring
        // destructuring of the option
        let mut ciglets = states.as_mut_slice();
        ciglets.next_if_op(|op| op == b'H');
        let mut soft_clip_at_start_ciglet = ciglets.next_if_op_mut(|op| op == b'S');
        ciglets.next_back_if_op(|op| op == b'H');
        let mut soft_clip_at_end_ciglet = ciglets.next_back_if_op_mut(|op| op == b'S');

        // Generate the query range: generate the starting value, then
        // compute/generate the end
        let query_range = if self.compatible_query_range_and_states {
            // Get the number of residues consumed, excluding clipping which may
            // be mutated when generating query_range. Validity: the number of
            // residues consumed in the query was clamped above when
            // compatible_query_range_and_states is set.
            let num_query_consumed_no_clipping = ciglets.num_query_consumed_checked().unwrap();

            // Decrease max_query_range_start as needed to ensure a random start
            // index can have an end index computed that does not exceed
            // max_query_range_end
            max_query_range_start = max_query_range_start.min(max_query_range_end - num_query_consumed_no_clipping);

            // Decrease max_query_range_start as needed to ensure a random start
            // index can have an end index computed such that the end index plus
            // the amount of soft clipping at the end does not overflow a usize
            let soft_clip_at_end = soft_clip_at_end_ciglet.as_ref().map_or(0, |c| c.inc);
            let len_with_clip_at_end = num_query_consumed_no_clipping + soft_clip_at_end;
            max_query_range_start = max_query_range_start.min(usize::MAX - len_with_clip_at_end);

            // Decrease the amount of clipping at the start to ensure a start
            // index exists such that subtracting the soft clipping does not
            // underflow
            if let Some(ref mut soft_clip_at_start_ciglet) = soft_clip_at_start_ciglet
                && soft_clip_at_start_ciglet.inc > max_query_range_start
            {
                soft_clip_at_start_ciglet.inc = max_query_range_start;
            }

            // Get the amount of soft clipping at the start
            let soft_clip_at_start = soft_clip_at_start_ciglet.as_ref().map_or(0, |c| c.inc);

            // Generate the random start index
            let start = u.int_in_range(soft_clip_at_start..=max_query_range_start)?;

            // Compute an end index, which will not exceed max_query_range_end
            let end = start + num_query_consumed_no_clipping;
            start..end
        } else {
            let start = u.int_in_range(0..=max_query_range_start)?;
            let end = u.int_in_range(start..=max_query_range_end)?;
            start..end
        };

        // Get the number of query residues consumed, including clipping
        let num_query_consumed = soft_clip_at_start_ciglet.as_ref().map_or(0, |c| c.inc)
            + ciglets.num_query_consumed_checked().unwrap()
            + soft_clip_at_end_ciglet.as_ref().map_or(0, |c| c.inc);

        // Generate the reference length if it isn't fixed
        let ref_len = if let Some(ref_len) = self.ref_len {
            ref_len
        } else {
            let mut min_ref_len = 0;
            if self.compatible_ref_len_and_states {
                // Validity: compatible_ref_len_and_states causes the states to
                // have a clamped match_len
                min_ref_len = min_ref_len.max(ciglets.num_ref_consumed_checked().unwrap());
            }
            if self.compatible_ref_len_and_range {
                min_ref_len = min_ref_len.max(ref_range.end);
            }

            u.int_in_range(min_ref_len..=max_ref_len)?
        };

        // Generate the query length if it isn't fixed
        let query_len = if let Some(query_len) = self.query_len {
            query_len
        } else {
            let mut min_query_len = 0;
            if self.compatible_query_len_and_states {
                min_query_len = min_query_len.max(num_query_consumed);
            }
            if self.compatible_query_len_and_range {
                min_query_len = min_query_len.max(query_range.end);
            }

            u.int_in_range(min_query_len..=max_query_len)?
        };

        // Annotate the soft clipping
        if self.annotate_clipping {
            // Separate out the prepending logic until after ciglet references
            // have been dropped to appease borrow checker
            let prepend_soft_clip = if let Some(ciglet) = soft_clip_at_start_ciglet {
                ciglet.inc = query_range.start;
                false
            } else {
                true
            };

            // If query_len and query_range aren't agreeing, there isn't much we
            // can do, so do nothing. This is noted in the docs for
            // annotate_clipping.
            if let Some(missing_clipping_at_end) = query_len.checked_sub(query_range.end) {
                if let Some(ref mut ciglet) = soft_clip_at_end_ciglet {
                    ciglet.inc = missing_clipping_at_end;
                } else {
                    states.soft_clip(missing_clipping_at_end);
                }
            }

            if prepend_soft_clip {
                states.prepend_soft_clip(query_range.start);
            }
        }

        Ok(Alignment {
            score: Score::arbitrary(u)?,
            ref_range,
            query_range,
            states,
            ref_len,
            query_len,
        })
    }
}

/// Specifications for generating an arbitrary alignment with compatible
/// sequences, as an [`AlignmentAndSeqs`].
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct AlignmentAndSeqsSpecs<Score, SeqSpecs> {
    /// The specifications for generating the query sequence.
    pub query_specs: SeqSpecs,

    /// The specifications for generating the reference sequence.
    pub reference_specs: SeqSpecs,

    /// The specifications for generating the [`Alignment`].
    ///
    /// Before generation, this will be mutated so that the true `query_len` and
    /// `ref_len` values are set in the specifications. However, the
    /// compatibility flags must still be set, as well as any other constraints.
    pub alignment_specs: AlignmentSpecs<Score>,
}

impl<'a, Score, SeqSpecs> ArbitrarySpecs<'a> for AlignmentAndSeqsSpecs<Score, SeqSpecs>
where
    Score: Arbitrary<'a> + Clone,
    SeqSpecs: ArbitrarySpecs<'a, Output: Len>,
{
    type Output = AlignmentAndSeqs<Score, SeqSpecs::Output>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        let reference = self.reference_specs.make_arbitrary(u)?;
        let query = self.query_specs.make_arbitrary(u)?;

        let mut specs = self.alignment_specs.clone();
        specs.ref_len = Some(reference.len());
        specs.query_len = Some(query.len());

        let alignment = specs.make_arbitrary(u)?;

        Ok(AlignmentAndSeqs {
            alignment,
            query,
            reference,
        })
    }
}

/// Specifications for generating an arbitrary alignment and a query sequence,
/// as an [`AlignmentAndQuery`] struct.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct AlignmentAndQuerySpecs<Score> {
    /// The specifications for generating the query sequence.
    pub query_specs: VecSpecs<ByteSpecs>,

    /// The specifications for generating the [`Alignment`].
    ///
    /// Before generation, this will be mutated so that the true `query_len`
    /// value is set in the specifications. However, the compatibility flags
    /// must still be set, as well as any other constraints.
    pub alignment_specs: AlignmentSpecs<Score>,
}

impl<'a, Score> ArbitrarySpecs<'a> for AlignmentAndQuerySpecs<Score>
where
    Score: Arbitrary<'a> + Clone,
{
    type Output = AlignmentAndQuery<Score>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        let query = self.query_specs.make_arbitrary(u)?;

        let mut specs = self.alignment_specs.clone();
        specs.query_len = Some(query.len());

        let alignment = specs.make_arbitrary(u)?;

        Ok(AlignmentAndQuery { alignment, query })
    }
}

/// A struct containing an [`Alignment`] along with the corresponding sequences.
#[derive(Clone, Eq, PartialEq, Debug, Default)]
pub struct AlignmentAndSeqs<Score, Seq> {
    /// The alignment struct.
    pub alignment: Alignment<Score>,
    /// The query sequence, which is the same length as `alignment.query_len`.
    pub query:     Seq,
    /// The reference sequence, which is the same length as `alignment.ref_len`.
    pub reference: Seq,
}

/// A struct containing an [`Alignment`] along with the corresponding query.
#[derive(Clone, Eq, PartialEq, Debug, Default)]
pub struct AlignmentAndQuery<Score> {
    /// The alignment struct.
    pub alignment: Alignment<Score>,
    /// The query sequence, which is the same length as `alignment.query_len`.
    pub query:     Vec<u8>,
}
