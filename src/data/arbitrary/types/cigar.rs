//! Arbitrary implementations and specification structs for [`Cigar`],
//! [`AlignmentStates`], and related structs.

use crate::{
    alignment::{AlignmentStates, CheckedCigar, ConsumingCheckedCigar, StatesSequence},
    data::{
        arbitrary::ArbitrarySpecs,
        cigar::{Cigar, Ciglet},
    },
};
use arbitrary::{Arbitrary, Result, Unstructured};
use std::num::NonZeroUsize;

impl<'a> Arbitrary<'a> for Ciglet {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Ciglet {
            inc: usize::arbitrary(u)?,
            op:  u8::arbitrary(u)?,
        })
    }
}

impl<'a> Arbitrary<'a> for Cigar {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Cigar(Vec::arbitrary(u)?))
    }
}

impl<'a> Arbitrary<'a> for AlignmentStates {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(AlignmentStates(Vec::arbitrary(u)?))
    }
}

/// Specifications for generating an arbitrary [`Ciglet`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct CigletSpecs {
    /// Whether to restrict the increment to be non-zero.
    pub nonzero_inc: bool,

    /// Whether to restrict the operation to `MIDNSHP=X`.
    pub valid_op: bool,
}

impl<'a> ArbitrarySpecs<'a> for CigletSpecs {
    type Output = Ciglet;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let inc = if self.nonzero_inc {
            NonZeroUsize::arbitrary(u)?.get()
        } else {
            usize::arbitrary(u)?
        };

        let op = if self.valid_op {
            *u.choose(b"MIDNSHP=X")?
        } else {
            u8::arbitrary(u)?
        };

        Ok(Ciglet { inc, op })
    }
}

/// Specifications for generating an arbitrary [`AlignmentStates`] struct.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct AlignmentStatesSpecs {
    /// Whether to avoid repeated operations in adjacent ciglets.
    pub avoid_repeated_ops: bool,

    /// Ensures the number of residues consumed in the query is at most the
    /// given value.
    pub max_query_inc: Option<usize>,

    /// Ensures the number of residues consumed in the reference is at most the
    /// given value.
    pub max_ref_inc: Option<usize>,

    /// Ensures the sum of all increments is at most the given value.
    pub max_total_inc: Option<usize>,

    /// The specifications for generating each [`Ciglet`].
    pub ciglet_specs: CigletSpecs,
}

impl<'a> ArbitrarySpecs<'a> for AlignmentStatesSpecs {
    type Output = AlignmentStates;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let mut vec = self.ciglet_specs.make_arbitrary_iter(u).collect::<Result<Vec<_>>>()?;

        // Adjust the operations so that adjacent ones are not the same if
        // needed, ensuring that valid operations remain valid
        if let Some(mut last_ciglet) = vec.first().copied()
            && self.avoid_repeated_ops
        {
            for next_ciglet in &mut vec[1..] {
                if next_ciglet.op == last_ciglet.op {
                    next_ciglet.op = match last_ciglet.op {
                        b'M' => b'I',
                        b'I' => b'D',
                        b'D' => b'N',
                        b'N' => b'S',
                        b'S' => b'H',
                        b'H' => b'P',
                        b'P' => b'=',
                        b'=' => b'X',
                        b'X' => b'M',
                        other => other.wrapping_add(1),
                    }
                }
                last_ciglet = *next_ciglet;
            }
        }

        let mut states = AlignmentStates(vec);

        if let Some(max_total) = self.max_total_inc {
            states.clamp_total(max_total);
        }

        if let Some(max_query_inc) = self.max_query_inc {
            states.clamp_query_len(max_query_inc);
        }

        if let Some(max_ref_inc) = self.max_ref_inc {
            states.clamp_match_len(max_ref_inc);
        }

        Ok(states)
    }
}

/// Specifications for generating an arbitrary [`Cigar`] string.
///
/// This ensures that the CIGAR string is derived from an [`AlignmentStates`],
/// rather than just being arbitrary bytes.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct CigarSpecs {
    /// Whether to randomly insert leading zeros in front of some ciglet
    /// increments.
    pub insert_leading_zeros: bool,

    /// The specifications for generating the underlying alignment states.
    pub alignment_states_specs: AlignmentStatesSpecs,
}

impl CigarSpecs {
    /// Generates an arbitrary [`Cigar`] string, using the provided `states`
    /// instead of generating them arbitrarily from `alignment_states_specs`.
    #[inline]
    #[allow(clippy::missing_errors_doc)]
    pub fn from_states(&self, u: &mut Unstructured<'_>, states: AlignmentStates) -> Result<Cigar> {
        let cigar = if self.insert_leading_zeros {
            let leading_zeros = u.arbitrary_iter::<u8>()?.flatten().chain(std::iter::repeat(0));
            let mut out = Vec::new();
            for (ciglet, leading_zeros) in states.into_iter().zip(leading_zeros) {
                out.extend(std::iter::repeat_n(b'0', leading_zeros as usize));
                out.extend_from_slice(ciglet.inc.to_string().as_bytes());
                out.push(ciglet.op);
            }
            Cigar::from_vec_unchecked(out)
        } else {
            states.to_cigar_unchecked()
        };

        Ok(cigar)
    }
}

impl<'a> ArbitrarySpecs<'a> for CigarSpecs {
    type Output = Cigar;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let states = self.alignment_states_specs.make_arbitrary(u)?;
        self.from_states(u, states)
    }
}

/// A trait providing the ability to truncate alignments so that sums of
/// increments stay below a provided cap.
pub(crate) trait ClampAlignment {
    /// Removes ciglets to ensure summing all the increments does not exceed
    /// `max_total`.
    fn clamp_total(&mut self, max_total: usize);

    /// Removes ciglets to ensure calculating the number of residues it consumes
    /// in the query does not exceed `max_query_len`.
    fn clamp_query_len(&mut self, max_query_len: usize);

    /// Removes ciglets to ensure calculating the number of residues it consumes
    /// in the query does not exceed `max_query_len`, excluding soft clipping at
    /// valid positions at the end
    ///
    /// The soft clipping is not counted if it is the last ciglet, or it is the
    /// second to last ciglet and the last ciglet is hard clipping.
    fn clamp_query_len_exclude_tail(&mut self, max_query_len: usize);

    /// Removes ciglets to ensure calculating the number of residues it consumes
    /// in the reference does not exceed `max_match_len`.
    fn clamp_match_len(&mut self, max_match_len: usize);
}

impl ClampAlignment for AlignmentStates {
    #[inline]
    fn clamp_total(&mut self, max_total: usize) {
        let needs_shrink = match self.total_increments_checked() {
            Some(total) => total > max_total,
            None => true,
        };

        if needs_shrink {
            let mut new_vec = Vec::new();
            let mut total = 0usize;
            for ciglet in self.iter().copied() {
                total = match total.checked_add(ciglet.inc) {
                    Some(query_len) => query_len,
                    None => break,
                };
                if total > max_total {
                    break;
                }

                new_vec.push(ciglet);
            }
            *self = AlignmentStates(new_vec);
        }
    }

    #[inline]
    fn clamp_query_len(&mut self, max_query_len: usize) {
        let needs_shrink = match self.num_query_consumed_checked() {
            Some(query_len) => query_len > max_query_len,
            None => true,
        };

        if needs_shrink {
            let mut new_vec = Vec::new();
            let mut query_len = 0usize;
            for ciglet in self.iter().copied() {
                if matches!(ciglet.op, b'M' | b'I' | b'S' | b'=' | b'X') {
                    query_len = match query_len.checked_add(ciglet.inc) {
                        Some(query_len) => query_len,
                        None => break,
                    };
                    if query_len > max_query_len {
                        break;
                    }
                }
                new_vec.push(ciglet);
            }
            *self = AlignmentStates(new_vec);
        }
    }

    #[inline]
    fn clamp_query_len_exclude_tail(&mut self, max_query_len: usize) {
        let mut ciglets = self.as_slice();

        // Remove clipping without risk of overflow when adding the amounts
        // clipped
        ciglets.next_back_if_op(|op| op == b'H');
        ciglets.next_back_if_op(|op| op == b'S');

        let needs_shrink = match ciglets.num_query_consumed_checked() {
            Some(query_len) => query_len > max_query_len,
            None => true,
        };

        if needs_shrink {
            let mut new_vec = Vec::new();
            // The current number of residues consumed in the query, excuding
            // clipping at the start and end
            let mut query_len = 0usize;
            // The amount of uncounted soft clipping in new_vec that appears at
            // the end
            let mut clipping_len = 0;

            let mut ciglets = self.as_slice();

            while let Some(ciglet) = ciglets.next_ciglet() {
                if ciglet.op == b'S' {
                    // Add any previous soft clipping to the count, since it is
                    // no longer at the end
                    query_len = match query_len.checked_add(clipping_len) {
                        Some(query_len) => query_len,
                        None => break,
                    };
                    if query_len > max_query_len {
                        break;
                    }
                    // Update clipping_len with the new amount of clipping
                    clipping_len = ciglet.inc;
                    // Push soft and hard clipping ciglets
                    new_vec.push(ciglet);
                    new_vec.extend(ciglets.next_if_op(|op| op == b'H'));
                } else {
                    // The number of residues in the query consumed by this ciglet
                    let ciglet_query_len = if matches!(ciglet.op, b'M' | b'I' | b'=' | b'X') {
                        ciglet.inc
                    } else {
                        0
                    };

                    // Add any previous soft clipping to the count, since it is
                    // no longer at the end, along with ciglet_query_len
                    query_len = match query_len
                        .checked_add(clipping_len)
                        .and_then(|l| l.checked_add(ciglet_query_len))
                    {
                        Some(query_len) => query_len,
                        None => break,
                    };
                    if query_len > max_query_len {
                        break;
                    }

                    // Clear clipping_len, since we no longer end in clipping
                    clipping_len = 0;
                    new_vec.push(ciglet);
                }
            }

            *self = AlignmentStates(new_vec);
        }
    }

    #[inline]
    fn clamp_match_len(&mut self, max_match_len: usize) {
        let needs_shrink = match self.num_ref_consumed_checked() {
            Some(match_len) => match_len > max_match_len,
            None => true,
        };

        if needs_shrink {
            let mut new_vec = Vec::new();
            let mut match_len = 0usize;
            for ciglet in self.iter().copied() {
                if matches!(ciglet.op, b'M' | b'D' | b'N' | b'=' | b'X') {
                    match_len = match match_len.checked_add(ciglet.inc) {
                        Some(match_len) => match_len,
                        None => break,
                    };
                    if match_len > max_match_len {
                        break;
                    }
                }
                new_vec.push(ciglet);
            }
            *self = AlignmentStates(new_vec);
        }
    }
}

impl ClampAlignment for Cigar {
    #[inline]
    fn clamp_total(&mut self, max_total: usize) {
        let needs_shrink = match self.total_increments_checked() {
            Some(total) => total > max_total,
            None => true,
        };

        if needs_shrink {
            let mut new_vec = Vec::new();
            let mut total = 0usize;
            for ciglet in &*self {
                total = match total.checked_add(ciglet.inc) {
                    Some(query_len) => query_len,
                    None => break,
                };
                if total > max_total {
                    break;
                }
                new_vec.extend_from_slice(format!("{ciglet}").as_bytes());
            }
            *self = Cigar::from_vec_unchecked(new_vec);
        }
    }

    #[inline]
    fn clamp_query_len(&mut self, max_query_len: usize) {
        let needs_shrink = match self.num_query_consumed_checked() {
            Some(query_len) => query_len > max_query_len,
            None => true,
        };

        if needs_shrink {
            let mut new_vec = Vec::new();
            let mut query_len = 0usize;
            for ciglet in &*self {
                if matches!(ciglet.op, b'M' | b'I' | b'S' | b'=' | b'X') {
                    query_len = match query_len.checked_add(ciglet.inc) {
                        Some(query_len) => query_len,
                        None => break,
                    };
                    if query_len > max_query_len {
                        break;
                    }
                }
                new_vec.extend_from_slice(ciglet.to_string().as_bytes());
            }
            *self = Cigar::from_vec_unchecked(new_vec);
        }
    }

    #[inline]
    fn clamp_query_len_exclude_tail(&mut self, max_query_len: usize) {
        let needs_shrink = {
            // Remove clipping
            let mut ciglets = self.iter();
            ciglets.remove_clipping_back();

            // Get count without clipping at front or back
            match ciglets.num_query_consumed_checked() {
                Some(query_len) => query_len > max_query_len,
                None => true,
            }
        };

        if needs_shrink {
            let mut new_vec: Vec<u8> = Vec::new();
            // The current number of residues consumed in the query, excuding
            // clipping at the start and end
            let mut query_len = 0usize;
            // The amount of uncounted soft clipping in new_vec that appears at
            // the end
            let mut clipping_len = 0;

            let mut ciglets = self.iter();

            while let Some(ciglet) = ciglets.next_ciglet() {
                if ciglet.op == b'S' {
                    // Add any previous soft clipping to the count, since it is
                    // no longer at the end
                    query_len = match query_len.checked_add(clipping_len) {
                        Some(query_len) => query_len,
                        None => break,
                    };
                    if query_len > max_query_len {
                        break;
                    }
                    // Update clipping_len with the new amount of clipping
                    clipping_len = ciglet.inc;
                    // Push soft and hard clipping ciglets
                    new_vec.extend_from_slice(ciglet.to_string().as_bytes());
                    if let Some(ciglet) = ciglets.next_if_op(|op| op == b'H') {
                        new_vec.extend_from_slice(ciglet.to_string().as_bytes());
                    }
                } else {
                    // The number of residues in the query consumed by this ciglet
                    let ciglet_query_len = if matches!(ciglet.op, b'M' | b'I' | b'=' | b'X') {
                        ciglet.inc
                    } else {
                        0
                    };

                    // Add any previous soft clipping to the count, since it is
                    // no longer at the end, along with ciglet_query_len
                    query_len = match query_len
                        .checked_add(clipping_len)
                        .and_then(|l| l.checked_add(ciglet_query_len))
                    {
                        Some(query_len) => query_len,
                        None => break,
                    };
                    if query_len > max_query_len {
                        break;
                    }

                    // Clear clipping_len, since we no longer end in clipping
                    clipping_len = 0;
                    new_vec.extend_from_slice(ciglet.to_string().as_bytes());
                }
            }

            *self = Cigar(new_vec);
        }
    }

    #[inline]
    fn clamp_match_len(&mut self, max_match_len: usize) {
        let needs_shrink = match self.num_ref_consumed_checked() {
            Some(match_len) => match_len > max_match_len,
            None => true,
        };

        if needs_shrink {
            let mut new_vec = Vec::new();
            let mut match_len = 0usize;
            for ciglet in &*self {
                if matches!(ciglet.op, b'M' | b'D' | b'N' | b'=' | b'X') {
                    match_len = match match_len.checked_add(ciglet.inc) {
                        Some(match_len) => match_len,
                        None => break,
                    };
                    if match_len > max_match_len {
                        break;
                    }
                }
                new_vec.extend_from_slice(format!("{ciglet}").as_bytes());
            }
            *self = Cigar::from_vec_unchecked(new_vec);
        }
    }
}
