//! Arbitrary implementations and specification structs for [`Cigar`],
//! [`AlignmentStates`], and related structs.

use crate::{
    alignment::{AlignmentStates, NextCiglet},
    data::{
        arbitrary::{ArbitrarySpecs, ByteSet},
        cigar::{Cigar, Ciglet, FormatCigletForCigarVec, LenInAlignment},
    },
};
use arbitrary::{Arbitrary, Result, Unstructured};
use std::cmp::Ordering;

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
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct CigletSpecs {
    /// Whether to restrict the increment to be non-zero.
    pub nonzero_inc: bool,

    /// The maximum increment to allow.
    pub max_inc: usize,

    /// Whether to restrict the operation to `MIDNSHP=X`.
    pub valid_op: bool,

    /// Whether to allow `M` as an operation.
    pub include_m: bool,

    /// Whether to allow `I` as an operation.
    pub include_i: bool,

    /// Whether to allow `D` as an operation.
    pub include_d: bool,

    /// Whether to allow `N` as an operation.
    pub include_n: bool,

    /// Whether to allow `S` as an operation.
    pub include_s: bool,

    /// Whether to allow `H` as an operation.
    pub include_h: bool,

    /// Whether to allow `P` as an operation.
    pub include_p: bool,

    /// Whether to allow `=` as an operation.
    pub include_eq: bool,

    /// Whether to allow `X` as an operation.
    pub include_x: bool,
}

impl Default for CigletSpecs {
    fn default() -> Self {
        Self {
            nonzero_inc: false,
            max_inc:     usize::MAX,
            valid_op:    false,
            include_m:   true,
            include_i:   true,
            include_d:   true,
            include_n:   true,
            include_s:   true,
            include_h:   true,
            include_p:   true,
            include_eq:  true,
            include_x:   true,
        }
    }
}

impl CigletSpecs {
    /// Builds the [`ByteSet`] corresponding to the operation.
    fn op_byte_set(&self) -> ByteSet {
        let include_op = [
            (b'M', self.include_m),
            (b'I', self.include_i),
            (b'D', self.include_d),
            (b'N', self.include_n),
            (b'S', self.include_s),
            (b'H', self.include_h),
            (b'P', self.include_p),
            (b'=', self.include_eq),
            (b'X', self.include_x),
        ];

        if self.valid_op {
            // This is the common case, so make it more efficient by
            // initializing a vector up front with specified capacity and
            // linearly building up ByteSet

            let mut options = Vec::with_capacity(include_op.len());

            for (op, include) in include_op {
                if include {
                    options.push(op);
                }
            }

            ByteSet::from(options)
        } else {
            // This is an uncommon case, if any include fields are false. In
            // that case, we run in quadratic time when building the ByteSet

            let mut byte_set = ByteSet::Any;

            for (byte, include) in include_op {
                if !include {
                    byte_set.remove(byte);
                }
            }

            byte_set
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for CigletSpecs {
    type Output = Ciglet;

    /// Generates an arbitrary [`Ciglet`] conforming to the given
    /// specifications.
    ///
    /// ## Errors
    ///
    /// - Any errors from the underlying [`arbitrary`] calls are propagated.
    /// - If `valid_op` is `true` and all `include_*` fields are `false`,
    ///   [`EmptyChoose`] is returned.
    ///
    /// ## Panics
    ///
    /// This panics if `nonzero_inc` is true and `max_inc` is 0.
    ///
    /// [`arbitrary`]: arbitrary::Arbitrary::arbitrary
    /// [`EmptyChoose`]: arbitrary::Error::EmptyChoose
    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let min_inc = usize::from(self.nonzero_inc);

        let inc = match min_inc.cmp(&self.max_inc) {
            Ordering::Less | Ordering::Equal => u.int_in_range(min_inc..=self.max_inc)?,
            Ordering::Greater => panic!("max_inc was 0 and nonzero_inc was true"),
        };

        let byte_set = self.op_byte_set();
        let op = byte_set.make_arbitrary(u)?;

        Ok(Ciglet { inc, op })
    }
}

/// Specifications for generating an arbitrary [`AlignmentStates`] struct.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct AlignmentStatesSpecs {
    /// Whether to avoid repeated operations in adjacent ciglets.
    pub avoid_repeated_ops: bool,

    /// Whether to avoid placing hard clipping in the middle of the alignment
    /// (except when solely flanked on a side by other hard clipping ciglets).
    pub avoid_middle_hard_clipping: bool,

    /// Whether to avoid placing soft clipping in the middle of the alignment
    /// (except when solely flanked on a side by other soft clipping ciglets,
    /// potentially followed by hard clipping ciglets).
    pub avoid_middle_soft_clipping: bool,

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

    /// Generates an arbitrary [`AlignmentStates`] conforming to the given
    /// specifications.
    ///
    /// ## Errors
    ///
    /// Any errors from the underlying [`arbitrary`] calls are propagated.
    ///
    /// ## Panics
    ///
    /// This panics if [`make_arbitrary`] for the `ciglet_specs` panics.
    ///
    /// [`arbitrary`]: arbitrary::Arbitrary::arbitrary
    /// [`make_arbitrary`]: ArbitrarySpecs::make_arbitrary
    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let mut vec = self.ciglet_specs.make_arbitrary_iter(u).collect::<Result<Vec<_>>>()?;

        if self.avoid_middle_hard_clipping {
            let mut iter = vec.as_slice().iter();

            let hard_clip_start = iter.by_ref().take_while(|ciglet| ciglet.op == b'H').count();
            let hard_clip_back = iter.by_ref().rev().take_while(|ciglet| ciglet.op == b'H').count();

            let start = hard_clip_start;
            let end = vec.len() - hard_clip_back;

            let mut i = 0;
            vec.retain(|ciglet| {
                let keep = !(start..end).contains(&i) || ciglet.op != b'H';
                i += 1;
                keep
            });
        }

        if self.avoid_middle_soft_clipping {
            let mut iter = vec.as_slice().iter();

            let clip_start = {
                let mut clipping_op = b'H';
                iter.by_ref()
                    .take_while(|ciglet| {
                        if ciglet.op == clipping_op {
                            true
                        } else if ciglet.op == b'S' && clipping_op == b'H' {
                            clipping_op = b'S';
                            true
                        } else {
                            false
                        }
                    })
                    .count()
            };

            let clip_back = {
                let mut clipping_op = b'H';
                iter.by_ref()
                    .rev()
                    .take_while(|ciglet| {
                        if ciglet.op == clipping_op {
                            true
                        } else if ciglet.op == b'S' && clipping_op == b'H' {
                            clipping_op = b'S';
                            true
                        } else {
                            false
                        }
                    })
                    .count()
            };

            let start = clip_start;
            let end = vec.len() - clip_back;

            let mut i = 0;
            vec.retain(|ciglet| {
                let keep = !(start..end).contains(&i) || ciglet.op != b'S';
                i += 1;
                keep
            });
        }

        if self.avoid_repeated_ops {
            let mut last_ciglet: Option<Ciglet> = None;

            vec.retain(|ciglet| {
                if last_ciglet.map(|ciglet| ciglet.op) == Some(ciglet.op) {
                    false
                } else {
                    last_ciglet = Some(*ciglet);
                    true
                }
            });
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
                out.push_formatted_ciglet(ciglet);
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
        let total_increments = self.iter().try_fold(0usize, |sum, ciglet| sum.checked_add(ciglet.inc));

        let needs_shrink = match total_increments {
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
        let needs_shrink = match self.query_len_in_alignment_checked() {
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
        ciglets.next_ciglet_back_if_op(|op| op == b'H');
        ciglets.next_ciglet_back_if_op(|op| op == b'S');

        let needs_shrink = match ciglets.query_len_in_alignment_checked() {
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
                    new_vec.extend(ciglets.next_ciglet_if_op(|op| op == b'H'));
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
        let needs_shrink = match self.ref_len_in_alignment_checked() {
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
        let total_increments = self.iter().try_fold(0usize, |sum, ciglet| sum.checked_add(ciglet.inc));

        let needs_shrink = match total_increments {
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
                new_vec.push_formatted_ciglet(ciglet);
            }
            *self = Cigar::from_vec_unchecked(new_vec);
        }
    }

    #[inline]
    fn clamp_query_len(&mut self, max_query_len: usize) {
        let needs_shrink = match self.query_len_in_alignment_checked() {
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
                new_vec.push_formatted_ciglet(ciglet);
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
            match ciglets.query_len_in_alignment_checked() {
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
                    new_vec.push_formatted_ciglet(ciglet);
                    if let Some(ciglet) = ciglets.next_ciglet_if_op(|op| op == b'H') {
                        new_vec.push_formatted_ciglet(ciglet);
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
                    new_vec.push_formatted_ciglet(ciglet);
                }
            }

            *self = Cigar(new_vec);
        }
    }

    #[inline]
    fn clamp_match_len(&mut self, max_match_len: usize) {
        let needs_shrink = match self.ref_len_in_alignment_checked() {
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
                new_vec.push_formatted_ciglet(ciglet);
            }
            *self = Cigar::from_vec_unchecked(new_vec);
        }
    }
}
