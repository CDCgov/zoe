//! A module defining [`AlignmentStates`] and methods on it.

use crate::{
    alignment::{AlignOp, AlignmentIndices},
    data::cigar::{Cigar, Ciglet},
};
use std::ops::Range;

/// A struct for storing alignment states.
///
/// This is similar to [`Cigar`] (and can be converted to one). However, instead
/// of storing the bytes for the CIGAR string, [`AlignmentStates`] stores
/// increment-operation pairs as a vector of [`Ciglet`] structs. This allows for
/// less parsing/checking during use.
///
/// ## Validity
///
/// This struct does not guarantee that the operations in each [`Ciglet`] are
/// valid.
///
/// [`AlignmentStates`] ensures that the increments are non-zero and that
/// adjacent [`Ciglet`] values have distinct operations. There are some
/// unchecked functions which may invalidate this, although those functions have
/// documented validity sections. Any arbitrary implementations may ignore these
/// assumptions as well.
///
/// ## Panics
///
/// For all methods adding additional operations to the [`AlignmentStates`], a
/// panic due to overflow will occur if more than [`usize::MAX`] of the same
/// operation appear in a row.
///
/// ## Limitations
///
/// Prepending methods do a simple [`Vec::insert`] on the first index, which
/// copies elements.
///
/// ## Conversion from String-Like Types
///
/// When converting from a string-like type to [`AlignmentStates`], `*` is
/// parsed as empty. Otherwise, the following assumptions must be met (otherwise
/// an error is returned):
///
/// - The CIGAR operations must be among `M, I, D, N, S, H, P, X, =`.
/// - Every operation in the CIGAR string must have a preceding increment.
/// - Every increment must be followed by an operation.
/// - The increment for each operation must be non-zero and less than or equal
///   to [`usize::MAX`].
/// - The CIGAR increment after merging consecutive ciglets with the same
///   operation must be less than or equal to [`usize::MAX`].
#[derive(Clone, Eq, PartialEq, Default)]
pub struct AlignmentStates(pub(crate) Vec<Ciglet>);

impl AlignmentStates {
    /// Initializes an empty alignment.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        AlignmentStates(Vec::new())
    }

    /// Initializes the states with capacity for `n` increment-operation pairs.
    #[inline]
    #[must_use]
    pub fn with_capacity(n: usize) -> Self {
        AlignmentStates(Vec::with_capacity(n))
    }

    /// Number of increment-operation pairs in the alignment.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::alignment::AlignmentStates;
    /// let states = AlignmentStates::try_from(b"3S10M1D9M").unwrap();
    /// assert_eq!(states.len(), 4);
    /// ```
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether the [`AlignmentStates`] contains no data (or zero
    /// alignment states).
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the [`Ciglet`] elements as a slice.
    #[inline]
    #[must_use]
    pub fn as_slice(&self) -> &[Ciglet] {
        self.0.as_slice()
    }

    /// Returns the [`Ciglet`] elements as a mutable slice.
    ///
    /// ## Validity
    ///
    /// Any mutations performed should ensure that increments are non-zero and
    /// that adjacent [`Ciglet`] values have distinct operations. Otherwise, the
    /// assumptions of [`AlignmentStates`](AlignmentStates#validity) may be
    /// invalidated.
    #[inline]
    #[must_use]
    pub fn as_mut_slice(&mut self) -> &mut [Ciglet] {
        self.0.as_mut_slice()
    }

    /// Returns a mutable reference to the vector of [`Ciglet`] elements.
    ///
    /// ## Validity
    ///
    /// Any mutations performed should ensure that increments are non-zero and
    /// that adjacent [`Ciglet`] values have distinct operations. Otherwise, the
    /// assumptions of [`AlignmentStates`](AlignmentStates#validity) may be
    /// invalidated.
    #[inline]
    #[must_use]
    pub fn as_mut_vec(&mut self) -> &mut Vec<Ciglet> {
        &mut self.0
    }

    /// Adds a state to the rightmost end of the alignment, merging state where
    /// appropriate.
    pub fn add_state(&mut self, op: AlignOp) {
        self.add_ciglet(Ciglet { inc: 1, op });
    }

    /// Adds a state to the leftmost end of the alignment, merging state where
    /// appropriate.
    pub fn prepend_state(&mut self, op: AlignOp) {
        self.prepend_ciglet(Ciglet { inc: 1, op });
    }

    /// Adds a [`Ciglet`] to the right end of the alignment, merging [`Ciglet`]s
    /// where appropriate.
    pub fn add_ciglet(&mut self, ciglet: Ciglet) {
        if ciglet.inc > 0 {
            if let Some(last) = self.0.last_mut()
                && last.op == ciglet.op
            {
                last.inc = last.inc.strict_add(ciglet.inc);
            } else {
                self.0.push(ciglet);
            }
        }
    }

    /// Adds a [`Ciglet`] to the leftmost end of the alignment, merging
    /// [`Ciglet`]s where appropriate.
    pub fn prepend_ciglet(&mut self, ciglet: Ciglet) {
        if ciglet.inc > 0 {
            if let Some(first) = self.0.first_mut()
                && first.op == ciglet.op
            {
                first.inc = first.inc.strict_add(ciglet.inc);
            } else {
                self.0.insert(0, ciglet);
            }
        }
    }

    /// Adds an increment-operation pair.
    ///
    /// This is equivalent to calling [`add_state`] `inc` times, or calling
    /// [`add_ciglet`] after combining the information into a [`Ciglet`].
    ///
    /// If the operation is the same as the rightmost operation, the ciglet is
    /// merged with the last one. If the increment is 0, no change occurs.
    ///
    /// [`add_state`]: AlignmentStates::add_state
    /// [`add_ciglet`]: AlignmentStates::add_ciglet
    #[inline]
    pub fn add_inc_op(&mut self, inc: usize, op: AlignOp) {
        self.add_ciglet(Ciglet { inc, op });
    }

    /// Prepends an increment-operation pair.
    ///
    /// This is equivalent to calling [`prepend_state`] `inc` times, or calling
    /// [`prepend_ciglet`] after combining the information into a [`Ciglet`].
    ///
    /// If the operation is the same as the leftmost operation, the ciglet is
    /// merged with the first one. If the increment is 0, no change occurs.
    ///
    /// [`prepend_state`]: AlignmentStates::prepend_state
    /// [`prepend_ciglet`]: AlignmentStates::prepend_ciglet
    #[inline]
    pub fn prepend_inc_op(&mut self, inc: usize, op: AlignOp) {
        self.prepend_ciglet(Ciglet { inc, op });
    }

    /// Creates an alignment consisting only of `M` (match states) without any
    /// gaps for the reference or query. Soft clipping is added as needed. Match
    /// states may include mismatches.
    pub(crate) fn new_no_gaps(aligned_range: Range<usize>, query_len: usize) -> Self {
        let mut states = AlignmentStates::with_capacity(3);

        states.soft_clip(aligned_range.start);
        states.add_inc_op(aligned_range.end - aligned_range.start, AlignOp::M);
        states.soft_clip(query_len - aligned_range.end);
        states
    }

    /// Extends the [`AlignmentStates`] with an iterator of [`Ciglet`] values.
    ///
    /// If the rightmost operation is the same as the first operation in the
    /// iterator, they are merged.
    ///
    /// ## Validity
    ///
    /// Adjacent ciglets in the iterator must have distinct operations and
    /// non-zero increments.
    pub(crate) fn extend_from_ciglets<I>(&mut self, ciglets: I)
    where
        I: IntoIterator<Item = Ciglet>, {
        let mut ciglets = ciglets.into_iter();
        let Some(first_ciglet) = ciglets.next() else { return };
        self.add_ciglet(first_ciglet);
        self.0.extend(ciglets);
    }

    /// Adds soft clipping `S` to the end of the alignment `inc` times.
    ///
    /// If the rightmost operation is `S`, `inc` is added to its `increment`. If
    /// `inc` is 0, no change occurs.
    pub fn soft_clip(&mut self, inc: usize) {
        self.add_ciglet(Ciglet { inc, op: AlignOp::S });
    }

    /// Adds soft clipping `S` to the start of the alignment `inc` times.
    ///
    /// If the leftmost operation is `S`, `inc` is added to its `increment`. If
    /// `inc` is 0, no change occurs.
    pub fn prepend_soft_clip(&mut self, inc: usize) {
        self.prepend_ciglet(Ciglet { inc, op: AlignOp::S });
    }

    /// Converts the [`AlignmentStates`] struct to a [`Cigar`] string, without
    /// checking for valid operations.
    #[must_use]
    pub fn to_cigar_unchecked(&self) -> Cigar {
        Cigar::from_ciglets_unchecked(self.0.iter().copied())
    }

    /// Collects an iterator of [`Ciglet`] values into an [`AlignmentStates`]
    /// struct without checking.
    ///
    /// ## Validity
    ///
    /// - `ciglets` must not contain adjacent operations that are equal
    /// - `ciglets` must not contain any increments that are zero
    /// - If `ciglets` contains any invalid operations, increment overflows, or
    ///   missing operations, the output may be truncated
    #[must_use]
    pub fn from_ciglets_unchecked<I: IntoIterator<Item = Ciglet>>(ciglets: I) -> Self {
        Self(ciglets.into_iter().collect())
    }

    /// Converts the [`Cigar`] string into an [`AlignmentStates`] struct without
    /// checking.
    ///
    /// ## Validity
    ///
    /// - `cigar` must not contain adjacent operations that are equal
    /// - `cigar` must not contain any increments that are zero
    /// - If `cigar` contains any invalid operations, increment overflows, or
    ///   missing operations, the output may be truncated
    #[must_use]
    pub fn from_cigar_unchecked(cigar: &Cigar) -> Self {
        Self(cigar.iter_unchecked().collect())
    }

    /// Reverses the order of the stored alignment states in-place.
    #[inline]
    pub fn make_reverse(&mut self) {
        self.0.reverse();
    }

    /// Returns an [`AlignmentStates`] with the order of the states reversed.
    #[inline]
    #[must_use]
    pub fn to_reverse(&self) -> Self {
        // Validity: reversing the ciglets will not alter the increments or
        // operations other than their order. Since the input does not have any
        // equal adjacent operations, the reversed ciglets also will not have
        // any
        Self::from_ciglets_unchecked(self.into_iter().rev())
    }

    /// Yields an iterator over the alignment states.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, Ciglet> {
        self.0.iter()
    }

    /// Generates a new alignment states replacing `M` with either `=` for
    /// matches and `X` for mismatches.
    ///
    /// Consider using [`Alignment::to_verbose_sequence_matching`] if you have
    /// an [`Alignment`] object.
    ///
    /// The `reference` should be the entire reference, and `ref_index` is the
    /// starting position of the alignment within the passed reference.
    ///
    /// ## Panics
    ///
    /// If either `reference` or `query` are of a shorter length than implied by
    /// the alignment, then this will panic due to out of bounds indexing.
    ///
    /// [`Alignment`]: super::Alignment
    /// [`Alignment::to_verbose_sequence_matching`]:
    ///     super::Alignment::to_verbose_sequence_matching
    #[must_use]
    pub fn to_verbose_sequence_matching(&self, reference: &[u8], query: &[u8], mut ref_index: usize) -> Self {
        let mut out = AlignmentStates::with_capacity(self.0.len());
        let mut query_index = 0;

        for ciglet in self {
            if ciglet.op == AlignOp::M {
                let query_slice = &query[query_index..query_index + ciglet.inc];
                let ref_slice = &reference[ref_index..ref_index + ciglet.inc];

                for (query_base, ref_base) in query_slice.iter().zip(ref_slice) {
                    out.add_state(if query_base == ref_base { AlignOp::Eq } else { AlignOp::X });
                }

                query_index += ciglet.inc;
                ref_index += ciglet.inc;
            } else {
                out.add_ciglet(ciglet);
                (query_index, ref_index) = (query_index, ref_index).increment_idxs_by(ciglet);
            }
        }

        out
    }
}
