use crate::{
    alignment::AlignmentIndices,
    data::types::cigar::{Cigar, Ciglet},
};

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

    /// Returns the [`Ciglet`] elements as a slice.
    #[inline]
    #[must_use]
    pub fn as_slice(&self) -> &[Ciglet] {
        self.0.as_slice()
    }

    /// Adds a state to the right end of the alignment.
    ///
    /// If the operation is the same as the rightmost operation, this state is
    /// included in the previous [`Ciglet`].
    pub fn add_state(&mut self, op: u8) {
        self.add_ciglet(Ciglet { inc: 1, op });
    }

    /// Adds a [`Ciglet`] to the right end of the alignment.
    ///
    /// If the operation is the same as the rightmost operation, the ciglet is
    /// merged with the last one. If the increment is 0, no change occurs.
    pub fn add_ciglet(&mut self, ciglet: Ciglet) {
        if ciglet.inc > 0 {
            if let Some(c) = self.0.last_mut()
                && c.op == ciglet.op
            {
                c.inc += ciglet.inc;
            } else {
                self.0.push(ciglet);
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
    pub fn add_inc_op(&mut self, inc: usize, op: u8) {
        self.add_ciglet(Ciglet { inc, op });
    }

    /// Creates an alignment consisting only of `M` (match states) without any
    /// gaps for the reference or query. Soft clipping is added as needed. Match
    /// states may include mismatches.
    #[cfg(feature = "dev-3pass")]
    pub(crate) fn new_no_gaps(aligned_range: std::ops::Range<usize>, query_len: usize) -> Self {
        let mut states = AlignmentStates::with_capacity(3);

        states.soft_clip(aligned_range.start);
        states.add_inc_op(aligned_range.end - aligned_range.start, b'M');
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
        if inc > 0 {
            if let Some(c) = self.0.last_mut()
                && c.op == b'S'
            {
                c.inc += inc;
            } else {
                self.0.push(Ciglet { inc, op: b'S' });
            }
        }
    }

    /// Adds soft clipping `S` to the start of the alignment `inc` times.
    ///
    /// If the leftmost operation is `S`, `inc` is added to its `increment`. If
    /// `inc` is 0, no change occurs.
    pub fn prepend_soft_clip(&mut self, inc: usize) {
        if inc > 0 {
            if let Some(c) = self.0.first_mut()
                && c.op == b'S'
            {
                c.inc += inc;
            } else {
                self.0.insert(0, Ciglet { inc, op: b'S' });
            }
        }
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
        Self(cigar.iter().collect())
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
            if ciglet.op == b'M' {
                let query_slice = &query[query_index..query_index + ciglet.inc];
                let ref_slice = &reference[ref_index..ref_index + ciglet.inc];

                for (query_base, ref_base) in query_slice.iter().zip(ref_slice) {
                    out.add_state(if query_base == ref_base { b'=' } else { b'X' });
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

/// Iterator yielding the aligned bases as specified by the given alignment
/// operations.
///
/// The first base is from the reference, and the second base is from the query.
/// Gaps are represented by `None`.
pub struct AlignmentIter<'a, I>
where
    I: Iterator<Item = Ciglet>, {
    reference_buffer: &'a [u8],
    query_buffer:     &'a [u8],
    ciglets:          I,
    inc:              usize,
    op:               u8,
}

impl<'a, I> AlignmentIter<'a, I>
where
    I: Iterator<Item = Ciglet>,
{
    /// Creates a new [`AlignmentIter`] from a reference, query, cigar
    /// string, and reference position.
    #[inline]
    #[must_use]
    pub(crate) fn new(
        reference: &'a [u8], query: &'a [u8], ciglets: impl IntoIterator<Item = Ciglet, IntoIter = I>, ref_index: usize,
    ) -> Self {
        let reference_buffer = &reference[ref_index..];
        let query_buffer = query;
        let mut ciglets = ciglets.into_iter();
        // If no valid ciglets, initialize with inc at 0 so that iterator is empty
        let Ciglet { inc, op } = ciglets.next().unwrap_or(Ciglet { inc: 0, op: b'M' });

        AlignmentIter {
            reference_buffer,
            query_buffer,
            ciglets,
            inc,
            op,
        }
    }

    /// Get the next operation from the iterator. This will advance the
    /// [`Ciglet`] iterator if necessary. `None` is returned if the iterator has
    /// reached its end.
    #[inline]
    #[must_use]
    fn get_next_op(&mut self) -> Option<u8> {
        if self.inc > 0 {
            self.inc -= 1;
            Some(self.op)
        } else {
            let Ciglet { inc, op } = self.ciglets.next()?;
            self.inc = inc;
            self.op = op;
            self.get_next_op()
        }
    }

    /// Forcibly skip to the next Ciglet in the iterator.
    #[inline]
    #[must_use]
    fn skip_to_next_ciglet(&mut self) -> Option<()> {
        let Ciglet { inc, op } = self.ciglets.next()?;
        self.inc = inc;
        self.op = op;
        Some(())
    }

    /// Remove a base from the beginning of the reference buffer.
    ///
    /// ## Panics
    ///
    /// The reference must contain at least one base.
    #[inline]
    #[must_use]
    fn advance_reference(&mut self) -> u8 {
        let out = self.reference_buffer[0];
        self.reference_buffer = &self.reference_buffer[1..];
        out
    }

    /// Remove a base from the beginning of the query buffer.
    ///
    /// ## Panics
    ///
    /// The query must contain at least one base.
    #[inline]
    #[must_use]
    fn advance_query(&mut self) -> u8 {
        let out = self.query_buffer[0];
        self.query_buffer = &self.query_buffer[1..];
        out
    }
}

impl<I> Iterator for AlignmentIter<'_, I>
where
    I: Iterator<Item = Ciglet>,
{
    type Item = (Option<u8>, Option<u8>);

    /// # Panics
    ///
    /// The reference and query must be at least as long as the length specified
    /// in the iterator's alignment operations. All CIGAR operations must be in
    /// `MIDNSHP=X`.
    fn next(&mut self) -> Option<Self::Item> {
        let op = self.get_next_op()?;

        match op {
            b'M' | b'=' | b'X' => Some((Some(self.advance_reference()), Some(self.advance_query()))),
            b'D' => Some((Some(self.advance_reference()), None)),
            b'I' => Some((None, Some(self.advance_query()))),
            b'S' => {
                // Skip this Ciglet. The query buffer is advanced by inc+1 since
                // inc was decremented in get_next_op
                self.query_buffer = &self.query_buffer[self.inc + 1..];
                self.skip_to_next_ciglet()?;
                self.next()
            }
            b'N' => Some((Some(self.advance_reference()), Some(b'N'))),
            b'H' | b'P' => {
                // Skip this Ciglet without modifying either buffer
                self.skip_to_next_ciglet()?;
                self.next()
            }
            _ => panic!("CIGAR op '{op}' not supported.\n"),
        }
    }
}

/// A trait representing a sequence of alignment states, where each element is
/// represented as a [`Ciglet`].
///
/// This provides iterator-like functionality, allowing the [`Ciglet`] elements
/// to be consumed (such as with [`next_ciglet`] and [`next_ciglet_back`]) or
/// peeked at (such as with [`peek_op`] or [`peek_back_op`]).
///
/// [`next_ciglet`]: StatesSequence::next_ciglet
/// [`next_ciglet_back`]: StatesSequence::next_ciglet_back
/// [`peek_op`]: StatesSequence::peek_op
/// [`peek_back_op`]: StatesSequence::peek_back_op
pub trait StatesSequence {
    /// Peeks at the operator for the next ciglet without consuming it.
    #[must_use]
    fn peek_op(&self) -> Option<u8>;

    /// Peeks at the operator for the last ciglet without consuming it.
    #[must_use]
    fn peek_back_op(&self) -> Option<u8>;

    /// Checks whether the sequence of alignment states is empty.
    ///
    /// This assumes that the states are valid and have a non-zero increment.
    #[must_use]
    fn is_empty(&self) -> bool;

    /// Retrieves the next [`Ciglet`] and removes it from the
    /// [`StatesSequence`], similar to [`Iterator::next`].
    fn next_ciglet(&mut self) -> Option<Ciglet>;

    /// Retrieves the next [`Ciglet`] from the end and removes it from the
    /// [`StatesSequence`], similar to [`DoubleEndedIterator::next_back`].
    fn next_ciglet_back(&mut self) -> Option<Ciglet>;

    /// Gets the next ciglet if the operator meets the specified predicate,
    /// otherwise the [`StatesSequence`] is not modified.
    #[inline]
    fn next_if_op(&mut self, f: impl FnOnce(u8) -> bool) -> Option<Ciglet> {
        if f(self.peek_op()?) { self.next_ciglet() } else { None }
    }

    /// Gets the last ciglet if the operator meets the specified predicate,
    /// otherwise the [`StatesSequence`] is not modified.
    #[inline]
    fn next_back_if_op(&mut self, f: impl FnOnce(u8) -> bool) -> Option<Ciglet> {
        if f(self.peek_back_op()?) {
            self.next_ciglet_back()
        } else {
            None
        }
    }

    /// Removes clipping from the start of the iterator.
    ///
    /// First, a hard clipping ciglet is removed if present. Then a soft
    /// clipping ciglet is removed if present. The total number of bases
    /// clipping is returned.
    #[inline]
    fn remove_clipping_front(&mut self) -> usize {
        if let Some(ciglet1) = self.next_if_op(|op| op == b'H' || op == b'S') {
            if ciglet1.op == b'H'
                && let Some(ciglet2) = self.next_if_op(|op| op == b'S')
            {
                ciglet1.inc + ciglet2.inc
            } else {
                ciglet1.inc
            }
        } else {
            0
        }
    }

    /// Removes clipping from the end of the iterator.
    ///
    /// First, a hard clipping ciglet is removed if present. Then a soft
    /// clipping ciglet is removed if present. The total number of bases
    /// clipping is returned.
    #[inline]
    fn remove_clipping_back(&mut self) -> usize {
        if let Some(ciglet1) = self.next_back_if_op(|op| op == b'H' || op == b'S') {
            if ciglet1.op == b'H'
                && let Some(ciglet2) = self.next_back_if_op(|op| op == b'S')
            {
                ciglet1.inc + ciglet2.inc
            } else {
                ciglet1.inc
            }
        } else {
            0
        }
    }
}

impl StatesSequence for &[Ciglet] {
    #[inline]
    fn peek_op(&self) -> Option<u8> {
        self.first().map(|ciglet| ciglet.op)
    }

    #[inline]
    fn peek_back_op(&self) -> Option<u8> {
        self.last().map(|ciglet| ciglet.op)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        (*self).is_empty()
    }

    #[inline]
    fn next_ciglet(&mut self) -> Option<Ciglet> {
        let (first, rest) = self.split_first()?;
        *self = rest;
        Some(*first)
    }

    #[inline]
    fn next_ciglet_back(&mut self) -> Option<Ciglet> {
        let (last, rest) = self.split_last()?;
        *self = rest;
        Some(*last)
    }
}

/// A trait for computing the number of residues consumed in the reference or
/// query for a given CIGAR string or [`AlignmentStates`], while checking for
/// overflow. This is helpful for randomly generated alignments.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *fuzzing* feature in your `Cargo.toml` to use this
/// trait.
///
/// </div>
#[cfg(feature = "fuzzing")]
pub trait CheckedCigar
where
    for<'a> &'a Self: IntoIterator<Item: std::borrow::Borrow<Ciglet>>, {
    #[inline]
    #[must_use]
    fn total_increments_checked(&self) -> Option<usize> {
        self.into_iter()
            .map(|ciglet| *std::borrow::Borrow::borrow(&ciglet))
            .try_fold(0usize, |sum, ciglet| sum.checked_add(ciglet.inc))
    }

    /// Sums the lengths for operations consuming the reference (`M`, `D`, `N`,
    /// `=`, and `X`), returning `None` for overflow.
    #[inline]
    #[must_use]
    fn num_ref_consumed_checked(&self) -> Option<usize> {
        self.into_iter()
            .map(|ciglet| *std::borrow::Borrow::borrow(&ciglet))
            .filter_map(|Ciglet { inc, op }| matches!(op, b'M' | b'D' | b'N' | b'=' | b'X').then_some(inc))
            .try_fold(0, usize::checked_add)
    }

    /// Sums the lengths for operations consuming the query (`M`, `I`, `S`, `=`,
    /// and `X`), returning `None` for overflow.
    #[must_use]
    fn num_query_consumed_checked(&self) -> Option<usize> {
        self.into_iter()
            .map(|ciglet| *std::borrow::Borrow::borrow(&ciglet))
            .filter_map(|Ciglet { inc, op }| matches!(op, b'M' | b'I' | b'S' | b'=' | b'X').then_some(inc))
            .try_fold(0, usize::checked_add)
    }
}

#[cfg(feature = "fuzzing")]
impl<T> CheckedCigar for T where for<'a> &'a T: IntoIterator<Item: std::borrow::Borrow<Ciglet>> {}
