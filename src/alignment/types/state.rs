use crate::{
    alignment::AlignmentIndices,
    data::types::cigar::{Cigar, Ciglet},
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
/// ## Limitations
///
/// Prepending methods do a simple [`Vec::insert`] on the first index, which
/// copies elements.
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
    pub fn add_state(&mut self, op: u8) {
        self.add_ciglet(Ciglet { inc: 1, op });
    }

    /// Adds a state to the leftmost end of the alignment, merging state where
    /// appropriate.
    pub fn prepend_state(&mut self, op: u8) {
        self.prepend_ciglet(Ciglet { inc: 1, op });
    }

    /// Adds a [`Ciglet`] to the right end of the alignment, merging [`Ciglet`]s
    /// where appropriate.
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

    /// Adds a [`Ciglet`] to the leftmost end of the alignment, merging
    /// [`Ciglet`]s where appropriate.
    pub fn prepend_ciglet(&mut self, ciglet: Ciglet) {
        if ciglet.inc > 0 {
            if let Some(c) = self.0.first_mut()
                && c.op == ciglet.op
            {
                c.inc += ciglet.inc;
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
    pub fn add_inc_op(&mut self, inc: usize, op: u8) {
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
    pub fn prepend_inc_op(&mut self, inc: usize, op: u8) {
        self.prepend_ciglet(Ciglet { inc, op });
    }

    /// Creates an alignment consisting only of `M` (match states) without any
    /// gaps for the reference or query. Soft clipping is added as needed. Match
    /// states may include mismatches.
    pub(crate) fn new_no_gaps(aligned_range: Range<usize>, query_len: usize) -> Self {
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
        self.prepend_ciglet(Ciglet { inc, op: b'S' });
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

/// An extension trait for alignment-like data that enables the next operation
/// to be peeked at (without consuming it, in the case of an iterator).
///
/// Mutable access to `self` is required solely for the purpose of fetching the
/// relevent data or advancing past empty ciglets.
pub trait PeekOp {
    /// Peeks at the operation for the next ciglet without consuming it. Empty
    /// ciglets are skipped and may be removed.
    #[must_use]
    fn peek_op(&mut self) -> Option<u8>;

    /// Peeks at the operation for the last ciglet without consuming it. Empty
    /// ciglets are skipped and may be removed.
    #[must_use]
    fn peek_op_back(&mut self) -> Option<u8>;
}

/// An extension trait for alignment-like data that enables a single increment
/// of the the next operation to be removed and returned.
///
/// This will mutate the underlying slice or iterator to contain one fewer
/// increment, and may result in an empty [`Ciglet`].
pub trait TakeOp: PeekOp {
    /// Removes and returns one increment of the next operation from the
    /// alignment states.
    #[must_use]
    fn take_op(&mut self) -> Option<u8>;

    /// Removes and returns one increment of the last operation from the
    /// alignment states.
    #[must_use]
    fn take_op_back(&mut self) -> Option<u8>;

    /// Removes and returns the next operation if it meets the specified
    /// predicate, otherwise `self` is not modified (aside from possibly
    /// removing empty ciglets).
    #[inline]
    fn take_op_if(&mut self, f: impl FnOnce(u8) -> bool) -> Option<u8> {
        if f(self.peek_op()?) { self.take_op() } else { None }
    }

    /// Removes and returns the last operation if it meets the specified
    /// predicate, otherwise `self` is not modified (aside from possibly
    /// removing empty ciglets).
    #[inline]
    fn take_op_back_if(&mut self, f: impl FnOnce(u8) -> bool) -> Option<u8> {
        if f(self.peek_op_back()?) { self.take_op_back() } else { None }
    }
}

/// An extension trait for alignment-like data that enables the next [`Ciglet`]
/// to be peeked at (without consuming it, in the case of an iterator).
///
/// Mutable access to `self` is required solely for the purpose of fetching the
/// relevent data or advancing past empty ciglets.
pub trait PeekCiglet {
    /// Peeks at the next ciglet without consuming it. Empty ciglets are
    /// skipped and may be removed.
    fn peek_ciglet(&mut self) -> Option<Ciglet>;

    /// Peeks at the last ciglet without consuming it. Empty ciglets are
    /// skipped and may be removed.
    fn peek_ciglet_back(&mut self) -> Option<Ciglet>;
}

/// Similar to [`PeekCiglet`] but providing mutable access to the peeked
/// [`Ciglet`].
pub trait PeekCigletMut {
    /// Peeks at the next ciglet without consuming it, as well as giving mutable
    /// access to it. Empty ciglets are skipped and may be removed.
    fn peek_ciglet_mut(&mut self) -> Option<&mut Ciglet>;

    /// Peeks at the last ciglet without consuming it, as well as giving mutable
    /// access to it. Empty ciglets are skipped and may be removed.
    fn peek_ciglet_back_mut(&mut self) -> Option<&mut Ciglet>;
}

/// An extension trait for alignment-like data that enables the next ciglet to
/// be consumed and returned.
///
/// This provides similar functionality to an iterator. When implemented on
/// slice-like types, consuming the ciglet involves shrinking the slice without
/// mutating the data it is refering to.
pub trait NextCiglet: PeekOp {
    /// Retrieves the next [`Ciglet`] and removes it from `self`, similar to
    /// [`Iterator::next`]. Empty ciglets are skipped.
    fn next_ciglet(&mut self) -> Option<Ciglet>;

    /// Retrieves the last [`Ciglet`] and removes it from `self`, similar to
    /// [`DoubleEndedIterator::next_back`]. Empty ciglets are skipped.
    fn next_ciglet_back(&mut self) -> Option<Ciglet>;

    /// Gets the next ciglet if the operation meets the specified predicate,
    /// otherwise `self` is not modified (aside from skipping empty ciglets).
    #[inline]
    fn next_ciglet_if_op(&mut self, f: impl FnOnce(u8) -> bool) -> Option<Ciglet> {
        if f(self.peek_op()?) { self.next_ciglet() } else { None }
    }

    /// Gets the last ciglet if the operation meets the specified predicate,
    /// otherwise `self` is not modified (aside from skipping empty ciglets).
    #[inline]
    fn next_ciglet_back_if_op(&mut self, f: impl FnOnce(u8) -> bool) -> Option<Ciglet> {
        if f(self.peek_op_back()?) {
            self.next_ciglet_back()
        } else {
            None
        }
    }

    /// Removes clipping from the start of the iterator.
    ///
    /// First, a hard clipping ciglet is removed if present. Then a soft
    /// clipping ciglet is removed if present. The total number of bases clipped
    /// is returned.
    #[inline]
    fn remove_clipping_front(&mut self) -> usize {
        let hard_clipping = self.next_ciglet_if_op(|op| op == b'H').map_or(0, |ciglet| ciglet.inc);
        let soft_clipping = self.next_ciglet_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
        hard_clipping + soft_clipping
    }

    /// Removes clipping from the end of the iterator.
    ///
    /// First, a hard clipping ciglet is removed if present. Then a soft
    /// clipping ciglet is removed if present. The total number of bases clipped
    /// is returned.
    #[inline]
    fn remove_clipping_back(&mut self) -> usize {
        let hard_clipping = self.next_ciglet_back_if_op(|op| op == b'H').map_or(0, |ciglet| ciglet.inc);
        let soft_clipping = self.next_ciglet_back_if_op(|op| op == b'S').map_or(0, |ciglet| ciglet.inc);
        hard_clipping + soft_clipping
    }
}

/// An extension trait for alignment-like data that enables the next ciglet to
/// be consumed and a mutable reference to it returned.
///
/// Consuming the ciglet causes the slice to be shrunk. The underlying data it
/// is refering to is not mutated, although any modifications performed to the
/// returned mutable reference will mutate it.
pub trait NextCigletMut<'a>: NextCiglet {
    /// Retrieves a mutable reference to the next [`Ciglet`] and removes it from
    /// `self`. Empty ciglets are skipped.
    fn next_ciglet_mut(&mut self) -> Option<&'a mut Ciglet>;

    /// Retrieves a mutable reference to the last [`Ciglet`] and removes it from
    /// `self`. Empty ciglets are skipped.
    fn next_ciglet_back_mut(&mut self) -> Option<&'a mut Ciglet>;

    /// Gets a mutable reference to the next ciglet if the operation meets the
    /// specified predicate, otherwise `self` is not mutated (aside from
    /// possibly removing empty ciglets).
    #[inline]
    fn next_ciglet_if_op_mut(&mut self, f: impl FnOnce(u8) -> bool) -> Option<&'a mut Ciglet> {
        if f(self.peek_op()?) { self.next_ciglet_mut() } else { None }
    }

    /// Gets a mutable reference to the last ciglet if the operation meets the
    /// specified predicate, otherwise `self` is not mutated (aside from
    /// possibly removing empty ciglets).
    #[inline]
    fn next_ciglet_back_if_op_mut(&mut self, f: impl FnOnce(u8) -> bool) -> Option<&'a mut Ciglet> {
        if f(self.peek_op_back()?) {
            self.next_ciglet_back_mut()
        } else {
            None
        }
    }
}

impl PeekOp for &[Ciglet] {
    #[inline]
    fn peek_op(&mut self) -> Option<u8> {
        self.peek_ciglet().map(|ciglet| ciglet.op)
    }

    #[inline]
    fn peek_op_back(&mut self) -> Option<u8> {
        self.peek_ciglet_back().map(|ciglet| ciglet.op)
    }
}

impl PeekCiglet for &[Ciglet] {
    fn peek_ciglet(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = self.first()?;

            if ciglet.inc == 0 {
                self.split_off_first();
            } else {
                return Some(*ciglet);
            }
        }
    }

    fn peek_ciglet_back(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = self.last()?;

            if ciglet.inc == 0 {
                self.split_off_last();
            } else {
                return Some(*ciglet);
            }
        }
    }
}

impl NextCiglet for &[Ciglet] {
    #[inline]
    fn next_ciglet(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = self.split_off_first()?;

            if ciglet.inc > 0 {
                return Some(*ciglet);
            }
        }
    }

    #[inline]
    fn next_ciglet_back(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = self.split_off_last()?;

            if ciglet.inc > 0 {
                return Some(*ciglet);
            }
        }
    }
}

impl PeekOp for &mut [Ciglet] {
    #[inline]
    fn peek_op(&mut self) -> Option<u8> {
        self.peek_ciglet().map(|ciglet| ciglet.op)
    }

    #[inline]
    fn peek_op_back(&mut self) -> Option<u8> {
        self.peek_ciglet_back().map(|ciglet| ciglet.op)
    }
}

impl TakeOp for &mut [Ciglet] {
    fn take_op(&mut self) -> Option<u8> {
        loop {
            let ciglet = self.first_mut()?;

            if let Some(new_inc) = ciglet.inc.checked_sub(1) {
                ciglet.inc = new_inc;
                return Some(ciglet.op);
            }

            self.split_off_first_mut();
        }
    }

    fn take_op_back(&mut self) -> Option<u8> {
        loop {
            let ciglet = self.last_mut()?;

            if let Some(new_inc) = ciglet.inc.checked_sub(1) {
                ciglet.inc = new_inc;
                return Some(ciglet.op);
            }

            self.split_off_last_mut();
        }
    }
}

impl PeekCiglet for &mut [Ciglet] {
    fn peek_ciglet(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = self.first()?;

            if ciglet.inc == 0 {
                self.split_off_first_mut();
            } else {
                return Some(*ciglet);
            }
        }
    }

    fn peek_ciglet_back(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = self.last()?;

            if ciglet.inc == 0 {
                self.split_off_last_mut();
            } else {
                return Some(*ciglet);
            }
        }
    }
}

impl NextCiglet for &mut [Ciglet] {
    #[inline]
    fn next_ciglet(&mut self) -> Option<Ciglet> {
        self.next_ciglet_mut().copied()
    }

    #[inline]
    fn next_ciglet_back(&mut self) -> Option<Ciglet> {
        self.next_ciglet_back_mut().copied()
    }
}

impl PeekCigletMut for &mut [Ciglet] {
    fn peek_ciglet_mut(&mut self) -> Option<&mut Ciglet> {
        let idx = if self.first()?.inc > 0 {
            0
        } else {
            self.iter().position(|ciglet| ciglet.inc > 0)?
        };

        // Range will be in bounds based on above
        let _ = self.split_off_mut(..idx);

        self.first_mut()
    }

    fn peek_ciglet_back_mut(&mut self) -> Option<&mut Ciglet> {
        let idx = if self.last()?.inc > 0 {
            self.len() - 1
        } else {
            self.iter().rposition(|ciglet| ciglet.inc > 0)?
        };

        // If range is out of bounds, then idx is already the last index
        let _ = self.split_off_mut(idx + 1..);

        self.last_mut()
    }
}

impl<'a> NextCigletMut<'a> for &'a mut [Ciglet] {
    fn next_ciglet_mut(&mut self) -> Option<&'a mut Ciglet> {
        loop {
            let ciglet = self.split_off_first_mut()?;

            if ciglet.inc > 0 {
                return Some(ciglet);
            }
        }
    }

    fn next_ciglet_back_mut(&mut self) -> Option<&'a mut Ciglet> {
        loop {
            let ciglet = self.split_off_last_mut()?;

            if ciglet.inc > 0 {
                return Some(ciglet);
            }
        }
    }
}

impl PeekOp for std::slice::Iter<'_, Ciglet> {
    fn peek_op(&mut self) -> Option<u8> {
        self.as_slice().peek_op()
    }

    fn peek_op_back(&mut self) -> Option<u8> {
        self.as_slice().peek_op_back()
    }
}

impl PeekCiglet for std::slice::Iter<'_, Ciglet> {
    fn peek_ciglet(&mut self) -> Option<Ciglet> {
        self.as_slice().peek_ciglet()
    }

    fn peek_ciglet_back(&mut self) -> Option<Ciglet> {
        self.as_slice().peek_ciglet_back()
    }
}

impl NextCiglet for std::slice::Iter<'_, Ciglet> {
    fn next_ciglet(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = *self.next()?;
            if ciglet.inc > 0 {
                return Some(ciglet);
            }
        }
    }

    fn next_ciglet_back(&mut self) -> Option<Ciglet> {
        loop {
            let ciglet = *self.next_back()?;
            if ciglet.inc > 0 {
                return Some(ciglet);
            }
        }
    }
}
