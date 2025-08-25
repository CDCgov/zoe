use crate::{
    alignment::Alignment,
    data::{
        cigar::{CigarError, CigletIteratorChecked},
        types::{
            amino_acids::AminoAcids,
            cigar::{Cigar, Ciglet},
            nucleotides::Nucleotides,
        },
    },
    private::Sealed,
};
use std::{
    fmt::Write,
    mem::MaybeUninit,
    simd::{LaneCount, SupportedLaneCount, prelude::*},
};

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

/// A struct for storing alignment states, which can be converted into a
/// [`Cigar`] string.
#[derive(Clone, PartialEq, Eq)]
pub struct AlignmentStates(pub(crate) Vec<Ciglet>);

impl AlignmentStates {
    /// Initializes an empty alignment.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        AlignmentStates(Vec::new())
    }

    /// Initializes the states with capacity for `n` (inc, op) pairs.
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
    pub fn add_state(&mut self, op: u8) {
        self.add_ciglet(Ciglet { inc: 1, op });
    }

    /// Adds a ciglet to the right end of the alignment. If the operation is the
    /// same as the rightmost operation, the ciglet is merged with the last one.
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

    /// Adds an alignment `op` of size `inc`.
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

    pub(crate) fn extend_from_ciglets<I>(&mut self, ciglets: I)
    where
        I: IntoIterator<Item = Ciglet>, {
        let mut ciglets = ciglets.into_iter();
        let Some(first_ciglet) = ciglets.next() else { return };
        self.add_ciglet(first_ciglet);
        self.0.extend(ciglets);
    }

    /// Adds soft clipping `S` to the end of the alignment `inc` times.
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

    /// Converts the [`AlignmentStates`] struct to a [`Cigar`] string. All
    /// operations should be valid for a CIGAR string.
    #[must_use]
    pub fn to_cigar_unchecked(&self) -> Cigar {
        Cigar::from_ciglets_unchecked(self.0.iter().copied())
    }

    /// Reverses the order of the stored alignment states in-place.
    #[inline]
    pub fn make_reverse(&mut self) {
        self.0.reverse();
    }

    /// Reverses the order of the stored alignment states in-place.
    #[inline]
    #[must_use]
    pub fn to_reverse(&self) -> Self {
        self.into_iter().rev().collect()
    }

    /// Yields an iterator over the alignment states.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, Ciglet> {
        self.0.iter()
    }
}

impl Default for AlignmentStates {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl PartialEq<Cigar> for AlignmentStates {
    #[inline]
    fn eq(&self, other: &Cigar) -> bool {
        let mut o = other.iter();
        let matches = self.iter().copied().eq(o.by_ref());
        matches && o.valid()
    }
}

impl<'a> IntoIterator for &'a AlignmentStates {
    type Item = Ciglet;
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, Ciglet>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter().copied()
    }
}

impl IntoIterator for AlignmentStates {
    type Item = Ciglet;
    type IntoIter = <Vec<Ciglet> as IntoIterator>::IntoIter;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl FromIterator<Ciglet> for AlignmentStates {
    #[inline]
    fn from_iter<T: IntoIterator<Item = Ciglet>>(iter: T) -> Self {
        AlignmentStates(iter.into_iter().collect::<Vec<_>>())
    }
}

impl std::fmt::Debug for AlignmentStates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.0.iter().map(|c| (c.inc, c.op as char))).finish()
    }
}

impl std::fmt::Display for AlignmentStates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buff = itoa::Buffer::new();
        for Ciglet { inc, op } in self {
            f.write_str(buff.format(inc))?;
            f.write_char(op as char)?;
        }
        Ok(())
    }
}

impl TryFrom<&[u8]> for AlignmentStates {
    type Error = CigarError;

    fn try_from(v: &[u8]) -> Result<Self, Self::Error> {
        if v == b"*" {
            Ok(AlignmentStates::new())
        } else {
            CigletIteratorChecked::new(v)
                .collect::<Result<Vec<_>, _>>()
                .map(AlignmentStates)
        }
    }
}

impl TryFrom<String> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(s: String) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(s.as_bytes())
    }
}

impl TryFrom<Vec<u8>> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(bytes: Vec<u8>) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(bytes.as_slice())
    }
}

impl<const N: usize> TryFrom<[u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: [u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl TryFrom<&mut [u8]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(&*v)
    }
}

impl<const N: usize> TryFrom<&[u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &[u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl<const N: usize> TryFrom<&mut [u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl TryFrom<&str> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(s.as_bytes())
    }
}

/// Produces aligned sequences using an iterator of [`Ciglet`] elements,
/// inserting gaps as needed.
///
/// This also takes a reference sequence, query sequence, and a reference index.
/// The first output is the reference, and the second is the query.
///
/// Only the portion of the sequences corresponding to the ciglets are included
/// in the output.
///
/// ## Panics
///
/// * All operations must be valid state, e.g, bytes in `MIDNSHP=X`.
/// * The reference and query must be at least as long as the length implied by
///   alignment operations.
#[inline]
#[must_use]
pub fn pairwise_align_with(
    reference: &[u8], query: &[u8], ciglets: impl IntoIterator<Item = Ciglet>, mut ref_index: usize,
) -> (Vec<u8>, Vec<u8>) {
    let mut query_index = 0;

    let mut ref_aln = Vec::with_capacity(reference.len() + (query.len() / 2));
    let mut query_aln = Vec::with_capacity(query.len() + (reference.len() / 2));

    for Ciglet { inc, op } in ciglets {
        match op {
            b'M' | b'=' | b'X' => {
                ref_aln.extend_from_slice(&reference[ref_index..ref_index + inc]);
                query_aln.extend_from_slice(&query[query_index..query_index + inc]);
                ref_index += inc;
                query_index += inc;
            }
            b'D' => {
                ref_aln.extend_from_slice(&reference[ref_index..ref_index + inc]);
                query_aln.extend(std::iter::repeat_n(b'-', inc));
                ref_index += inc;
            }
            b'I' => {
                ref_aln.extend(std::iter::repeat_n(b'-', inc));
                query_aln.extend_from_slice(&query[query_index..query_index + inc]);
                query_index += inc;
            }
            b'S' => query_index += inc,
            b'N' => {
                ref_aln.extend_from_slice(&reference[ref_index..ref_index + inc]);
                query_aln.extend(std::iter::repeat_n(b'N', inc));
                ref_index += inc;
            }
            b'H' | b'P' => {}
            // TODO: Could ignore if we allow only valid CIGAR by default.
            _ => panic!("CIGAR op '{op}' not supported.\n"),
        }
    }

    (ref_aln, query_aln)
}

/// Enables sequence expansion based on alignment information.
pub trait PairwiseSequence: AsRef<[u8]> + Sealed {
    type Output: From<Vec<u8>>;

    /// Aligns two sequences using a [`Ciglet`] iterator starting at the given
    /// reference position. See [`pairwise_align_with`] for more details.
    ///
    /// ## Panics
    ///
    /// * All opcodes in the CIGAR must be characters in `MIDNSHP=X`.
    /// * The reference and query must be at least as long as the length
    ///   specified in the alignment states.
    #[inline]
    fn align_and_collect(
        &self, query: &Self, ciglets: impl IntoIterator<Item = Ciglet>, position: usize,
    ) -> (Self::Output, Self::Output) {
        let (r, q) = pairwise_align_with(self.as_ref(), query.as_ref(), ciglets, position - 1);
        (r.into(), q.into())
    }

    /// Returns an iterator over the bases in `query` aligned to `self`, as
    /// specified by the [`Ciglet`] iterator and the starting `position`. The
    /// first base in the output is from the reference, and the second base is
    /// from the query. Gaps are represented by `None`.
    ///
    /// This has equivalent behavior as [`align_and_collect`] but as an
    /// iterator. [`align_and_collect`] may offer faster performance since it
    /// processes each ciglet at once, but [`align_iter`] avoids allocations and
    /// may offer a more convenient syntax for handling gaps.
    ///
    /// ## Panics
    ///
    /// * All opcodes must be characters in `MIDNSHP=X`.
    /// * The reference and query must be at least as long as the length
    ///   specified in the alignment operations.
    ///
    /// [`align_and_collect`]: PairwiseSequence::align_and_collect
    /// [`align_iter`]: PairwiseSequence::align_iter
    #[inline]
    fn align_iter<'a, I>(&'a self, query: &'a Self, ciglets: I, position: usize) -> AlignmentIter<'a, I>
    where
        I: Iterator<Item = Ciglet>, {
        AlignmentIter::new(self.as_ref(), query.as_ref(), ciglets, position - 1)
    }
}

impl PairwiseSequence for Vec<u8> {
    type Output = Self;
}

impl PairwiseSequence for &[u8] {
    type Output = Vec<u8>;
}

impl PairwiseSequence for Nucleotides {
    type Output = Self;
}

impl PairwiseSequence for AminoAcids {
    type Output = Self;
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
    /// in the iterator's alignment operations.
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

/// Experimental backtracking matrix for alignment with affine gap scores.
#[derive(Clone, Debug)]
pub(crate) struct BacktrackMatrix {
    pub data: Vec<u8>,
    cols:     usize,
    cursor:   usize,
}

impl BacktrackMatrix {
    const UP: u8 = 1;
    const UP_EXTENDING: u8 = 2;
    const LEFT: u8 = 4;
    const LEFT_EXTENDING: u8 = 8;
    const STOP: u8 = 16;

    pub(crate) fn new(rows: usize, cols: usize) -> Self {
        BacktrackMatrix {
            data: vec![0u8; rows * cols],
            cols,
            cursor: 0,
        }
    }

    #[inline]
    pub(crate) fn up(&mut self) {
        self.data[self.cursor] |= Self::UP;
    }

    #[inline]
    pub(crate) fn left(&mut self) {
        self.data[self.cursor] |= Self::LEFT;
    }

    #[inline]
    pub(crate) fn up_extending(&mut self) {
        self.data[self.cursor] |= Self::UP_EXTENDING;
    }

    #[inline]
    pub(crate) fn left_extending(&mut self) {
        self.data[self.cursor] |= Self::LEFT_EXTENDING;
    }

    #[inline]
    pub(crate) fn stop(&mut self) {
        self.data[self.cursor] = Self::STOP;
    }
}

#[derive(Clone, Debug)]
pub(crate) struct BacktrackMatrixStriped<const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    pub data:  Vec<Simd<u8, N>>,
    num_vecs:  usize,
    v_cursor:  usize,
    curr_lane: usize,
}

impl<const N: usize> BacktrackMatrixStriped<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Allocates lazily-initialized backtrack data which can later be converted
    /// to a [`BacktrackMatrixStriped`].
    ///
    /// The total number of SIMD vectors is given by `size`.
    #[inline]
    #[must_use]
    pub(crate) fn make_uninit_data(size: usize) -> Vec<MaybeUninit<Simd<u8, N>>> {
        let mut data = Vec::with_capacity(size);
        data.resize_with(size, MaybeUninit::uninit);
        data
    }

    /// Wraps a vector of SIMD vectors in a [`BacktrackMatrixStriped`] to
    /// facilitate backtracking.
    #[inline]
    #[must_use]
    pub(crate) fn new(data: Vec<Simd<u8, N>>, num_vecs: usize) -> Self {
        assert!(data.len().is_multiple_of(num_vecs));
        BacktrackMatrixStriped {
            data,
            num_vecs,
            v_cursor: 0,
            curr_lane: 0,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn debug_cell(&self, r: usize, c: usize) -> String {
        let v = c % self.num_vecs;
        let lane = (c - v) / self.num_vecs;
        let v = self.data[self.num_vecs * r + v];
        v.debug_cell(lane)
    }

    #[allow(dead_code)]
    pub(crate) fn print_row(&self, r: usize) {
        let mut bt = self.clone();
        let nc = self.num_vecs * N;

        print!("{r:02}: x");
        for c in 0..nc {
            bt.move_to(r, c);
            if bt.is_stop() {
                print!("o");
            } else if bt.is_up() {
                print!("^");
            } else if bt.is_left() {
                print!("<");
            } else if bt.is_up_extending() {
                print!(":");
            } else if bt.is_left_extending() {
                print!("-");
            } else {
                print!("\\");
            }
        }
        println!();
    }
}

/// Trait for SIMD backtracking operations on packed u8 values
pub(crate) trait SimdBacktrackFlags<const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    const UP: Simd<u8, N>;
    const UP_EXTENDING: Simd<u8, N>;
    const LEFT: Simd<u8, N>;
    const LEFT_EXTENDING: Simd<u8, N>;
    const STOP: Simd<u8, N>;

    /// Creates a new match variable
    fn simd_match() -> Simd<u8, N> {
        Simd::splat(0)
    }

    /// Add Up direction to lanes selected by mask
    fn simd_up(&mut self, mask: Mask<i8, N>);

    /// Add Left direction to lanes selected by mask
    fn simd_left(&mut self, mask: Mask<i8, N>);

    fn simd_correct_and_set_left(&mut self, mask: Mask<i8, N>);

    /// Apply Up Extending direction to lanes selected by mask
    fn simd_up_extending(&mut self, mask: Mask<i8, N>);

    /// Add Left Extending direction to lanes selected by mask
    fn simd_left_extending(&mut self, mask: Mask<i8, N>);

    /// Set Stop state to lanes selected by mask
    fn simd_stop(&mut self, mask: Mask<i8, N>);

    /// Debug representation of a cell at the given lane
    fn debug_cell(&self, lane: usize) -> String;
}

impl<const N: usize> SimdBacktrackFlags<N> for Simd<u8, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    const UP: Simd<u8, N> = Simd::splat(BacktrackMatrix::UP);
    const UP_EXTENDING: Simd<u8, N> = Simd::splat(BacktrackMatrix::UP_EXTENDING);
    const LEFT: Simd<u8, N> = Simd::splat(BacktrackMatrix::LEFT);
    const LEFT_EXTENDING: Simd<u8, N> = Simd::splat(BacktrackMatrix::LEFT_EXTENDING);
    const STOP: Simd<u8, N> = Simd::splat(BacktrackMatrix::STOP);

    #[inline]
    fn simd_up(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(*self | Self::UP, *self);
    }

    #[inline]
    fn simd_left(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(*self | Self::LEFT, *self);
    }

    #[inline]
    fn simd_correct_and_set_left(&mut self, mask: Mask<i8, N>) {
        *self = mask.select((*self & Self::UP_EXTENDING) | Self::LEFT, *self);
    }

    #[inline]
    fn simd_up_extending(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(*self | Self::UP_EXTENDING, *self);
    }

    #[inline]
    fn simd_left_extending(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(*self | Self::LEFT_EXTENDING, *self);
    }

    #[inline]
    fn simd_stop(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(Self::STOP, *self);
    }

    fn debug_cell(&self, lane: usize) -> String {
        let cell = self[lane];

        let mut states = String::new();
        if cell & BacktrackMatrix::STOP > 0 {
            states.push('o');
        }

        if cell & BacktrackMatrix::UP > 0 {
            states.push('^');
        }

        if cell & BacktrackMatrix::LEFT > 0 {
            states.push('<');
        }

        if cell & BacktrackMatrix::UP_EXTENDING > 0 {
            states.push(':');
        }

        if cell & BacktrackMatrix::LEFT_EXTENDING > 0 {
            states.push('-');
        }

        if cell == 0 {
            states.push('\\');
        }

        states
    }
}

/// Trait for backtracking operations on a backtracking matrix
pub(crate) trait BackTrackable {
    /// Moves the cursor to the specified row and column
    fn move_to(&mut self, r: usize, c: usize);

    /// Checks if the current cell is in the **up** state
    fn is_up(&self) -> bool;

    /// Checks if the current cell is in the **left** state
    fn is_left(&self) -> bool;

    /// Checks if the current cell is in the **left extending** state
    fn is_left_extending(&self) -> bool;

    /// Checks if the current cell is in the **up extending** state
    fn is_up_extending(&self) -> bool;

    /// Checks if the current cell is in the **stop** state
    fn is_stop(&self) -> bool;

    /// Creates an [`Alignment`] from a backtrack matrix and other information.
    fn to_alignment<T>(
        &mut self, score: T, mut r_end: usize, mut c_end: usize, ref_len: usize, query_len: usize,
    ) -> Alignment<T> {
        let mut states = AlignmentStates::new();
        let mut op = 0;

        self.move_to(r_end, c_end); // 0-based move to max

        // 1-based as though we had a padded matrix
        r_end += 1;
        c_end += 1;

        let (mut r, mut c) = (r_end, c_end);

        // soft clip 3'
        states.soft_clip(query_len - c);

        while !self.is_stop() && r > 0 && c > 0 {
            if op == b'D' && self.is_up_extending() {
                op = b'D';
                r -= 1;
            } else if op == b'I' && self.is_left_extending() {
                op = b'I';
                c -= 1;
            } else if self.is_up() {
                op = b'D';
                r -= 1;
            } else if self.is_left() {
                op = b'I';
                c -= 1;
            } else {
                op = b'M';
                r -= 1;
                c -= 1;
            }
            states.add_state(op);
            self.move_to(r.saturating_sub(1), c.saturating_sub(1)); // 0-based next position
        }

        // soft clip 5'
        states.soft_clip(c);
        states.make_reverse();

        // r and c are decremented and becomes 0-based
        Alignment {
            score,
            ref_range: r..r_end,
            query_range: c..c_end,
            states,
            ref_len,
            query_len,
        }
    }

    /// Prints the backtracking matrix in a human-readable format
    #[allow(dead_code)]
    fn print(&self);
}

impl BackTrackable for BacktrackMatrix {
    #[inline]
    fn move_to(&mut self, r: usize, c: usize) {
        self.cursor = self.cols * r + c;
    }

    #[inline]
    fn is_up(&self) -> bool {
        self.data[self.cursor] & Self::UP > 0
    }

    #[inline]
    fn is_left(&self) -> bool {
        self.data[self.cursor] & Self::LEFT > 0
    }

    #[inline]
    fn is_left_extending(&self) -> bool {
        self.data[self.cursor] & Self::LEFT_EXTENDING > 0
    }

    #[inline]
    fn is_up_extending(&self) -> bool {
        self.data[self.cursor] & Self::UP_EXTENDING > 0
    }

    #[inline]
    fn is_stop(&self) -> bool {
        self.data[self.cursor] & Self::STOP > 0
    }

    fn print(&self) {
        let mut bt = self.clone();
        let nr = self.data.len() / self.cols;
        let nc = self.cols;

        print!("    ");
        for _ in 0..=nc {
            print!("x");
        }
        println!();
        for r in 0..nr {
            print!("{r:02}: x");
            for c in 0..nc {
                bt.move_to(r, c);
                if bt.is_stop() {
                    print!("o");
                } else if bt.is_up() {
                    print!("^");
                } else if bt.is_left() {
                    print!("<");
                } else if bt.is_up_extending() {
                    print!(":");
                } else if bt.is_left_extending() {
                    print!("-");
                } else {
                    print!("\\");
                }
            }
            println!();
        }
    }
}

impl<const N: usize> BackTrackable for BacktrackMatrixStriped<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn move_to(&mut self, r: usize, c: usize) {
        let v = c % self.num_vecs;
        self.curr_lane = (c - v) / self.num_vecs;
        self.v_cursor = self.num_vecs * r + v;
    }

    #[inline]
    fn is_up(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::UP) > 0
    }

    #[inline]
    fn is_left(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::LEFT) > 0
    }

    #[inline]
    fn is_left_extending(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::LEFT_EXTENDING) > 0
    }

    #[inline]
    fn is_up_extending(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::UP_EXTENDING) > 0
    }

    #[inline]
    fn is_stop(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::STOP) > 0
    }

    fn print(&self) {
        let mut bt = self.clone();
        let nr = self.data.len() / self.num_vecs;
        let nc = self.num_vecs * N;

        print!("    ");
        for _ in 0..=nc {
            print!("x");
        }
        println!();
        for r in 0..nr {
            print!("{r:02}: x");
            for c in 0..nc {
                bt.move_to(r, c);
                if bt.is_stop() {
                    print!("o");
                } else if bt.is_up() {
                    print!("^");
                } else if bt.is_left() {
                    print!("<");
                } else if bt.is_up_extending() {
                    print!(":");
                } else if bt.is_left_extending() {
                    print!("-");
                } else {
                    print!("\\");
                }
            }
            println!();
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn states_sequence() {
        // The clipping is not valid here, but we are using it to exercise all
        // options
        let cigar = Cigar::try_from("3H4S4S3H10M9I3M1D5M3H5S10S1H").unwrap();
        let alignment_states = cigar.iter().collect::<AlignmentStates>();
        let mut cigar_slice = alignment_states.as_slice();
        let mut ciglet_iterator = cigar.into_iter();

        assert_eq!(cigar_slice.peek_op(), Some(b'H'));
        assert_eq!(ciglet_iterator.peek_op(), Some(b'H'));
        assert_eq!(cigar_slice.peek_back_op(), Some(b'H'));
        assert_eq!(ciglet_iterator.peek_back_op(), Some(b'H'));

        // Shouldn't have changed, since we didn't advance
        assert_eq!(cigar_slice.peek_op(), Some(b'H'));
        assert_eq!(ciglet_iterator.peek_op(), Some(b'H'));
        assert_eq!(cigar_slice.peek_back_op(), Some(b'H'));
        assert_eq!(ciglet_iterator.peek_back_op(), Some(b'H'));

        assert_eq!(cigar_slice.remove_clipping_front(), 7);
        assert_eq!(ciglet_iterator.remove_clipping_front(), 7);
        assert_eq!(cigar_slice.remove_clipping_back(), 11);
        assert_eq!(ciglet_iterator.remove_clipping_back(), 11);

        assert_eq!(cigar_slice.peek_op(), Some(b'S'));
        assert_eq!(ciglet_iterator.peek_op(), Some(b'S'));
        assert_eq!(cigar_slice.peek_back_op(), Some(b'S'));
        assert_eq!(ciglet_iterator.peek_back_op(), Some(b'S'));

        assert_eq!(cigar_slice.remove_clipping_front(), 4);
        assert_eq!(ciglet_iterator.remove_clipping_front(), 4);
        assert_eq!(cigar_slice.remove_clipping_back(), 5);
        assert_eq!(ciglet_iterator.remove_clipping_back(), 5);

        assert_eq!(cigar_slice.remove_clipping_front(), 3);
        assert_eq!(ciglet_iterator.remove_clipping_front(), 3);
        assert_eq!(cigar_slice.remove_clipping_back(), 3);
        assert_eq!(ciglet_iterator.remove_clipping_back(), 3);

        assert_eq!(cigar_slice.remove_clipping_front(), 0);
        assert_eq!(ciglet_iterator.remove_clipping_front(), 0);
        assert_eq!(cigar_slice.remove_clipping_back(), 0);
        assert_eq!(ciglet_iterator.remove_clipping_back(), 0);

        assert_eq!(cigar_slice.next_if_op(|op| op == b'I'), None);
        assert_eq!(ciglet_iterator.next_if_op(|op| op == b'I'), None);
        assert_eq!(cigar_slice.next_back_if_op(|op| op == b'I'), None);
        assert_eq!(ciglet_iterator.next_back_if_op(|op| op == b'I'), None);

        assert_eq!(cigar_slice.next_if_op(|op| op == b'M'), Some(Ciglet { inc: 10, op: b'M' }));
        assert_eq!(
            ciglet_iterator.next_if_op(|op| op == b'M'),
            Some(Ciglet { inc: 10, op: b'M' })
        );
        assert_eq!(
            cigar_slice.next_back_if_op(|op| op == b'M'),
            Some(Ciglet { inc: 5, op: b'M' })
        );
        assert_eq!(
            ciglet_iterator.next_back_if_op(|op| op == b'M'),
            Some(Ciglet { inc: 5, op: b'M' })
        );

        assert_eq!(cigar_slice.next_ciglet(), Some(Ciglet { inc: 9, op: b'I' }));
        assert_eq!(ciglet_iterator.next_ciglet(), Some(Ciglet { inc: 9, op: b'I' }));
        assert_eq!(cigar_slice.next_ciglet_back(), Some(Ciglet { inc: 1, op: b'D' }));
        assert_eq!(ciglet_iterator.next_ciglet_back(), Some(Ciglet { inc: 1, op: b'D' }));

        assert!(!cigar_slice.is_empty());
        assert!(!ciglet_iterator.is_empty());

        assert_eq!(cigar_slice.next_ciglet(), Some(Ciglet { inc: 3, op: b'M' }));
        assert_eq!(ciglet_iterator.next_ciglet(), Some(Ciglet { inc: 3, op: b'M' }));
        assert_eq!(cigar_slice.next_ciglet_back(), None);
        assert_eq!(ciglet_iterator.next_ciglet_back(), None);
        assert_eq!(cigar_slice.next_ciglet(), None);
        assert_eq!(ciglet_iterator.next_ciglet(), None);

        assert!(cigar_slice.is_empty());
        assert!(ciglet_iterator.is_empty());
    }

    #[test]
    fn align_with_cigar() {
        let data: [(usize, Cigar, [AminoAcids; 4]); 1] = [(
            2,
            Cigar::from_slice_unchecked("4M"),
            [b"PLEASANTLY".into(), b"MEANLY".into(), b"LEAS".into(), b"MEAN".into()],
        )];

        for (rpos, cig, [ref_input, query_input, ref_output, query_ouptut]) in data {
            assert_eq!(
                ref_input.align_and_collect(&query_input, &cig, rpos),
                (ref_output, query_ouptut)
            );
        }
    }

    #[test]
    fn test_eq_states_cigar() {
        let mut states = AlignmentStates::new();
        states.add_ciglet(Ciglet { inc: 10, op: b'M' });

        let cigar = Cigar::from_slice_unchecked(b"10M3D");
        assert_ne!(states, cigar);
        assert_ne!(cigar, states);

        let cigar = Cigar::from_slice_unchecked(b"10M$$");
        assert_ne!(states, cigar);
        assert_ne!(cigar, states);

        let cigar = Cigar::from_slice_unchecked(b"10M");
        assert_eq!(states, cigar);
        assert_eq!(cigar, states);

        let cigar = Cigar::new();
        let states = AlignmentStates::new();
        assert_eq!(states, cigar);
        assert_eq!(cigar, states);
    }
}
