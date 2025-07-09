use crate::{
    data::types::{
        amino_acids::AminoAcids,
        cigar::{Cigar, Ciglet},
        nucleotides::Nucleotides,
    },
    private::Sealed,
};

/// A struct for storing alignment states and incrementally building a [`Cigar`]
/// string.
#[derive(Clone, PartialEq, Eq)]
pub struct AlignmentStates(pub(crate) Vec<Ciglet>);

impl AlignmentStates {
    /// Initializes an empty alignment
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        AlignmentStates(Vec::new())
    }

    /// Adds a state to the right end of the alignment
    pub fn add_state(&mut self, op: u8) {
        self.add_ciglet(Ciglet { inc: 1, op });
    }

    /// Adds a ciglet to the right end of the alignment. If the operation is the
    /// same as the rightmost operation, the ciglet is merged with the last one.
    pub fn add_ciglet(&mut self, ciglet: Ciglet) {
        if let Some(c) = self.0.last_mut()
            && c.op == ciglet.op
        {
            c.inc += 1;
        } else {
            self.0.push(ciglet);
        }
    }

    /// Adds ciglets to the alignment. It is assumed that consecutive ciglets
    /// have different operations.
    pub(crate) fn extend_from_ciglets<I>(&mut self, ciglets: I)
    where
        I: IntoIterator<Item = Ciglet>, {
        let mut ciglets = ciglets.into_iter();
        let Some(first_ciglet) = ciglets.next() else { return };
        self.add_ciglet(first_ciglet);
        self.0.extend(ciglets);
    }

    /// Adds soft clipping `S` to the end of the alignment `inc` times
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

    /// Converts the [`AlignmentStates`] struct to a [`Cigar`] string. All
    /// operations should be valid for a CIGAR string.
    #[must_use]
    pub fn to_cigar_unchecked(&self) -> Cigar {
        Cigar::from_ciglets_unchecked(self.0.iter().copied())
    }

    /// Reverses the order of the stored alignment states
    #[inline]
    pub fn reverse(&mut self) -> &mut Self {
        self.0.reverse();
        self
    }

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

impl std::fmt::Debug for AlignmentStates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.0.iter().map(|c| (c.inc, c.op as char))).finish()
    }
}

/// Aligns two sequences using anything that can be converted to an iterator
/// over [`Ciglet`], inserting gaps as needed. Also takes a reference sequence,
/// query sequence, and a reference index, and returns two aligned sequences
/// with gaps inserted according to the alignment state operations. The first
/// output is the reference, and the second is the query.
///
/// Only the portion of the sequences corresponding to the [`Cigar`] /
/// [`AlignmentStates`] are included in the output.
///
/// ## Panics
///
/// * All operations must be valid state, e.g, bytes in `MIDNSHP=X`.
/// * The reference and query must be at least as long as the length implied
///   by alignment operations.
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
/// operations. The first base is from the reference, and the second base is
/// from the query. Gaps are represented by `None`.
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

/// Experimental score cell for alignment.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct ScoreCell {
    pub endgap_up: i32,
    pub matching:  i32,
}

/// Experimental scoring struct for alignment.
pub(crate) struct BestScore {
    row:   usize,
    col:   usize,
    score: i32,
}

impl BestScore {
    pub fn new() -> Self {
        Self {
            row:   0,
            col:   0,
            score: 0,
        }
    }

    #[inline]
    pub fn add_score(&mut self, r: usize, c: usize, new_score: i32) {
        if new_score >= self.score {
            self.row = r;
            self.col = c;
            self.score = new_score;
        }
    }

    pub fn get_best_score(&self) -> (usize, usize, i32) {
        (self.row, self.col, self.score)
    }
}

impl Default for BestScore {
    fn default() -> Self {
        Self::new()
    }
}

/// Experimental backtracking matrix for alignment with affine gap scores.
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
    pub(crate) fn move_to(&mut self, r: usize, c: usize) {
        self.cursor = self.cols * r + c;
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

    #[inline]
    pub(crate) fn is_up(&self) -> bool {
        self.data[self.cursor] & Self::UP > 0
    }

    #[inline]
    pub(crate) fn is_left(&self) -> bool {
        self.data[self.cursor] & Self::LEFT > 0
    }

    #[inline]
    pub(crate) fn is_left_extending(&self) -> bool {
        self.data[self.cursor] & Self::LEFT_EXTENDING > 0
    }

    #[inline]
    pub(crate) fn is_up_extending(&self) -> bool {
        self.data[self.cursor] & Self::UP_EXTENDING > 0
    }

    #[inline]
    pub(crate) fn is_stop(&self) -> bool {
        self.data[self.cursor] & Self::STOP > 0
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::data::types::cigar::Cigar;

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
