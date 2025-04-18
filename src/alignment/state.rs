use crate::{
    data::{
        cigar::CigletIterator,
        types::{
            amino_acids::AminoAcids,
            cigar::{Cigar, Ciglet},
            nucleotides::Nucleotides,
        },
    },
    private::Sealed,
};

#[derive(Clone, Debug)]
pub(crate) struct AlignmentStates(Vec<Ciglet>);

impl AlignmentStates {
    #[must_use]
    #[inline]
    pub fn new() -> Self {
        AlignmentStates(Vec::new())
    }

    pub(crate) fn add_state(&mut self, op: u8) {
        if let Some(c) = self.0.last_mut()
            && c.op == op
        {
            c.inc += 1;
        } else {
            self.0.push(Ciglet { inc: 1, op });
        }
    }

    pub(crate) fn soft_clip(&mut self, inc: usize) {
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

    #[must_use]
    pub(crate) fn to_cigar(&self) -> Cigar {
        let mut condensed = Vec::new();
        let mut format_buffer = itoa::Buffer::new();

        for Ciglet { inc, op } in self.0.iter().copied() {
            condensed.extend_from_slice(format_buffer.format(inc).as_bytes());
            condensed.push(op);
        }

        Cigar(condensed)
    }

    #[must_use]
    pub(crate) fn reverse(mut self) -> Self {
        self.0.reverse();
        self
    }

    #[allow(dead_code)]
    #[must_use]
    pub(crate) fn invert(mut self) -> Self {
        for cig in &mut self.0 {
            if cig.op == b'I' {
                cig.op = b'D';
            } else if cig.op == b'D' {
                cig.op = b'I';
            }
        }
        self
    }
}

impl Default for AlignmentStates {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

/// Aligns two sequences using a CIGAR string, inserting gaps as needed. Takes a
/// reference sequence, query sequence, CIGAR string, and a reference position,
/// and returns two aligned sequences with gaps inserted according to the CIGAR
/// operations. The first output is the reference, and the second is the query.
///
/// Only the portion of the sequences corresponding to the cigar string are
/// included in the output.
///
/// ## Panics
///
/// * All operations in the CIGAR string must be characters in `MIDNSHP=X`.
/// * The reference and query must be at least as long as the length specified
///   in the CIGAR string.
#[inline]
#[must_use]
pub fn pairwise_align_with_cigar(reference: &[u8], query: &[u8], cigar: &Cigar, ref_position: usize) -> (Vec<u8>, Vec<u8>) {
    let mut ref_index = ref_position - 1;
    let mut query_index = 0;

    let mut ref_aln = Vec::with_capacity(reference.len() + (query.len() / 2));
    let mut query_aln = Vec::with_capacity(query.len() + (reference.len() / 2));

    for Ciglet { inc, op } in cigar {
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

    /// Aligns two sequences using a CIGAR string starting at the given
    /// reference position. See [`pairwise_align_with_cigar`] for more details.
    ///
    /// ## Panics
    ///
    /// * All opcodes in the CIGAR string must be characters in `MIDNSHP=X`.
    /// * The reference and query must be at least as long as the length
    ///   specified in the CIGAR string.
    #[inline]
    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        let (r, q) = pairwise_align_with_cigar(self.as_ref(), query.as_ref(), cigar, position);
        (r.into(), q.into())
    }

    /// Returns an iterator over the bases in `query` aligned to `self`, as
    /// specified by the CIGAR string `cigar` and the starting `position`. The
    /// first base in the output is from the reference, and the second base is
    /// from the query. Gaps are represented by `None`.
    ///
    /// This has equivalent behavior as [`align_with_cigar`] but as an iterator.
    /// [`align_with_cigar`] may offer faster performance since it processes
    /// each ciglet at once, but [`align_with_cigar_iter`] avoids allocations
    /// and may offer a more convenient syntax for handling gaps.
    ///
    /// ## Panics
    ///
    /// * All opcodes in the CIGAR string must be characters in `MIDNSHP=X`.
    /// * The reference and query must be at least as long as the length
    ///   specified in the CIGAR string.
    ///
    /// [`align_with_cigar`]: PairwiseSequence::align_with_cigar
    /// [`align_with_cigar_iter`]: PairwiseSequence::align_with_cigar_iter
    #[inline]
    fn align_with_cigar_iter<'a>(&'a self, query: &'a Self, cigar: &'a Cigar, position: usize) -> AlignWithCigarIter<'a> {
        AlignWithCigarIter::new(self.as_ref(), query.as_ref(), cigar, position)
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

/// Iterator yielding the aligned bases as specified by a CIGAR string. The
/// first base is from the reference, and the second base is from the query.
/// Gaps are represented by `None`.
pub struct AlignWithCigarIter<'a> {
    reference_buffer: &'a [u8],
    query_buffer:     &'a [u8],
    ciglets:          CigletIterator<'a>,
    inc:              usize,
    op:               u8,
}

impl<'a> AlignWithCigarIter<'a> {
    /// Creates a new [`AlignWithCigarIter`] from a reference, query, cigar
    /// string, and reference position.
    #[inline]
    #[must_use]
    fn new(reference: &'a [u8], query: &'a [u8], cigar: &'a Cigar, ref_position: usize) -> Self {
        let ref_index = ref_position - 1;
        let reference_buffer = &reference[ref_index..];
        let query_buffer = query;
        let mut ciglets = cigar.into_iter();
        // If no valid ciglets, initialize with inc at 0 so that iterator is empty
        let Ciglet { inc, op } = ciglets.next().unwrap_or(Ciglet { inc: 0, op: b'M' });

        AlignWithCigarIter {
            reference_buffer,
            query_buffer,
            ciglets,
            inc,
            op,
        }
    }

    /// Get the next operation from the iterator. This will advance the Ciglet
    /// iterator if necessary. `None` is returned if the CIGAR string has
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

impl Iterator for AlignWithCigarIter<'_> {
    type Item = (Option<u8>, Option<u8>);

    /// # Panics
    ///
    /// The reference and query must be at least as long as the length specified
    /// in the CIGAR string.
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
                ref_input.align_with_cigar(&query_input, &cig, rpos),
                (ref_output, query_ouptut)
            );
        }
    }
}
