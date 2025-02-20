use crate::data::types::{
    amino_acids::AminoAcids,
    cigar::{Cigar, Ciglet},
    nucleotides::Nucleotides,
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
    fn default() -> Self {
        Self::new()
    }
}

/// Aligns two sequences using a CIGAR string, inserting gaps as needed. Takes a
/// reference sequence, query sequence, CIGAR stringl reference position, and
/// returns two aligned sequences with gaps inserted according to the CIGAR
/// operations.
///
/// ## Panics
///
/// Panics if an invalid cigar string character is provided. Valid characters
/// are: MIDNSHP=X
#[inline]
#[must_use]
pub fn pairwise_align_with_cigar(reference: &[u8], query: &[u8], cigar: &Cigar, ref_position: usize) -> (Vec<u8>, Vec<u8>) {
    let mut ref_index = ref_position - 1;
    let mut query_index = 0;

    let mut ref_aln: Vec<u8> = Vec::with_capacity(reference.len() + (query.len() / 2));
    let mut query_aln: Vec<u8> = Vec::with_capacity(query.len() + (reference.len() / 2));

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
pub trait PairwiseSequence {
    type Output;

    /// Aligns two sequences using a CIGAR string starting at the
    /// given reference position.
    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output);
}

impl PairwiseSequence for Vec<u8> {
    type Output = Self;

    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        pairwise_align_with_cigar(self, query, cigar, position)
    }
}

impl PairwiseSequence for &[u8] {
    type Output = Vec<u8>;

    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        pairwise_align_with_cigar(self, query, cigar, position)
    }
}

impl PairwiseSequence for Nucleotides {
    type Output = Self;

    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        let (r, q) = pairwise_align_with_cigar(self.as_ref(), query.as_ref(), cigar, position);
        (r.into(), q.into())
    }
}

impl PairwiseSequence for AminoAcids {
    type Output = Self;

    fn align_with_cigar(&self, query: &Self, cigar: &Cigar, position: usize) -> (Self::Output, Self::Output) {
        let (r, q) = pairwise_align_with_cigar(self.as_ref(), query.as_ref(), cigar, position);
        (r.into(), q.into())
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
            "4M".into(),
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
