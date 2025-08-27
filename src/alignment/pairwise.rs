use crate::{
    alignment::AlignmentIter,
    data::{
        cigar::Ciglet,
        types::{amino_acids::AminoAcids, nucleotides::Nucleotides},
    },
    private::Sealed,
};

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
