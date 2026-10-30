use crate::{alignment::AlignOp, data::cigar::Ciglet};

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
/// - All operations must be valid state, e.g, bytes in `MIDNSHP=X`.
/// - The reference and query must be at least as long as the length implied by
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
            AlignOp::M | AlignOp::Eq | AlignOp::X => {
                ref_aln.extend_from_slice(&reference[ref_index..ref_index + inc]);
                query_aln.extend_from_slice(&query[query_index..query_index + inc]);
                ref_index += inc;
                query_index += inc;
            }
            AlignOp::D => {
                ref_aln.extend_from_slice(&reference[ref_index..ref_index + inc]);
                query_aln.extend(std::iter::repeat_n(b'-', inc));
                ref_index += inc;
            }
            AlignOp::I => {
                ref_aln.extend(std::iter::repeat_n(b'-', inc));
                query_aln.extend_from_slice(&query[query_index..query_index + inc]);
                query_index += inc;
            }
            AlignOp::S => query_index += inc,
            AlignOp::N => {
                ref_aln.extend_from_slice(&reference[ref_index..ref_index + inc]);
                query_aln.extend(std::iter::repeat_n(b'N', inc));
                ref_index += inc;
            }
            AlignOp::H | AlignOp::P => {}
        }
    }

    (ref_aln, query_aln)
}
