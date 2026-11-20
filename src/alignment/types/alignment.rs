//! The [`Alignment`] type and its methods.

use crate::{
    alignment::{AlignOp, AlignmentIter, AlignmentStates, ByteSeq, pairwise_align_with},
    data::{cigar::Ciglet, views::IndexAdjustable},
};
use std::ops::Range;

// For the `Alignment` struct below, both ranges are 0-based and end-exclusive.
// For a global alignment, the ranges will each be encompass the full length of
// the sequences. For local alignment, `query_range` does NOT include hard or
// soft clipped bases, even though `states` does contain this information.

/// A struct representing the information for an alignment, such as its score
/// and where in the sequences it occurs.
///
/// ## Validity and Panics
///
/// [`Alignment`] contains many fields that are related to each other. The
/// fields should generally not be mutated or overwritten, unless care is taken
/// to ensure all state remains consistent. API involving [`Alignment`] may
/// panic if assumptions are invalidated. The particular requirements are:
///
/// - `ref_range` must have `start <= end` and `end <= ref_len`
/// - `query_range` must have `start <= end` and `end <= query_len`
/// - The length of `ref_range` must equal [`ref_len_in_alignment`] called on
///   `states`
/// - The length of `query_range` must equal [`query_len_in_alignment`] called
///   on `states`, excluding soft clipping
/// - Only hard clipping is allowed to exist between soft clipping and an end of
///   `states`
///
/// ## Conversion to SAM format
///
/// To convert this to a [`SamData`] record, see [`SamData::from_alignment`].
/// Some additional information will need to be provided, such as any flags, the
/// query sequence and quality scores if present, etc.
///
/// ## Reversibility
///
/// [`Alignment`] and `MaybeAligned<Alignment<_>>` both offer methods for
/// reversing an alignment. For example, given an alignment between `query` and
/// `reference`, calling [`Alignment::to_reverse`] will return a new optimal
/// alignment between `query.to_reverse()` and `reference.to_reverse()`, without
/// recomputing the alignment. [`Alignment::make_reverse`] does the same thing
/// in-place.
///
/// Notably, since there may be many optimal alignments with the same score,
/// this is not guaranteed to equal the same alignment that an alignment routine
/// would return when provided the reversed sequences.
///
/// This can also be interpreted as an alignment between
/// `query.to_reverse_complement()` and `reference.to_reverse_complement()` if
/// the [`WeightMatrix`] is compatible with this. For example, the score for
/// `A→A` must be the same as `T→T`, `A→G` must be the same as `T→C`, and so on.
/// This property holds for simple weight matrices with a fixed match and
/// mismatch score, but it also holds for many other custom matrices.
///
/// One use-case for these functions is when the reverse complements of the
/// references are precomputed. In this case, it may be necessary to reverse the
/// alignment, so that it corresponds to the alignment between a
/// reverse-complemented query and the original reference.
///
/// [`SamData`]: crate::data::records::sam::SamData
/// [`SamData::from_alignment`]:
///     crate::data::records::sam::SamData::from_alignment
/// [`WeightMatrix`]: crate::data::WeightMatrix
/// [`ref_len_in_alignment`]:
///     crate::data::types::cigar::LenInAlignment::ref_len_in_alignment
/// [`query_len_in_alignment`]:
///     crate::data::types::cigar::LenInAlignment::query_len_in_alignment
#[non_exhaustive]
#[derive(Clone, Eq, PartialEq, Debug, Default)]
pub struct Alignment<T> {
    /// The score of the alignment
    pub score:       T,
    /// A range for the indices of the reference that are included in the
    /// alignment
    pub ref_range:   Range<usize>,
    /// A range for the indices of the query that are included in the alignment,
    /// already **excludes** clipped bases
    pub query_range: Range<usize>,
    /// States mapping reference and query.
    pub states:      AlignmentStates,
    /// The length of the reference sequence
    pub ref_len:     usize,
    /// The length of the query sequence
    pub query_len:   usize,
}

impl<T> Alignment<T> {
    /// Returns an [`Alignment`] struct representing a global alignment.
    ///
    /// This sets `ref_range` to `0..ref_len` and `query_range` to
    /// `0..query_len`.
    #[inline]
    #[must_use]
    pub fn new_global(score: T, states: AlignmentStates, ref_len: usize, query_len: usize) -> Self {
        Self {
            score,
            ref_range: 0..ref_len,
            query_range: 0..query_len,
            states,
            ref_len,
            query_len,
        }
    }

    /// Given the output of an alignment algorithm, generate the aligned
    /// sequences, using `-` as a gap character.
    ///
    /// The arguments `reference` and `query` should be the full sequences
    /// originally passed to the alignment algorithm.
    ///
    /// The first output is the aligned reference, and the second is the aligned
    /// query. Portions of the reference that were not matched are not included
    /// in the output, nor are clipped portions of the query.
    ///
    /// ## Panics
    ///
    /// In addition to the requirements in [Validity and
    /// Panics](Alignment#validity-and-panics), the query and reference must be
    /// the proper length for the alignment.
    #[inline]
    #[must_use]
    pub fn get_aligned_seqs<R: ByteSeq, Q: ByteSeq>(&self, reference: R, query: Q) -> (R::Owned, Q::Owned) {
        let reference = reference.as_ref();
        let query = query.as_ref();

        assert!(
            reference.len() == self.ref_len,
            "Expected reference length {expected}, but found {actual}",
            expected = self.ref_len,
            actual = reference.len()
        );

        assert!(
            query.len() == self.query_len,
            "Expected query length {expected}, but found {actual}",
            expected = self.query_len,
            actual = query.len()
        );

        let (reference, query) = pairwise_align_with(reference, query, &self.states, self.ref_range.start);
        (reference.into(), query.into())
    }

    /// Returns an iterator over the pairs of aligned bases in `reference` and
    /// `query`. `None` is used to represent gaps.
    ///
    /// This has equivalent behavior as [`get_aligned_seqs`] but as an iterator.
    /// [`get_aligned_seqs`] may offer faster performance, but
    /// [`get_aligned_iter`] avoids allocations and may offer a more convenient
    /// syntax for handling gaps.
    ///
    /// ## Panics
    ///
    /// In addition to the requirements in [Validity and
    /// Panics](Alignment#validity-and-panics), the query and reference must be
    /// the proper length for the alignment.
    ///
    /// [`get_aligned_seqs`]: Alignment::get_aligned_seqs
    /// [`get_aligned_iter`]: Alignment::get_aligned_iter
    #[inline]
    #[must_use]
    pub fn get_aligned_iter<R: ByteSeq, Q: ByteSeq>(
        &self, reference: R, query: Q,
    ) -> AlignmentIter<'_, R::IntoIter, Q::IntoIter> {
        let reference_bytes = reference.as_ref();
        let query_bytes = query.as_ref();

        assert!(
            reference_bytes.len() == self.ref_len,
            "Expected reference length {expected}, but found {actual}",
            expected = self.ref_len,
            actual = reference_bytes.len()
        );

        assert!(
            query_bytes.len() == self.query_len,
            "Expected query length {expected}, but found {actual}",
            expected = self.query_len,
            actual = query_bytes.len()
        );

        AlignmentIter::new(
            reference.into_bytes_iter(),
            query.into_bytes_iter(),
            &self.states,
            self.ref_range.start,
        )
    }

    /// Given the output of an alignment algorithm, generate the aligned query,
    /// using `-` as a gap character.
    ///
    /// This is similar to [`get_aligned_seqs`], but only requires passing the
    /// query (and only returns the aligned query). This is useful when the
    /// alignment is against a model or family of references (e.g., a pHMM).
    ///
    /// Portions of the query that were clipped are not included.
    ///
    /// ## Panics
    ///
    /// In addition to the requirements in [Validity and
    /// Panics](Alignment#validity-and-panics), the query must be the proper
    /// length for the alignment.
    ///
    /// [`get_aligned_seqs`]: Alignment::get_aligned_seqs
    pub fn get_aligned_query<Q: ByteSeq>(&self, query: Q) -> Q::Owned {
        let query = query.as_ref();

        assert!(
            query.len() == self.query_len,
            "Expected query length {expected}, but found {actual}",
            expected = self.query_len,
            actual = query.len()
        );

        let mut query_index = 0;
        let mut query_aln = Vec::with_capacity(query.len() + (self.ref_range.len() / 2));

        for Ciglet { inc, op } in &self.states {
            match op {
                AlignOp::M | AlignOp::Eq | AlignOp::X | AlignOp::I => {
                    query_aln.extend_from_slice(&query[query_index..query_index + inc]);
                    query_index += inc;
                }
                AlignOp::D => {
                    query_aln.extend(std::iter::repeat_n(b'-', inc));
                }

                AlignOp::S => query_index += inc,
                AlignOp::N => {
                    query_aln.extend(std::iter::repeat_n(b'N', inc));
                }
                AlignOp::H | AlignOp::P => {}
            }
        }

        query_aln.into()
    }

    /// Returns the number of bases at the start of the query that did not
    /// align.
    ///
    /// If the entire alignment is unmapped, then the query length is split
    /// between [`unaligned_query_start`] and [`unaligned_query_end`] based on
    /// the `query_range` field.
    ///
    /// [`unaligned_query_start`]: Alignment::unaligned_query_start
    /// [`unaligned_query_end`]: Alignment::unaligned_query_end
    pub fn unaligned_query_start(&self) -> usize {
        self.query_range.start
    }

    /// Returns the number of bases at the end of the query that did not align.
    ///
    /// If the entire alignment is unmapped, then the query length is split
    /// between [`unaligned_query_start`] and [`unaligned_query_end`] based on
    /// the `query_range` field.
    ///
    /// [`unaligned_query_start`]: Alignment::unaligned_query_start
    /// [`unaligned_query_end`]: Alignment::unaligned_query_end
    pub fn unaligned_query_end(&self) -> usize {
        self.query_len - self.query_range.end
    }

    /// Returns the number of bases at the start of the reference that did not
    /// align.
    ///
    /// If the entire alignment is unmapped, then the reference length is split
    /// between [`unaligned_ref_start`] and [`unaligned_ref_end`] based on the
    /// `ref_range` field.
    ///
    /// [`unaligned_ref_start`]: Alignment::unaligned_ref_start
    /// [`unaligned_ref_end`]: Alignment::unaligned_ref_end
    pub fn unaligned_ref_start(&self) -> usize {
        self.ref_range.start
    }

    /// Returns the number of bases at the end of the reference that did not
    /// align.
    ///
    /// If the entire alignment is unmapped, then the reference length is split
    /// between [`unaligned_ref_start`] and [`unaligned_ref_end`] based on the
    /// `ref_range` field.
    ///
    /// [`unaligned_ref_start`]: Alignment::unaligned_ref_start
    /// [`unaligned_ref_end`]: Alignment::unaligned_ref_end
    pub fn unaligned_ref_end(&self) -> usize {
        self.ref_len - self.ref_range.end
    }

    /// Returns whether the alignment is unmapped.
    ///
    /// This is based on whether the `ref_range` field is empty.
    pub fn is_unmapped(&self) -> bool {
        self.ref_range.is_empty()
    }
}

impl<T: Copy> Alignment<T> {
    /// Gets the alignment for when the query and reference are swapped.
    ///
    /// ## Validity
    ///
    /// The alignment should not contain any [`AlignOp::N`] operations, since
    /// there is not a corresponding operation for skipped regions from a query.
    #[must_use]
    pub fn invert(&self) -> Self {
        let mut states = AlignmentStates::new();
        states.soft_clip(self.ref_range.start);
        let inverted_ciglets = self.states.0.iter().filter_map(|&ciglet| match ciglet.op {
            AlignOp::S | AlignOp::H => None,
            AlignOp::D => Some(Ciglet {
                inc: ciglet.inc,
                op:  AlignOp::I,
            }),
            AlignOp::I => Some(Ciglet {
                inc: ciglet.inc,
                op:  AlignOp::D,
            }),
            _ => Some(ciglet),
        });
        // Validity: No adjacent ciglets have the same operation since this is
        // true of AlignmentStates and the filter_map call is injective. The
        // increments are non-zero since the same is true of AlignmentStates.
        states.extend_from_ciglets(inverted_ciglets);
        states.soft_clip(self.ref_len - self.ref_range.end);

        Self {
            score: self.score,
            ref_range: self.query_range.clone(),
            query_range: self.ref_range.clone(),
            states,
            ref_len: self.query_len,
            query_len: self.ref_len,
        }
    }

    /// Gets the alignment for when the query and reference are both reversed.
    ///
    /// See [`Reversibility`](Alignment#reversibility) for more details.
    #[must_use]
    pub fn to_reverse(&self) -> Self {
        let ref_range = (self.ref_len - self.ref_range.end)..(self.ref_len - self.ref_range.start);
        let query_range = (self.query_len - self.query_range.end)..(self.query_len - self.query_range.start);
        let states = self.states.to_reverse();

        Alignment {
            score: self.score,
            ref_range,
            query_range,
            states,
            ref_len: self.ref_len,
            query_len: self.query_len,
        }
    }

    /// Converts an alignment in-place so that it represents the alignment
    /// between the reversed query and reference.
    ///
    /// See [`Reversibility`](Alignment#reversibility) for more details.
    pub fn make_reverse(&mut self) {
        self.ref_range = (self.ref_len - self.ref_range.end)..(self.ref_len - self.ref_range.start);
        self.query_range = (self.query_len - self.query_range.end)..(self.query_len - self.query_range.start);
        self.states.make_reverse();
    }

    /// Generates a new alignment replacing `M` in `states` with either `=` for
    /// matches and `X` for mismatches.
    ///
    /// The first argument should contain the full reference that was used to
    /// generate the alignment.
    ///
    /// ## Panics
    ///
    /// If either `reference` or `query` are of a shorter length than implied by
    /// the alignment, then this will panic due to out of bounds indexing.
    #[inline]
    #[must_use]
    pub fn to_verbose_sequence_matching<R: ByteSeq, Q: ByteSeq>(&self, reference: R, query: Q) -> Self {
        Self {
            score:       self.score,
            ref_range:   self.ref_range.clone(),
            query_range: self.query_range.clone(),
            states:      self
                .states
                .to_verbose_sequence_matching(reference.as_ref(), query.as_ref(), self.ref_range.start),
            ref_len:     self.ref_len,
            query_len:   self.query_len,
        }
    }

    /// Clips an [`Alignment`] so that it corresponds to the provided reference
    /// range.
    ///
    /// If the reference range is out of bounds, then `None` is returned. The
    /// `score`, `ref_len`, and `query_len` fields are not changed.
    ///
    /// ## Panics
    ///
    /// This function requires the following assumptions to avoid panicking.
    /// Note that any [`Alignment`] struct returned by a *Zoe* alignment
    /// function will satisfy these.
    ///
    /// - The `query_len` field must be greater than or equal to the number of
    ///   residues that the `states` field consumes in the query.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     alignment::{AlignmentStates, ScalarProfile, ScoreParams, SeqSrc},
    /// #     prelude::*,
    /// # };
    /// // match, mismatch, gap open, gap extend, and ignored base
    /// let params = ScoreParams::new_dna_params(1, -1, -3, -1, None).unwrap();
    ///
    /// let reference = Nucleotides::from(b"GATAATCACATGTGTTGCACGTTGTAAGGTAGCATGCCTTGA");
    /// let query = Nucleotides::from(b"GATATCAAATGCGTGCTTGCACGATGTAAGGTAGC");
    ///
    /// let profile = ScalarProfile::new(&query, &params);
    /// let alignment = profile.sw_align(SeqSrc::Reference(&reference)).unwrap();
    ///
    /// assert_eq!(alignment.states, AlignmentStates::try_from(b"4M1D10M3I18M").unwrap());
    ///
    /// let coding_alignment = alignment.slice_to_ref_range(9..27).unwrap();
    /// assert_eq!(coding_alignment.states, AlignmentStates::try_from(b"8S6M3I12M6S").unwrap());
    /// assert_eq!(coding_alignment.query_range, 8..29);
    /// ```
    pub fn slice_to_ref_range(&self, ref_range: Range<usize>) -> Option<Self> {
        // Adjust the requested ref_range to be an offset from the start of
        // self.ref_range
        let rel_ref_range = if ref_range.start >= self.ref_range.start {
            ref_range.sub(self.ref_range.start)
        } else {
            // We are trying to access a portion of the reference range that was
            // not aligned to
            return None;
        };

        let mut query_idx = 0;
        let mut ref_idx = 0;

        let mut ciglets = self.states.iter().copied();
        let mut states = AlignmentStates::new();

        // Find the first ciglet overlapping the reference range. Strict
        // inequality ensures that any hard or soft clipping ciglets are
        // skipped.
        let Some(mut first_ciglet) = ciglets.find_map(|ciglet| {
            (query_idx, ref_idx) = (query_idx, ref_idx).increment_idxs_by(ciglet);
            let num_in_range = ref_idx.saturating_sub(rel_ref_range.start);
            (num_in_range > 0).then_some(Ciglet {
                op:  ciglet.op,
                inc: num_in_range,
            })
        }) else {
            // If the requested reference range is empty AND is using a number
            // in the appropriate range, then this function will return an empty
            // alignment instead of `None`. The specific range was chosen to
            // mimic `get`.
            if ref_range.is_empty() && (self.ref_range.start..=self.ref_range.end).contains(&ref_range.start) {
                let mut states = AlignmentStates::with_capacity(1);
                states.soft_clip(self.query_len);
                return Some(Alignment {
                    score: self.score,
                    ref_range: ref_range.clone(),
                    query_range: 0..0,
                    states,
                    ref_len: self.ref_len,
                    query_len: self.query_len,
                });
            }

            if !self.ref_range.is_empty() {
                debug_assert!(ref_range.start >= self.ref_range.end);
            }

            return None;
        };

        // Calculates where the sliced alignment starts in the query
        let (query_start, ref_start) = (query_idx, ref_idx).decrement_idxs_by(first_ciglet);
        debug_assert_eq!(ref_start, rel_ref_range.start);

        // Add soft clipping at start
        states.soft_clip(query_start);

        // Check whether the sliced alignment is only a single ciglet
        if ref_idx >= rel_ref_range.end {
            first_ciglet.inc = first_ciglet.inc.saturating_sub(ref_idx - rel_ref_range.end);
            states.add_ciglet(first_ciglet);
            let (query_end, ref_end) = (query_start, ref_start).increment_idxs_by(first_ciglet);
            states.soft_clip(self.query_len - query_end);

            debug_assert_eq!(ref_end, rel_ref_range.end);

            return Some(Self {
                score: self.score,
                ref_range,
                query_range: query_start..query_end,
                states,
                ref_len: self.ref_len,
                query_len: self.query_len,
            });
        }

        // Add the first ciglet
        states.add_ciglet(first_ciglet);

        for mut ciglet in ciglets {
            // Add the contribution of the next ciglet, but don't update idxs
            // yet
            let (new_query_idx, new_ref_idx) = (query_idx, ref_idx).increment_idxs_by(ciglet);

            // Check whether ciglet caused us to meet or pass the end of
            // rel_ref_range
            if let Some(num_past) = new_ref_idx.checked_sub(rel_ref_range.end) {
                ciglet.inc = ciglet.inc.saturating_sub(num_past);
                states.add_ciglet(ciglet);
                let (query_end, ref_end) = (query_idx, ref_idx).increment_idxs_by(ciglet);
                states.soft_clip(self.query_len - query_end);

                debug_assert!(ciglet.inc > 0);
                debug_assert_eq!(ref_end, rel_ref_range.end);

                return Some(Self {
                    score: self.score,
                    ref_range,
                    query_range: query_start..query_end,
                    states,
                    ref_len: self.ref_len,
                    query_len: self.query_len,
                });
            }

            // Add ciglet and update idxs
            states.add_ciglet(ciglet);
            (query_idx, ref_idx) = (new_query_idx, new_ref_idx);
        }

        if ref_idx < rel_ref_range.end {
            debug_assert!(ref_range.end > self.ref_range.end);
            return None;
        }

        debug_assert_eq!(ref_idx, rel_ref_range.end);
        let query_range = query_start..query_idx;
        states.soft_clip(self.query_len - query_idx);
        Some(Self {
            score: self.score,
            ref_range,
            query_range,
            states,
            ref_len: self.ref_len,
            query_len: self.query_len,
        })
    }
}

/// An extension trait for adding or subtracting a [`Ciglet`] to a tuple
/// representing the query index and reference index.
pub(crate) trait AlignmentIndices {
    /// Adds the contribution of the [`Ciglet`] to a query and reference index.
    fn increment_idxs_by(self, ciglet: Ciglet) -> Self;

    /// Subtracts the contribution of the [`Ciglet`] to a query and reference
    /// index.
    fn decrement_idxs_by(self, ciglet: Ciglet) -> Self;
}

impl AlignmentIndices for (usize, usize) {
    #[inline]
    fn increment_idxs_by(self, ciglet: Ciglet) -> Self {
        let (query_idx, ref_idx) = self;
        match ciglet.op {
            AlignOp::M | AlignOp::Eq | AlignOp::X => (query_idx + ciglet.inc, ref_idx + ciglet.inc),
            AlignOp::D | AlignOp::N => (query_idx, ref_idx + ciglet.inc),
            AlignOp::I | AlignOp::S => (query_idx + ciglet.inc, ref_idx),
            AlignOp::H | AlignOp::P => (query_idx, ref_idx),
        }
    }

    #[inline]
    fn decrement_idxs_by(self, ciglet: Ciglet) -> Self {
        let (query_idx, ref_idx) = self;
        match ciglet.op {
            AlignOp::M | AlignOp::Eq | AlignOp::X => (query_idx - ciglet.inc, ref_idx - ciglet.inc),
            AlignOp::D | AlignOp::N => (query_idx, ref_idx - ciglet.inc),
            AlignOp::I | AlignOp::S => (query_idx - ciglet.inc, ref_idx),
            AlignOp::H | AlignOp::P => (query_idx, ref_idx),
        }
    }
}
