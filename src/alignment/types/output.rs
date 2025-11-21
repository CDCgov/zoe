use crate::{
    alignment::{AlignmentIter, AlignmentStates, pairwise_align_with},
    data::{cigar::Ciglet, views::IndexAdjustable},
};
use std::ops::Range;

/// The output of an alignment algorithm, allowing unmapped and overflowed
/// alignments to be represented.
///
/// Often this is some [`Alignment`], but if no portion of the query mapped to
/// the reference, [`Unmapped`] is used. [`Overflowed`] is used when the score
/// exceeded numeric type range.
///
/// [`Some`]: MaybeAligned::Some
/// [`Overflowed`]: MaybeAligned::Overflowed
/// [`Unmapped`]: MaybeAligned::Unmapped
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum MaybeAligned<T> {
    /// Contains a successful alignment result
    Some(T),
    /// Indicates the alignment score overflowed the numeric type
    Overflowed,
    /// Indicates the sequence could not be mapped/aligned
    Unmapped,
}

impl<T> MaybeAligned<T> {
    /// Unwraps the alignment, consuming the [`MaybeAligned<T>`] and returning
    /// the contained [`Alignment<T>`].
    ///
    /// # Panics
    ///
    /// Panics if the value is [`MaybeAligned::Overflowed`] or
    /// [`MaybeAligned::Unmapped`].
    #[inline]
    #[must_use]
    pub fn unwrap(self) -> T {
        match self {
            MaybeAligned::Some(aln) => aln,
            MaybeAligned::Overflowed => panic!("Alignment score overflowed!"),
            MaybeAligned::Unmapped => panic!("Sequence could not be mapped!"),
        }
    }

    /// Gets the contained [`Alignment`] as an [`Option`].
    #[inline]
    #[must_use]
    pub fn get(self) -> Option<T> {
        match self {
            MaybeAligned::Some(aln) => Some(aln),
            _ => None,
        }
    }

    /// Returns the alignment if it did not overflow, otherwise calls `f` and
    /// returns the result.
    #[inline]
    #[must_use]
    pub fn or_else_overflowed(self, f: impl FnOnce() -> Self) -> Self {
        if let MaybeAligned::Overflowed = self { f() } else { self }
    }

    /// Maps a [`MaybeAligned<T>`] to [`MaybeAligned<U>`] by applying a function
    /// to the contained value. Leaves [`Overflowed`] and [`Unmapped`] variants
    /// unchanged.
    ///
    /// [`Overflowed`]: MaybeAligned::Overflowed
    /// [`Unmapped`]: MaybeAligned::Unmapped
    #[inline]
    #[must_use]
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> MaybeAligned<U> {
        match self {
            MaybeAligned::Some(value) => MaybeAligned::Some(f(value)),
            MaybeAligned::Overflowed => MaybeAligned::Overflowed,
            MaybeAligned::Unmapped => MaybeAligned::Unmapped,
        }
    }

    /// Returns [`Unmapped`] if the [`MaybeAligned`] is [`Unmapped`],
    /// [`Overflowed`] if it is [`Overflowed`], otherwise calls `f` with the
    /// wrapped value and returns the result.
    ///
    /// [`Overflowed`]: MaybeAligned::Overflowed
    /// [`Unmapped`]: MaybeAligned::Unmapped
    #[inline]
    #[must_use]
    pub fn and_then<U, F>(self, f: F) -> MaybeAligned<U>
    where
        F: FnOnce(T) -> MaybeAligned<U>, {
        match self {
            MaybeAligned::Some(x) => f(x),
            MaybeAligned::Overflowed => MaybeAligned::Overflowed,
            MaybeAligned::Unmapped => MaybeAligned::Unmapped,
        }
    }
}

// For the `Alignment` struct below, both ranges are 0-based and end-exclusive.
// For a global alignment, the ranges will each be encompass the full length of
// the sequences. For local alignment, `query_range` does NOT include hard or
// soft clipped bases, even though `states` does contain this information.

/// A struct representing the information for an alignment, such as its score
/// and where in the sequences it occurs.
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
    /// `reference` and `query` should be the full sequences originally passed
    /// to the alignment algorithm.
    ///
    /// The first output is the aligned reference, and the second is the aligned
    /// query. Portions of the reference that were not matched are not included
    /// in the output, nor are clipped portions of the query.
    ///
    /// ## Panics
    ///
    /// This function requires the following assumptions to avoid panicking.
    /// Note that any [`Alignment`] struct returned by a *Zoe* alignment
    /// function will satisfy these.
    ///
    /// - All operations must be valid state, e.g, bytes in `MIDNSHP=X`.
    /// - The query and reference must be at least as long as the length implied
    ///   by the alignment operations.
    #[inline]
    #[must_use]
    pub fn get_aligned_seqs(&self, reference: &[u8], query: &[u8]) -> (Vec<u8>, Vec<u8>) {
        pairwise_align_with(reference, query, &self.states, self.ref_range.start)
    }

    /// Returns an iterator over the pairs of aligned bases in `reference` and
    /// `query`. `None` is used to represent gaps.
    ///
    /// This has equivalent behavior as [`get_aligned_seqs`] but as an iterator.
    /// [`get_aligned_seqs`] may offer faster performance, but
    /// [`get_aligned_iter`] avoids allocations and may offer a more convenient
    /// syntax for handling gaps.
    ///
    /// The same conditions are required as [`get_aligned_seqs`] in order for
    /// the iterator to be non-panicking.
    ///
    /// [`get_aligned_seqs`]: Alignment::get_aligned_seqs
    /// [`get_aligned_iter`]: Alignment::get_aligned_iter
    #[inline]
    #[must_use]
    pub fn get_aligned_iter<'a>(
        &'a self, reference: &'a [u8], query: &'a [u8],
    ) -> AlignmentIter<'a, std::iter::Copied<std::slice::Iter<'a, Ciglet>>> {
        AlignmentIter::new(reference, query, &self.states, self.ref_range.start)
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
    /// This function requires the following assumptions to avoid panicking.
    /// Note that any [`Alignment`] struct returned by a *Zoe* alignment
    /// function will satisfy these.
    ///
    /// - All operations must be valid state, e.g, bytes in `MIDNSHP=X`.
    /// - The query must be at least as long as the length implied by the
    ///   alignment operations.
    ///
    /// [`get_aligned_seqs`]: Alignment::get_aligned_seqs
    #[allow(clippy::missing_panics_doc)]
    pub fn get_aligned_query(&self, query: &[u8]) -> Vec<u8> {
        let mut query_index = 0;
        let mut query_aln = Vec::with_capacity(query.len() + (self.ref_range.len() / 2));

        for Ciglet { inc, op } in &self.states {
            match op {
                b'M' | b'=' | b'X' | b'I' => {
                    query_aln.extend_from_slice(&query[query_index..query_index + inc]);
                    query_index += inc;
                }
                b'D' => {
                    query_aln.extend(std::iter::repeat_n(b'-', inc));
                }

                b'S' => query_index += inc,
                b'N' => {
                    query_aln.extend(std::iter::repeat_n(b'N', inc));
                }
                b'H' | b'P' => {}
                // TODO: Could ignore if we allow only valid CIGAR by default.
                _ => panic!("CIGAR op '{op}' not supported.\n"),
            }
        }

        query_aln
    }
}

impl<T: Copy> Alignment<T> {
    /// Gets the alignment for when the query and reference are swapped.
    #[must_use]
    pub fn invert(&self) -> Self {
        let mut states = AlignmentStates::new();
        states.soft_clip(self.ref_range.start);
        let inverted_ciglets = self.states.0.iter().filter_map(|&ciglet| match ciglet.op {
            b'S' | b'H' => None,
            b'D' => Some(Ciglet {
                inc: ciglet.inc,
                op:  b'I',
            }),
            b'I' => Some(Ciglet {
                inc: ciglet.inc,
                op:  b'D',
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

    /// Gets the alignment for when the query and reference are both reversed
    /// (or reverse-complemented).
    #[must_use]
    pub fn to_reverse(&self) -> Self {
        let ref_range = (self.ref_len - self.ref_range.end)..(self.ref_len - self.ref_range.start - 1);
        let query_range = (self.query_len - self.query_range.end)..(self.query_len - self.query_range.start - 1);
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
    /// between the reversed (or reverse-complemented) query and reference.
    pub fn make_reverse(&mut self) {
        self.ref_range = (self.ref_len - self.ref_range.end)..(self.ref_len - self.ref_range.start - 1);
        self.query_range = (self.query_len - self.query_range.end)..(self.query_len - self.query_range.start - 1);
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
    pub fn to_verbose_sequence_matching(&self, reference: &[u8], query: &[u8]) -> Self {
        Self {
            score:       self.score,
            ref_range:   self.ref_range.clone(),
            query_range: self.query_range.clone(),
            states:      self
                .states
                .to_verbose_sequence_matching(reference, query, self.ref_range.start),
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
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     alignment::{AlignmentStates, ScalarProfile},
    /// #     data::WeightMatrix,
    /// #     prelude::*,
    /// # };
    /// const MATRIX: WeightMatrix<'_, i8, 5> = WeightMatrix::new_dna_matrix(1, -1, None);
    /// const GAP_OPEN: i8 = -3;
    /// const GAP_EXTEND: i8 = -1;
    ///
    /// let reference = Nucleotides::from(b"GATAATCACATGTGTTGCACGTTGTAAGGTAGCATGCCTTGA");
    /// let query = Nucleotides::from(b"GATATCAAATGCGTGCTTGCACGATGTAAGGTAGC");
    ///
    /// let profile = ScalarProfile::new(&query, &MATRIX, GAP_OPEN, GAP_EXTEND).unwrap();
    /// let alignment = profile.smith_waterman_alignment(reference.as_bytes()).unwrap();
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

        // Find the first ciglet overlapping the reference range
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
            b'M' | b'=' | b'X' => (query_idx + ciglet.inc, ref_idx + ciglet.inc),
            b'D' | b'N' => (query_idx, ref_idx + ciglet.inc),
            b'I' | b'S' => (query_idx + ciglet.inc, ref_idx),
            _ => (query_idx, ref_idx),
        }
    }

    #[inline]
    fn decrement_idxs_by(self, ciglet: Ciglet) -> Self {
        let (query_idx, ref_idx) = self;
        match ciglet.op {
            b'M' | b'=' | b'X' => (query_idx - ciglet.inc, ref_idx - ciglet.inc),
            b'D' | b'N' => (query_idx, ref_idx - ciglet.inc),
            b'I' | b'S' => (query_idx - ciglet.inc, ref_idx),
            _ => (query_idx, ref_idx),
        }
    }
}
