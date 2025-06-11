use crate::{
    data::{
        byte_types::IsBase,
        cigar::{Ciglet, ExpandedCigar},
        sam::SamData,
    },
    prelude::QualityScores,
    search::ByteSubstring,
};
use std::{
    iter::repeat_n,
    ops::{Add, AddAssign, Range},
};

/// Represents a single insertion in a [`SamAligned`] record.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
struct SamInsertion {
    /// The reference index (0-based) before which the insertion occurs.
    ref_index:   usize,
    /// The starting index (0-based) of the insertion in the query sequence.
    query_start: usize,
    /// The ending index (0-based, non-inclusive) of the insertion in the query
    ///   sequence.
    query_end:   usize,
}

impl SamInsertion {
    #[inline]
    #[must_use]
    pub fn query_range(&self) -> Range<usize> {
        self.query_start..self.query_end
    }
}

/// Struct holding the alignment information of a SAM query sequence to some
/// reference.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
struct SamAligned {
    /// Query bases aligned to the reference (insertions removed, deletions
    /// added).
    aligned:           Vec<u8>,
    /// Query quality scores aligned to the reference (insertions removed,
    /// deletions added).
    qaligned:          Vec<u8>,
    /// Start of range of reference indices. Start is inclusive.
    ref_start:         usize,
    /// End of range of reference indices. The `ref_end` is *exclusive*, thus
    /// the final aligned query base will be at `ref_end - 1`.
    ref_end:           usize,
    /// Vector of query insertion relative to reference. Ordered by reference
    /// indices.
    insertions:        Vec<SamInsertion>,
    /// Number of bases hard/soft clipped at the start
    num_clipped_start: usize,
    /// Number of bases hard/soft clipped at the end
    num_clipped_end:   usize,
}

impl SamAligned {
    /// Initializes a new [`SamAligned`] object from the fields.
    #[inline]
    fn new(
        aligned: Vec<u8>, qaligned: Vec<u8>, ref_start: usize, ref_end: usize, insertions: Vec<SamInsertion>,
        num_clipped_start: usize, num_clipped_end: usize,
    ) -> Self {
        SamAligned {
            aligned,
            qaligned,
            ref_start,
            ref_end,
            insertions,
            num_clipped_start,
            num_clipped_end,
        }
    }

    /// At the given reference index, provides the reference-aligned query's
    /// nucleotide and encoded (ASCII) quality score as an optional tuple of
    /// `(base, qs)`, otherwise, `None` is returned.
    #[must_use]
    fn get_base_and_quality(&self, at_index: usize) -> Option<(u8, u8)> {
        // Safety: reference range is the same length as `aligned` and
        // `qaligned`; thus, subtracting the start of the range will render
        // appropriate 0-based index. If the reference index it outside the
        // range, then `None` is returned.
        if self.ref_range().contains(&at_index) {
            Some((
                self.aligned[at_index - self.ref_start],
                self.qaligned[at_index - self.ref_start],
            ))
        } else {
            None
        }
    }

    /// At the given reference index, provides the reference-aligned query's
    /// nucleotide as an `Option`, otherwise, `None` is returned.
    #[must_use]
    fn get_base(&self, at_index: usize) -> Option<u8> {
        // Safety: reference range is the same length as `aligned` and
        // `qaligned`; thus, subtracting the start of the range will render
        // appropriate 0-based index. If the reference index it outside the
        // range, then `None` is returned.
        if self.ref_range().contains(&at_index) {
            Some(self.aligned[at_index - self.ref_start])
        } else {
            None
        }
    }

    /// Return a half-open reference range as `Range<uize>` for the aligned query.
    #[must_use]
    fn ref_range(&self) -> Range<usize> {
        self.ref_start..self.ref_end
    }

    /// Merges two aligned queries' reference ranges such that the combined
    /// range spans both aligned regions.
    #[must_use]
    fn merge_ref_range(&self, other: &SamAligned) -> Range<usize> {
        std::cmp::min(self.ref_start, other.ref_start)..std::cmp::max(self.ref_end, other.ref_end)
    }

    /// Checks if the `SamAligned` contains an insertion after the 0-based
    /// reference index and returns the query's unaligned range.
    #[must_use]
    fn get_insert_after(&self, reference_index: usize) -> Option<SamInsertion> {
        self.insertions
            .iter()
            .find(|insert| insert.ref_index == reference_index)
            .cloned()
    }
}

impl From<&SamData> for SamAligned {
    #[inline]
    fn from(row: &SamData) -> Self {
        row.get_aligned()
    }
}

/// Makes a merged qname for paired-end reads based on a provided qname plus
/// adding or replacing `3` as the read side. Used in IRMA.
pub(crate) fn make_merged_qname(s: &str) -> String {
    let mut merged = String::with_capacity(s.len() + 2);
    if let Some(index) = s.find(' ') {
        let (head, tail) = s.split_at(index);

        if !(head.starts_with("SRR") || head.starts_with("DRR") || head.starts_with("ERR")) || !head.contains('.') {
            // Illumina format
            merged.push_str(head);
            merged.push_str(" 3");
            if let Some(index) = tail.find(':') {
                merged.push_str(&tail[index..]);
            }
        } else if let Some(index) = head.match_indices('.').nth(1).map(|(i, _)| i) {
            // SRA format, read side included
            let (head, _) = head.split_at(index);
            merged.push_str(head);
            merged.push_str(".3");
            merged.push_str(tail);
        } else {
            // SRA format, no read side
            merged.push_str(head);
            merged.push_str(".3");
            merged.push_str(tail);
        }
    } else if let Some(index) = s.find('/') {
        // Legacy Illumina
        let (id, _) = s.split_at(index);
        merged.push_str(id);
        merged.push_str("/3");
    } else if (s.starts_with("SRR") || s.starts_with("DRR") || s.starts_with("ERR")) && s.contains('.') {
        let mut pieces = s.split('_');
        let id = pieces.next().unwrap_or_default();

        if let Some(index) = id.match_indices('.').nth(1).map(|(i, _)| i) {
            // SRA with read side
            let (new_id, _) = id.split_at(index);
            merged.push_str(new_id);
            merged.push_str(".3");
            for piece in pieces {
                merged.push('_');
                merged.push_str(piece);
            }
        } else {
            // SRA, no read side
            merged.push_str(id);
            merged.push_str(".3");
            for piece in pieces {
                merged.push('_');
                merged.push_str(piece);
            }
        }
    } else {
        // IRMA Illumina legacy output
        let mut indices = s.match_indices(':');
        let (left, right) = (indices.nth(5), indices.next());
        if let (Some((start, _)), Some((stop, _))) = (left, right)
            && let Some(us) = s[start..stop].find('_')
        {
            let underscore_index = start + us;
            merged.push_str(&s[..underscore_index]);
            merged.push_str("_3");
            merged.push_str(&s[stop..]);
        }
    }
    merged
}

/// [`PairedMergeStats`] holds statistics related to read pair merging operations.
#[derive(Copy, Clone, Debug, Default)]
pub struct PairedMergeStats {
    /// Total number of overlapping or paired bases
    pub observations:    u64,
    /// Paired bases that agree with each other but disagree with consensus
    pub true_variations: u64,
    /// Paired bases that disagree with each other
    pub variant_errors:  u64,
    /// Paired bases where one is a deletion and one is not
    pub deletion_errors: u64,
    /// Total number of paired insertions, in agreement or otherwise
    pub insert_obs:      u64,
    /// Total number of mismatching paired insertions, including disagreement in
    /// insertion presence
    pub insert_errors:   u64,
}

impl Add for PairedMergeStats {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self {
            observations:    self.observations + other.observations,
            true_variations: self.true_variations + other.true_variations,
            variant_errors:  self.variant_errors + other.variant_errors,
            deletion_errors: self.deletion_errors + other.deletion_errors,
            insert_obs:      self.insert_obs + other.insert_obs,
            insert_errors:   self.insert_errors + other.insert_errors,
        }
    }
}

impl AddAssign for PairedMergeStats {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            observations:    self.observations + other.observations,
            true_variations: self.true_variations + other.true_variations,
            variant_errors:  self.variant_errors + other.variant_errors,
            deletion_errors: self.deletion_errors + other.deletion_errors,
            insert_obs:      self.insert_obs + other.insert_obs,
            insert_errors:   self.insert_errors + other.insert_errors,
        }
    }
}

impl SamData {
    /// Provides a struct [`SamAligned`], which contains the aligned sequence
    /// and quality scores based on the [`Cigar`] as well as the alignment range
    /// for the reference.
    ///
    /// ## Panics
    ///
    /// Currently panics on invalid cigar states. To be removed in the future
    /// for either a `Result` or type-state validated Cigars.
    #[must_use]
    fn get_aligned(&self) -> SamAligned {
        // Convert SAM 1-based to 0-based
        let mut ref_index = self.pos - 1;
        let mut query_index = 0;

        let mut aln: Vec<u8> = Vec::new();
        let mut q_aln: Vec<u8> = Vec::new();
        let mut insertions = Vec::new();

        let mut ciglets = self.cigar.iter();

        let num_clipped_start = ciglets.remove_clipping();
        let num_clipped_end = ciglets.remove_clipping_back();

        for Ciglet { inc, op } in &self.cigar {
            match op {
                b'M' | b'=' | b'X' => {
                    for _ in 0..inc {
                        q_aln.push(self.qual[query_index]);
                        aln.push(self.seq[query_index]);
                        query_index += 1;
                        ref_index += 1;
                    }
                }
                b'D' => {
                    for _ in 0..inc {
                        // We use the minimum value for deletions so that it can
                        // still be parsed by other programs.
                        q_aln.push(b'!');
                        aln.push(b'-');
                    }
                    ref_index += inc;
                }
                b'I' => {
                    let ins = SamInsertion {
                        // Insertion reference index uses 5' site (0-based)
                        // by convention.
                        ref_index:   ref_index - 1,
                        query_start: query_index,
                        query_end:   query_index + inc,
                    };
                    insertions.push(ins);
                    query_index += inc;
                }
                b'S' => query_index += inc,
                b'N' => {
                    for _ in 0..inc {
                        q_aln.push(b'!');
                        aln.push(b'N');
                    }
                    ref_index += inc;
                }
                b'H' => {}
                // TODO: Replace with a continue if type-state is adopted.
                _ => panic!("Extended CIGAR {op} not yet supported.\n"),
            }
        }

        // Convert start position from 1-based to 0-based, then use 0-based ref_index
        // with half-open range: (self.pos - 1)..ref_index;

        SamAligned::new(
            aln,
            q_aln,
            self.pos - 1,
            ref_index,
            insertions,
            num_clipped_start,
            num_clipped_end,
        )
    }

    /// Merges SAM read pairs using the reference alignment to parsimoniously
    /// detect and correct errors. Based on the work by Shepard et al. 2016 for
    /// IRMA.
    #[allow(clippy::too_many_lines)]
    #[must_use]
    pub fn merge_pair_using_reference(
        &self, other: &SamData, reference: &[u8], bowtie_format: bool,
    ) -> (SamData, PairedMergeStats) {
        let mut stats = PairedMergeStats::default();

        let m_qname = if bowtie_format {
            self.qname.clone()
        } else {
            // IRMA merged style: set to 3
            make_merged_qname(&self.qname)
        };

        match (self.is_unmapped(), other.is_unmapped()) {
            (true, true) => return (SamData::unmapped(&m_qname, &self.rname), stats),
            (_, true) => {
                let mut new = self.clone();
                new.qname = m_qname;
                return (new, stats);
            }
            (true, _) => {
                let mut new = other.clone();
                new.qname = m_qname;
                return (new, stats);
            }
            _ => {}
        }

        let a1 = self.get_aligned();
        let a2 = other.get_aligned();

        let paired_range = a1.merge_ref_range(&a2);

        let m_num_clipped_start = {
            // How far is each read from the start of the merged read?
            let a1_offset_from_merged = a1.ref_start - paired_range.start;
            let a2_offset_from_merged = a2.ref_start - paired_range.start;
            // How many clipped bases extend left of paired_range.start?
            let a1_clipping_for_paired = a1.num_clipped_start.saturating_sub(a1_offset_from_merged);
            let a2_clipping_for_paired = a2.num_clipped_start.saturating_sub(a2_offset_from_merged);
            // Take the maximum of these (whichever has more clipped bases)
            a1_clipping_for_paired.max(a2_clipping_for_paired)
        };

        let m_num_clipped_end = {
            // How far is each read from the end of the merged read?
            let a1_offset_from_merged = paired_range.end - a1.ref_end;
            let a2_offset_from_merged = paired_range.end - a2.ref_end;
            // How many clipped bases extend right of paired_range.end?
            let a1_clipping_for_paired = a1.num_clipped_end.saturating_sub(a1_offset_from_merged);
            let a2_clipping_for_paired = a2.num_clipped_end.saturating_sub(a2_offset_from_merged);
            // Take the maximum of these (whichever has more clipped bases)
            a1_clipping_for_paired.max(a2_clipping_for_paired)
        };

        let m_flag = 0;
        let m_rname = self.rname.clone();

        let (mapq1, mapq2) = (self.mapq, other.mapq);
        let m_mapq = mapq1.midpoint(mapq2);
        let m_pos = paired_range.start + 1;

        let mut merged_cigars = Vec::with_capacity(paired_range.len());
        let mut merged_seq = Vec::with_capacity(paired_range.len());
        let mut merged_quals = Vec::with_capacity(paired_range.len());

        merged_cigars.extend(std::iter::repeat_n(b'H', m_num_clipped_start));

        // 0-based index relative to reference
        for ref_index in paired_range {
            let p1 = a1.get_base_and_quality(ref_index);
            let p2 = a2.get_base_and_quality(ref_index);
            let r = reference[ref_index];

            match (p1, p2) {
                (Some((x, qx)), Some((y, qy))) => {
                    stats.observations += 1;
                    if x == y {
                        if x == b'-' {
                            stats.true_variations += 1;
                            merged_cigars.push(b'D');
                        } else {
                            if x != r {
                                stats.true_variations += 1;
                            }
                            merged_cigars.push(b'M');

                            merged_seq.push(x);
                            merged_quals.push(std::cmp::max(qx, qy));
                        }
                    } else {
                        // x ≠ y
                        stats.variant_errors += 1;
                        merged_cigars.push(b'M');

                        if x == r {
                            if y == b'-' {
                                stats.deletion_errors += 1;
                            }

                            merged_seq.push(x);
                            merged_quals.push(qx);
                        } else if y == r {
                            if x == b'-' {
                                stats.deletion_errors += 1;
                            }

                            merged_seq.push(y);
                            merged_quals.push(qy);
                        // x ≠ y
                        } else if y == b'-' {
                            stats.deletion_errors += 1;

                            merged_seq.push(x);
                            merged_quals.push(qx);
                        } else if x == b'-' {
                            stats.deletion_errors += 1;

                            merged_seq.push(y);
                            merged_quals.push(qy);

                        // x, y ≠ r, x ≠ y
                        } else if x.is_base_acgt() && y.is_not_base_acgt() {
                            merged_seq.push(x);
                            merged_quals.push(qx);
                        } else if y.is_base_acgt() && x.is_not_base_acgt() {
                            merged_seq.push(y);
                            merged_quals.push(qy);

                        // Saturating isn't strictly necessary since the
                        // upper range would be 252 = (93 + 33) * 2
                        //
                        // NB: x, y must be bases
                        } else if qx > qy.saturating_add(4) {
                            merged_seq.push(x);
                            merged_quals.push(qx);
                        } else if qy > qx.saturating_add(4) {
                            merged_seq.push(y);
                            merged_quals.push(qy);
                        } else {
                            merged_seq.push(b'N');
                            merged_quals.push(qx.midpoint(qy));
                        }
                    }
                }
                (Some((base, quality)), None) | (None, Some((base, quality))) => {
                    if base == b'-' {
                        merged_cigars.push(b'D');
                    } else {
                        merged_cigars.push(b'M');
                        merged_seq.push(base);
                        merged_quals.push(quality); // ensure this is an encoded value during testing
                    }
                }
                (None, None) => merged_cigars.push(b'N'),
            }

            let range1 = a1.get_insert_after(ref_index);
            let range2 = a2.get_insert_after(ref_index);

            match (range1, range2) {
                (Some(r1), Some(r2)) => {
                    stats.insert_obs += 1;
                    let insert1 = &self.seq[r1.query_range()];
                    let quals1 = &self.qual[r1.query_range()];

                    let insert2 = &other.seq[r2.query_range()];
                    let quals2 = &other.qual[r2.query_range()];

                    if insert1 == insert2 {
                        merged_seq.extend_from_slice(insert1.to_ascii_lowercase().as_slice());
                        merged_quals.extend(quals1.iter().zip(quals2).map(|(q1, q2)| std::cmp::max(*q1, *q2)));
                        merged_cigars.extend(repeat_n(b'I', insert1.len()));
                    } else if insert2.contains_substring(insert1) {
                        merged_seq.extend_from_slice(insert1.to_ascii_lowercase().as_slice());
                        merged_quals.extend_from_slice(quals1);
                        merged_cigars.extend(repeat_n(b'I', insert1.len()));

                        stats.insert_errors += 1;
                    } else if insert1.contains_substring(insert2) {
                        merged_seq.extend_from_slice(insert2.to_ascii_lowercase().as_slice());
                        merged_quals.extend_from_slice(quals2);
                        merged_cigars.extend(repeat_n(b'I', insert2.len()));

                        stats.insert_errors += 1;
                    } else {
                        stats.insert_errors += 1;
                    }
                }
                (Some(r1), None) => {
                    // Was an insertion possible along seq2?
                    let pp2 = a2.get_base(ref_index + 1);

                    // If so, then that's a discrepancy we will ignore
                    if p2.is_some() && pp2.is_some() {
                        stats.insert_obs += 1;
                        stats.insert_errors += 1;
                    // If unmated at this locus, then just add the insert
                    } else {
                        let insert = &self.seq[r1.query_range()];
                        let quals = &self.qual[r1.query_range()];

                        merged_seq.extend_from_slice(insert.to_ascii_lowercase().as_slice());
                        merged_quals.extend_from_slice(quals);
                        merged_cigars.extend(repeat_n(b'I', insert.len()));
                    }
                }
                (None, Some(r2)) => {
                    // Was an insertion possible along seq1?
                    let pp1 = a1.get_base(ref_index + 1);

                    // If so, then that's a discrepancy we will ignore
                    if p1.is_some() && pp1.is_some() {
                        stats.insert_obs += 1;
                        stats.insert_errors += 1;
                    } else {
                        let insert = &other.seq[r2.query_range()].to_ascii_lowercase();
                        let quals = &other.qual[r2.query_range()];

                        // Remove lower-casing if it does nothing downstream.
                        merged_seq.extend_from_slice(insert.to_ascii_lowercase().as_slice());
                        merged_quals.extend_from_slice(quals);
                        merged_cigars.extend(repeat_n(b'I', insert.len()));
                    }
                }
                _ => {}
            }
        }

        merged_cigars.extend(std::iter::repeat_n(b'H', m_num_clipped_end));

        let merged_cigars = ExpandedCigar::from(merged_cigars).condense_to_cigar();

        (
            SamData::new(
                m_qname,
                m_flag,
                m_rname,
                m_pos,
                m_mapq,
                merged_cigars,
                merged_seq.into(),
                // Safety: QualityScores guarantee being in range. The merge quality
                // scores must be from either sequence OR the integer average, which
                // must still be valid state. Even deletions, which have no quality
                // score, is initialized to the minimum QS encoded value `!`.
                unsafe { QualityScores::from_vec_unchecked(merged_quals) },
            ),
            stats,
        )
    }
}
