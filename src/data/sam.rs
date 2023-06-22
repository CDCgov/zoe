use crate::data::byte_types::{ByteAvg, IsBase};
use crate::data::types::{
    cigar::{Cigar, ExpandedCigar},
    nucleotides::Nucleotides,
    phred::QualityScores,
};
use crate::data::ContainsSubsequence;
use lazy_static::lazy_static;
use regex::Regex;
use std::io::{BufRead, Error as IOError, ErrorKind};
use std::iter::repeat;
use std::ops::{Add, AddAssign, Range};
use std::{fs::File, path::Path};

// We define `index` to be 1-based and `position` to be 0-based to avoid off-by-one errors and encourage better semantics

pub enum SamRow {
    Header(String),
    Data(Box<SamData>),
}

#[derive(Debug, Clone)]
pub struct SamData {
    /// Query name.
    pub qname: String,
    /// SAM flag: strandedness, etc.
    pub flag:  u16,
    /// Reference name.
    pub rname: String,
    /// Query's 1-based reference position after alignment.
    pub pos:   usize,
    /// Mystical map quality value.
    pub mapq:  u8,
    /// Old style cigar format that does not include match and mismatch as separate values.
    pub cigar: Cigar,
    /// Reference name of the mate / next read. Currently not implemented and set to `*`.
    rnext:     char,
    /// Position of the mate / next read. Currently not implemented and st to `0`.
    pnext:     u32,
    /// So-called "observed template length." Currently not implemented and always set to `0`.
    tlen:      i32,
    /// Query sequence.
    pub seq:   Nucleotides,
    /// Query quality scores in ASCII-encoded format with Phred Quality of +33.
    pub qual:  QualityScores,
}

impl SamData {
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        qname: String, flag: u16, rname: String, pos: usize, mapq: u8, cigar: Cigar, seq: Nucleotides, qual: QualityScores,
    ) -> Self {
        SamData {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext: '*',
            pnext: 0,
            tlen: 0,
            seq,
            qual,
        }
    }

    /// Provides a struct `SamAligned`, which contains the aligned sequence and
    /// quality scores based on the `Cigar` as well as the alignment range for
    /// the reference.
    ///
    /// # Panics
    ///
    /// Currently panics on invalid cigar states. To be removed in the future
    /// for either a `Result` or type-state validated Cigars.
    #[must_use]
    pub fn get_aligned(&self) -> SamAligned {
        // Convert SAM 1-based to 0-based
        let mut ref_index = self.pos - 1;
        let mut query_index = 0;

        let mut aln: Vec<u8> = Vec::new();
        let mut q_aln: Vec<u8> = Vec::new();
        let mut insertions = Vec::new();

        for (inc, op) in self.cigar.into_iter_tuple() {
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
                        // TO-DO: consider what is correct for q_aln deletion states
                        q_aln.push(b' ');
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
                        q_aln.push(b' ');
                        aln.push(b'N');
                    }
                    ref_index += inc;
                }
                b'H' => continue,
                // TO-DO: Replace with a continue if type-state is adopted.
                _ => panic!("Extended CIGAR {op} not yet supported.\n"),
            }
        }

        // Convert start position from 1-based to 0-based, then use 0-based ref_index
        // with half-open range: (self.pos - 1)..ref_index;

        SamAligned::new(aln, q_aln, self.pos - 1, ref_index, insertions)
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

        lazy_static! {
            static ref UPDATE_PAIR: Regex =
                Regex::new(r"(.+?[_ ])[12](:.+)").expect("REGEX for changing 'qname' didn't compile.");
        }

        let a1: SamAligned = self.get_aligned();
        let a2: SamAligned = other.get_aligned();

        let paired_range = a1.merge_ref_range(&a2);

        let m_qname = if bowtie_format {
            self.qname.clone()
        // IRMA merged style: set to 3
        } else {
            UPDATE_PAIR.replace(&self.qname, r"${1}3${2}").into_owned()
        };

        let m_flag = 0;
        let m_rname = self.rname.clone();

        let (mapq1, mapq2) = (self.mapq, other.mapq);
        let m_mapq = mapq1.avg(mapq2);
        let m_pos = paired_range.start + 1;

        let mut merged_cigars = Vec::with_capacity(paired_range.len());
        let mut merged_seq = Vec::with_capacity(paired_range.len());
        let mut merged_quals = Vec::with_capacity(paired_range.len());

        use std::cmp::max;

        // 0-based index relateive to reference
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
                            merged_quals.push(max(qx, qy));
                        }
                    } else {
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
                        } else if x.is_base() && y.is_not_base() {
                            if y == b'-' {
                                stats.deletion_errors += 1;
                            }

                            merged_seq.push(x);
                            merged_quals.push(qx);
                        } else if y.is_base() && x.is_not_base() {
                            if x == b'-' {
                                stats.deletion_errors += 1;
                            }

                            merged_seq.push(y);
                            merged_quals.push(qy);

                        // Saturating isn't strictly necessary since the
                        // upper range would be 252 = (93 + 33) * 2
                        } else if qx > qy.saturating_add(4) {
                            // Only if "N" vs. "-"
                            if y == b'-' {
                                stats.deletion_errors += 1;
                            }

                            merged_seq.push(x);
                            merged_quals.push(qx);
                        } else if qy > qx.saturating_add(4) {
                            if x == b'-' {
                                stats.deletion_errors += 1;
                            }

                            merged_seq.push(y);
                            merged_quals.push(qy);
                        } else {
                            merged_seq.push(b'N');
                            merged_quals.push(qx.avg(qy));
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

                        quals1
                            .iter()
                            .zip(quals2)
                            .map(|(q1, q2)| max(*q1, *q2))
                            .collect_into(&mut merged_quals);

                        repeat(b'I').take(insert1.len()).collect_into(&mut merged_cigars);
                    } else if insert2.contains_subsequence(insert1) {
                        merged_seq.extend_from_slice(insert1.to_ascii_lowercase().as_slice());
                        merged_quals.extend_from_slice(quals1);
                        repeat(b'I').take(insert1.len()).collect_into(&mut merged_cigars);

                        stats.insert_errors += 1;
                    } else if insert1.contains_subsequence(insert2) {
                        merged_seq.extend_from_slice(insert2.to_ascii_lowercase().as_slice());
                        merged_quals.extend_from_slice(quals2);
                        repeat(b'I').take(insert2.len()).collect_into(&mut merged_cigars);

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
                        repeat(b'I').take(insert.len()).collect_into(&mut merged_cigars);
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
                        repeat(b'I').take(insert.len()).collect_into(&mut merged_cigars);
                    }
                }
                _ => {}
            }
        }

        let merged_cigars = ExpandedCigar::from(merged_cigars).condense_cigar();

        (
            SamData::new(
                m_qname,
                m_flag,
                m_rname,
                m_pos,
                m_mapq,
                merged_cigars,
                merged_seq.into(),
                merged_quals.into(),
            ),
            stats,
        )
    }

    /// Gets an optional slice of query nucleotides using 0-based query indices.
    #[inline]
    pub fn get_bases<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.seq.get(index)
    }

    /// Takes an index or range and return an optional slice of the quality scores for the SAM record.
    #[inline]
    pub fn get_qscores<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.qual.get(index)
    }
}

#[derive(Default, Clone, Debug, Copy)]
pub struct PairedMergeStats {
    pub observations:    u64,
    pub true_variations: u64,
    pub variant_errors:  u64,
    pub deletion_errors: u64,
    pub insert_obs:      u64,
    pub insert_errors:   u64,
}

impl Add for PairedMergeStats {
    type Output = Self;

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

#[derive(Debug, Clone)]
pub struct SamInsertion {
    pub ref_index:   usize,
    pub query_start: usize,
    pub query_end:   usize,
}

impl SamInsertion {
    #[must_use]
    pub fn query_range(&self) -> Range<usize> {
        self.query_start..self.query_end
    }
}

/// Struct holding the alignment information of a SAM query sequence to some reference.
#[derive(Debug, Clone)]
pub struct SamAligned {
    /// Query bases aligned to the reference (insertions removed, deletions added).
    pub aligned:    Vec<u8>,
    /// Query quality scores aligned to the reference (insertions removed, deletions added).
    pub qaligned:   Vec<u8>,
    /// Start of range of reference indices. Start is inclusive.
    pub ref_start:  usize,
    /// End of range of reference indices. The `ref_end` is *exclusive*, thus
    /// the final aligned query base will be at `ref_end - 1`.
    pub ref_end:    usize,
    /// Vector of query insertion relative to reference. Ordered by reference
    /// indices.
    pub insertions: Vec<SamInsertion>,
}

impl SamAligned {
    fn new(
        aligned: Vec<u8>, qaligned: Vec<u8>, ref_start: usize, ref_end: usize, insertions: Vec<SamInsertion>,
    ) -> SamAligned {
        SamAligned {
            aligned,
            qaligned,
            ref_start,
            ref_end,
            insertions,
        }
    }

    /// At the given reference index, provides the reference-aligned query's
    /// nucleotide and encoded (ASCII) quality score as a Optional tuple of
    /// `(base, qs)`, otherwise, `None` is returned.
    #[must_use]
    pub fn get_base_and_quality(&self, at_index: usize) -> Option<(u8, u8)> {
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
    pub fn get_base(&self, at_index: usize) -> Option<u8> {
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
    pub fn ref_range(&self) -> Range<usize> {
        self.ref_start..self.ref_end
    }

    /// Merges two aligned queries' reference ranges such that the combined
    /// range spans both aligned regions.
    #[must_use]
    pub fn merge_ref_range(&self, other: &SamAligned) -> Range<usize> {
        use std::cmp::{max, min};

        min(self.ref_start, other.ref_start)..max(self.ref_end, other.ref_end)
    }

    /// Checks if the `SamAligned` contains an insertion after the 0-based
    /// reference index and returns the query's unaligned range.
    #[must_use]
    pub fn get_insert_after(&self, reference_index: usize) -> Option<SamInsertion> {
        self.insertions
            .iter()
            .find(|insert| insert.ref_index == reference_index)
            .cloned()
    }
}

impl From<&SamData> for SamAligned {
    fn from(row: &SamData) -> Self {
        row.get_aligned()
    }
}

pub struct SAMReader<R: std::io::Read> {
    pub sam_reader: std::io::Lines<std::io::BufReader<R>>,
}

impl<R: std::io::Read> SAMReader<R> {
    pub fn new(inner: R) -> Self {
        SAMReader {
            sam_reader: std::io::BufReader::new(inner).lines(),
        }
    }
}

impl<R: std::io::Read> Iterator for SAMReader<R> {
    type Item = std::io::Result<SamRow>;

    fn next(&mut self) -> Option<Self::Item> {
        let Some(line) = self.sam_reader.next() else {
            return None;
        };

        match line {
            Ok(s) if s.starts_with('@') => Some(Ok(SamRow::Header(s))),
            Ok(s) => {
                let r: Vec<&str> = s.split('\t').collect();
                // TO-DO: account for optional fields
                if r.len() < 11 {
                    return Some(Err(IOError::new(
                        ErrorKind::InvalidData,
                        "SAM file did not have at least 11 fields!",
                    )));
                }

                let read_bit_flags = match r[1].parse::<u16>() {
                    Ok(f) => f,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM bit flag '{}'", &r[1]),
                        )))
                    }
                };

                // Per the SAM spec, the 1-based reference position
                let aligned_reference_position = match r[3].parse::<usize>() {
                    Ok(p) => p,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM reference position '{}'", &r[3]),
                        )))
                    }
                };

                let read_mapping_quality = match r[4].parse::<u8>() {
                    Ok(q) => q,
                    Err(e) => {
                        return Some(Err(IOError::new(
                            ErrorKind::InvalidData,
                            format!("Error: {e}, invalid SAM mapq '{}'", &r[4]),
                        )))
                    }
                };

                let row = SamData {
                    qname: r[0].to_owned(),
                    flag:  read_bit_flags,
                    rname: r[2].to_owned(),
                    pos:   aligned_reference_position,
                    mapq:  read_mapping_quality,
                    cigar: r[5].into(),
                    rnext: '*',
                    pnext: 0,
                    tlen:  0,
                    // TO-DO: Where should casing be handled?
                    seq:   r[9].as_bytes().to_ascii_uppercase().into(),
                    qual:  r[10].as_bytes().into(),
                };
                Some(Ok(SamRow::Data(Box::new(row))))
            }
            Err(e) => Some(Err(e)),
        }
    }
    // Push down predicates?
}

impl SAMReader<std::fs::File> {
    /// Reads a SAM text file into an iterator backed by a buffered reader.
    ///
    /// # Errors
    ///
    /// Will return `Err` if file or permissions or the like do not exist.
    pub fn from_filename<P>(filename: P) -> Result<SAMReader<File>, std::io::Error>
    where
        P: AsRef<Path>, {
        let file = File::open(filename)?;
        Ok(SAMReader::new(file))
    }
}

impl std::fmt::Display for SamRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SamRow::Data(d) => write!(f, "{d}"),
            SamRow::Header(h) => write!(f, "{h}"),
        }
    }
}

impl std::fmt::Display for SamData {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let SamData {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
        } = self;

        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
        )
    }
}
