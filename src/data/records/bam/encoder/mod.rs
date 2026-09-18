//! Encode SAM records into binary alignment blocks for BAM output.
//!
//! This layer validates the SAM-side values that have stricter BAM
//! representations, prepares the variable-width record payloads, computes BAM
//! bins, normalizes BAM flag bits, and handles BAM's long-CIGAR representation
//! when the inline CIGAR field would overflow.

use crate::data::{
    bam::{
        encoder::{
            binning::compute_bin,
            fields::{encode_aux_fields, encode_cigar, encode_qual, encode_read_name, encode_seq},
        },
        error::{BamEncodingError, BamError, BamRecordError, NumberSizeTarget},
        header::Header,
    },
    cigar::ToCigletIterator,
    sam::{Flag, GetSamFields, is_missing_sam_field},
    views::Len,
};
use std::ops::Range;

mod binning;
mod fields;

/// Maximum operation length that fits in BAM's 28-bit CIGAR increment field.
const MAX_CIGAR_INC: u32 = 0x0FFF_FFFF;

/// Fully prepared BAM alignment record ready for binary serialization.
pub(super) struct PreparedBamRecord {
    /// BAM reference ID for `RNAME`, or `-1` for the SAM sentinel `*`.
    pub ref_id:       i32,
    /// Zero-based leftmost mapping position, or `-1` when `RNAME` is `*` or
    /// SAM `POS` is `0`.
    pub pos0:         i32,
    /// Mapping quality copied directly from the SAM record.
    mapq:             u8,
    /// Hierarchical BAM bin computed from the reference interval, or the
    /// reserved unplaced-unmapped bin when no coordinate is available.
    pub bin:          u16,
    /// Checked half-open reference interval used for BAM/BAI indexing.
    pub ref_interval: Option<Range<u32>>,
    /// Number of inline CIGAR operations stored in `cigar_field`.
    n_cigar_op:       u16,
    /// BAM flag word after normalizing SAM flags for *Zoe*'s BAM output.
    flag:             Flag,
    /// Query sequence length stored in the BAM core record.
    l_seq:            u32,
    /// NUL-terminated BAM read name payload.
    read_name:        Vec<u8>,
    /// BAM-encoded CIGAR words, or the 2-op placeholder used for long CIGARs.
    cigar_field:      Vec<u32>,
    /// Query sequence encoded in BAM's packed 4-bit nucleotide representation.
    seq:              Vec<u8>,
    /// Query quality scores encoded as raw Phred bytes, or `0xFF` for missing
    /// quality scores.
    qual:             Vec<u8>,
    /// BAM auxiliary fields, including synthesized `CG:B:I` when needed.
    aux:              Vec<u8>,
}

impl PreparedBamRecord {
    /// Validates a [`GetSamFields`] record and prepares all BAM-encoded fields.
    ///
    /// ## Errors
    ///
    /// Returns an error if the SAM record cannot be represented in BAM, if its
    /// CIGAR or auxiliary fields cannot be parsed, if an option field contains
    /// the reserved long CIGAR `CG` tag, or if any encoded field would overflow
    /// BAM's size limits.
    pub(super) fn new(header: &Header, data: impl GetSamFields) -> Result<Self, BamRecordError> {
        /// The maximum inclusive position allowed by the SAM file format.
        const MAX_SAM_POS: usize = i32::MAX as usize;

        let seq = data.seq().filter(|s| !is_missing_sam_field(s));
        let qual = data.qual().filter(|q| !is_missing_sam_field(q));
        let cigar = data.cigar().filter(|c| c.to_ciglet_iterator_checked().next().is_some());

        let ref_id = header.get_ref_id(data.rname())?;
        let pos = data.pos().unwrap_or(0);

        // 0-indexed
        let pos0 = match pos {
            0 => -1,
            1..=MAX_SAM_POS => i32::try_from(pos - 1).expect("range validated"),
            _ => {
                return Err(BamEncodingError::SizeOverflow {
                    field:  "SAM POS",
                    target: NumberSizeTarget::MaxInclusive(i32::MAX as u64),
                }
                .into());
            }
        };

        let read_name = encode_read_name(data.qname())?;
        let encoded_seq = encode_seq(seq);
        let encoded_qual = encode_qual(qual, seq)?;
        let (encoded_cigar, spans) = encode_cigar(cigar.as_ref())?;
        if let Some(spans) = &spans
            && let Some(seq) = seq
            && spans.query_span != seq.len()
        {
            return Err(BamEncodingError::other(format!(
                "SEQ length ({seq_len}) does not match query-consuming CIGAR length ({query_span})",
                seq_len = seq.len(),
                query_span = spans.query_span
            ))
            .into());
        }

        let flag = normalize_flags(data.flag(), cigar.is_none());

        let l_seq = seq.map_or(Ok(0), |seq| {
            u32::try_from(seq.len()).map_err(|_| BamEncodingError::SizeOverflow {
                field:  "SEQ",
                target: NumberSizeTarget::MaxInclusive(u64::from(u32::MAX)),
            })
        })?;

        let long_cigar = encoded_cigar.len() > u16::MAX as usize;

        let (aux, cigar_field, n_cigar_op) = if let Some(spans) = &spans
            && long_cigar
        {
            if l_seq > MAX_CIGAR_INC {
                return Err(BamEncodingError::SizeOverflow {
                    field:  "long-CIGAR placeholder sequence length",
                    target: NumberSizeTarget::MaxExclusive(1 << 28),
                }
                .into());
            }

            if spans.ref_span > MAX_CIGAR_INC {
                return Err(BamEncodingError::SizeOverflow {
                    field:  "long-CIGAR placeholder reference span",
                    target: NumberSizeTarget::MaxExclusive(1 << 28),
                }
                .into());
            }

            let aux = encode_aux_fields(data.opt_fields(), Some(&encoded_cigar))?;
            let cigar_field = vec![(l_seq << 4) | 4, (spans.ref_span << 4) | 3];
            let n_cigar_op = 2;

            (aux, cigar_field, n_cigar_op)
        } else {
            let aux = encode_aux_fields(data.opt_fields(), None)?;
            let n_cigar_op = u16::try_from(encoded_cigar.len()).map_err(|_| BamEncodingError::SizeOverflow {
                field:  "number of CIGAR ops",
                target: NumberSizeTarget::MaxInclusive(u64::from(u16::MAX)),
            })?;

            (aux, encoded_cigar, n_cigar_op)
        };

        let ref_interval = checked_interval(pos0, spans.as_ref().map(|spans| spans.ref_span), flag)?;
        let bin = compute_bin(ref_interval.as_ref())?;

        Ok(Self {
            ref_id,
            pos0,
            mapq: data.mapq().unwrap_or(255),
            bin,
            ref_interval,
            n_cigar_op,
            flag,
            l_seq,
            read_name,
            cigar_field,
            seq: encoded_seq,
            qual: encoded_qual,
            aux,
        })
    }

    /// Encodes the prepared record as one complete BAM alignment block.
    ///
    /// BAM mate/reference-next fields are currently emitted with Zoe's
    /// placeholder values: `RNEXT = -1`, `PNEXT = -1`, and `TLEN = 0`.
    ///
    /// ## Errors
    ///
    /// Returns an error if the total block size overflows BAM's representable
    /// limits.
    pub(super) fn encode(&self, qname: Option<&str>) -> Result<Vec<u8>, BamError> {
        /// Placeholder mate reference ID written because *Zoe* does not
        /// currently preserve SAM `RNEXT`.
        const RNEXT: i32 = -1;
        /// Placeholder mate position written because *Zoe* does not currently
        /// preserve SAM `PNEXT`.
        const PNEXT: i32 = -1;
        /// Placeholder template length written because *Zoe* does not currently
        /// preserve SAM `TLEN`.
        const TLEN: i32 = 0;
        /// Total length of BAM fields with fixed length `refID` through `tlen`.
        const BAM_RECORD_CORE_SIZE: usize = 32;

        let cigar_size = self
            .cigar_field
            .len()
            .checked_mul(std::mem::size_of::<u32>())
            .ok_or_else(|| {
                BamError::record(
                    qname,
                    BamEncodingError::SizeOverflow {
                        field:  "BAM CIGAR",
                        target: NumberSizeTarget::MaxInclusive(usize::MAX as u64),
                    },
                )
            })?;
        let block_size = bam_block_size(&[
            BAM_RECORD_CORE_SIZE,
            self.read_name.len(),
            cigar_size,
            self.seq.len(),
            self.qual.len(),
            self.aux.len(),
        ])
        .map_err(|source| BamError::record(qname, source))?;

        let block_size_u32 = u32::try_from(block_size).map_err(|_| {
            BamError::record(
                qname,
                BamEncodingError::SizeOverflow {
                    field:  "BAM block",
                    target: NumberSizeTarget::MaxInclusive(u64::from(u32::MAX)),
                },
            )
        })?;
        let capacity = block_size.checked_add(4).ok_or_else(|| {
            BamError::record(
                qname,
                BamEncodingError::SizeOverflow {
                    field:  "BAM block buffer",
                    target: NumberSizeTarget::MaxInclusive(usize::MAX as u64),
                },
            )
        })?;
        let mut buf = Vec::with_capacity(capacity);

        buf.extend_from_slice(&block_size_u32.to_le_bytes());
        buf.extend_from_slice(&self.ref_id.to_le_bytes());
        buf.extend_from_slice(&self.pos0.to_le_bytes());
        buf.push(u8::try_from(self.read_name.len()).map_err(|_| {
            BamError::record(
                qname,
                BamEncodingError::SizeOverflow {
                    field:  "QNAME length",
                    target: NumberSizeTarget::MaxInclusive(u64::from(u8::MAX)),
                },
            )
        })?);
        buf.push(self.mapq);
        buf.extend_from_slice(&self.bin.to_le_bytes());
        buf.extend_from_slice(&self.n_cigar_op.to_le_bytes());
        buf.extend_from_slice(&self.flag.0.to_le_bytes());
        buf.extend_from_slice(&self.l_seq.to_le_bytes());
        // RNEXT, PNEXT, and TLEN are not populated yet.
        buf.extend_from_slice(&RNEXT.to_le_bytes());
        buf.extend_from_slice(&PNEXT.to_le_bytes());
        buf.extend_from_slice(&TLEN.to_le_bytes());
        buf.extend_from_slice(&self.read_name);
        for cig in &self.cigar_field {
            buf.extend_from_slice(&cig.to_le_bytes());
        }
        buf.extend_from_slice(&self.seq);
        buf.extend_from_slice(&self.qual);
        buf.extend_from_slice(&self.aux);

        Ok(buf)
    }
}

/// Normalizes the SAM flag word for the BAM record *Zoe* can faithfully emit.
///
/// Unsupported or mate-dependent bits are cleared, and records with an empty
/// CIGAR are marked as unmapped so the serialized BAM flag is consistent with
/// *Zoe*'s treatment of empty-CIGAR records.
fn normalize_flags(flag: Option<Flag>, cigar_is_missing: bool) -> Flag {
    let mut flag = flag.unwrap_or_default();

    flag = filter_unsupported_flags(flag);

    if cigar_is_missing {
        flag.set_unmapped();
    }

    flag
}

/// Removes SAM/BAM flag bits whose meaning *Zoe* cannot preserve in BAM output.
///
/// Cleared bits include paired-template and mate-dependent flags such as `0x1`
/// (multiple segments), `0x2` (properly aligned), `0x8` (next segment
/// unmapped), `0x20` (next segment reverse-complemented), `0x40` (first
/// segment), and `0x80` (last segment). Reserved or otherwise unknown bits are
/// also cleared.
fn filter_unsupported_flags(mut flag: Flag) -> Flag {
    flag.unset_segmented();
    flag.unset_properly_segmented();
    flag.unset_unmapped_next_segment();
    flag.unset_revcomp_next_segment();
    flag.unset_first_template();
    flag.unset_last_template();
    flag.standardize();

    flag
}

/// Computes the checked total BAM block size from a list of encoded field
/// sizes.
fn bam_block_size(component_sizes: &[usize]) -> Result<usize, BamRecordError> {
    component_sizes.iter().try_fold(0usize, |total, size| {
        total.checked_add(*size).ok_or_else(|| {
            BamRecordError::from(BamEncodingError::SizeOverflow {
                field:  "BAM block",
                target: NumberSizeTarget::MaxInclusive(usize::MAX as u64),
            })
        })
    })
}

/// Computes the effective reference interval used by BAM binning and indexing.
///
/// Records without a coordinate return `Ok(None)`. Coordinate-bearing records
/// are checked against the classic BAI coordinate limit.
fn checked_interval(pos0: i32, ref_span: Option<u32>, flag: Flag) -> Result<Option<Range<u32>>, BamRecordError> {
    /// Exclusive upper bound for coordinates representable by the classic BAI
    /// index.
    const BAI_COORDINATE_LIMIT: u32 = 1 << 29;
    let beg = match pos0 {
        ..=-2 => return Err(BamRecordError::BinningOutOfRange),
        -1 => return Ok(None),
        pos0 @ 0..=i32::MAX => pos0.cast_unsigned(),
    };

    let effective_span = if flag.is_unmapped() { 1 } else { ref_span.unwrap_or(1).max(1) };

    let end = beg.checked_add(effective_span).ok_or(BamRecordError::BinningOutOfRange)?;

    if end > BAI_COORDINATE_LIMIT {
        return Err(BamRecordError::BinningOutOfRange);
    }

    Ok(Some(beg..end))
}
