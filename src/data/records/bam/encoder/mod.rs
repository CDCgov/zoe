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
            fields::{CigarSpans, encode_aux_fields, encode_cigar, encode_qual, encode_read_name, encode_seq},
        },
        error::{BamEncodingError, BamError, BamRecordError, NumberSizeTarget},
        header::Header,
    },
    sam::{SamData, is_missing_sam_field},
    views::Len,
};

mod binning;
mod fields;

/// Maximum operation length that fits in BAM's 28-bit CIGAR increment field.
const MAX_CIGAR_INC: u32 = 0x0FFF_FFFF;

/// SAM/BAM flag bit indicating that the read itself is unmapped.
const READ_UNMAPPED: u16 = 0x4;

/// Fully prepared BAM alignment record ready for binary serialization.
pub(super) struct PreparedBamRecord {
    /// BAM reference ID for `RNAME`, or `-1` for the SAM sentinel `*`.
    ref_id:      i32,
    /// Zero-based leftmost mapping position, or `-1` when `RNAME` is `*` or
    /// SAM `POS` is `0`.
    pos0:        i32,
    /// Mapping quality copied directly from the SAM record.
    mapq:        u8,
    /// Hierarchical BAM bin computed from the reference interval, or the
    /// reserved unplaced-unmapped bin when no coordinate is available.
    bin:         u16,
    /// Number of inline CIGAR operations stored in `cigar_field`.
    n_cigar_op:  u16,
    /// BAM flag word after normalizing SAM flags for *Zoe*'s BAM output.
    flag:        u16,
    /// Query sequence length stored in the BAM core record.
    l_seq:       u32,
    /// NUL-terminated BAM read name payload.
    read_name:   Vec<u8>,
    /// BAM-encoded CIGAR words, or the 2-op placeholder used for long CIGARs.
    cigar_field: Vec<u32>,
    /// Query sequence encoded in BAM's packed 4-bit nucleotide representation.
    seq:         Vec<u8>,
    /// Query quality scores encoded as raw Phred bytes, or `0xFF` for missing
    /// quality scores.
    qual:        Vec<u8>,
    /// BAM auxiliary fields, including synthesized `CG:B:I` when needed.
    aux:         Vec<u8>,
}

impl PreparedBamRecord {
    /// Validates a [`SamData`] record and prepares all BAM-encoded fields.
    ///
    /// ## Errors
    ///
    /// Returns an error if the SAM record cannot be represented in BAM, if its
    /// CIGAR or auxiliary fields cannot be parsed, if an option field contains
    /// the reserved long CIGAR `CG` tag, or if any encoded field would overflow
    /// BAM's size limits.
    #[allow(clippy::too_many_lines)]
    pub(super) fn new(header: &Header, data: &SamData) -> Result<Self, BamRecordError> {
        /// The maximum inclusive position allowed by the SAM file format.
        const MAX_SAM_POS: usize = i32::MAX as usize;

        let ref_id = header.get_ref_id(&data.rname)?;
        // 0-indexed
        let pos0 = match data.pos {
            0 => -1,
            1..=MAX_SAM_POS => i32::try_from(data.pos - 1).expect("range validated"),
            _ => {
                return Err(BamEncodingError::SizeOverflow {
                    field:  "SAM POS",
                    target: NumberSizeTarget::MaxInclusive(i32::MAX as usize),
                }
                .into());
            }
        };
        let seq_missing = is_missing_sam_field(&data.seq);
        let qual_missing = is_missing_sam_field(&data.qual);
        let l_seq = if seq_missing { 0 } else { data.seq.len() };

        let cigar_is_missing = data.cigar.is_empty();

        if seq_missing && !qual_missing {
            return Err(BamEncodingError::other("QUAL must be missing when SEQ is missing").into());
        }
        if !qual_missing && data.qual.len() != l_seq {
            return Err(BamEncodingError::other(format!(
                "QUAL length ({qual_len}) does not match SEQ length ({l_seq})",
                qual_len = data.qual.len(),
            ))
            .into());
        }

        let read_name = encode_read_name(&data.qname)?;
        let seq = encode_seq(&data.seq);
        let qual = encode_qual(&data.qual, l_seq)?;
        let (encoded_cigar, CigarSpans { query_span, ref_span }) = encode_cigar(&data.cigar)?;

        if !cigar_is_missing && !seq_missing && query_span != l_seq {
            return Err(BamEncodingError::other(format!(
                "SEQ length ({l_seq}) does not match query-consuming CIGAR length ({query_span})"
            ))
            .into());
        }

        let long_cigar = encoded_cigar.len() > u16::MAX as usize;
        let flag = normalize_flags(data.flag, cigar_is_missing);
        let l_seq = u32::try_from(l_seq).map_err(|_| BamEncodingError::SizeOverflow {
            field:  "SEQ",
            target: NumberSizeTarget::MaxInclusive(u32::MAX as usize),
        })?;

        if long_cigar {
            if l_seq > MAX_CIGAR_INC {
                return Err(BamEncodingError::SizeOverflow {
                    field:  "long-CIGAR placeholder sequence length",
                    target: NumberSizeTarget::MaxExclusive(1usize << 28),
                }
                .into());
            }

            if ref_span > MAX_CIGAR_INC {
                return Err(BamEncodingError::SizeOverflow {
                    field:  "long-CIGAR placeholder reference span",
                    target: NumberSizeTarget::MaxExclusive(1usize << 28),
                }
                .into());
            }
        }

        let aux = encode_aux_fields(&data.opt_fields, if long_cigar { Some(&encoded_cigar) } else { None })?;

        let (cigar_field, n_cigar_op) = if long_cigar {
            (vec![(l_seq << 4) | 4, (ref_span << 4) | 3], 2)
        } else {
            let n = u16::try_from(encoded_cigar.len()).map_err(|_| BamEncodingError::SizeOverflow {
                field:  "number of CIGAR ops",
                target: NumberSizeTarget::MaxInclusive(u16::MAX as usize),
            })?;
            (encoded_cigar, n)
        };

        let bin = compute_bin(pos0, ref_span, flag)?;

        Ok(Self {
            ref_id,
            pos0,
            mapq: data.mapq,
            bin,
            n_cigar_op,
            flag,
            l_seq,
            read_name,
            cigar_field,
            seq,
            qual,
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
    pub(super) fn encode(&self, qname: &str) -> Result<Vec<u8>, BamError> {
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
                        target: NumberSizeTarget::MaxInclusive(usize::MAX),
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
                    target: NumberSizeTarget::MaxInclusive(u32::MAX as usize),
                },
            )
        })?;
        let capacity = block_size.checked_add(4).ok_or_else(|| {
            BamError::record(
                qname,
                BamEncodingError::SizeOverflow {
                    field:  "BAM block buffer",
                    target: NumberSizeTarget::MaxInclusive(usize::MAX),
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
                    target: NumberSizeTarget::MaxInclusive(u8::MAX as usize),
                },
            )
        })?);
        buf.push(self.mapq);
        buf.extend_from_slice(&self.bin.to_le_bytes());
        buf.extend_from_slice(&self.n_cigar_op.to_le_bytes());
        buf.extend_from_slice(&self.flag.to_le_bytes());
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
fn normalize_flags(flag: u16, cigar_is_missing: bool) -> u16 {
    let mut flag = filter_unsupported_flags(flag);

    if cigar_is_missing {
        flag |= READ_UNMAPPED;
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
fn filter_unsupported_flags(flag: u16) -> u16 {
    const READ_REVERSE_COMPLEMENTED: u16 = 0x10;
    const SECONDARY_ALIGNMENT: u16 = 0x100;
    const FAILED_QUALITY_CHECKS: u16 = 0x200;
    const DUPLICATE: u16 = 0x400;
    const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

    const SUPPORTED_FLAGS: u16 = READ_UNMAPPED
        | READ_REVERSE_COMPLEMENTED
        | SECONDARY_ALIGNMENT
        | FAILED_QUALITY_CHECKS
        | DUPLICATE
        | SUPPLEMENTARY_ALIGNMENT;

    flag & SUPPORTED_FLAGS
}

/// Computes the checked total BAM block size from a list of encoded field
/// sizes.
fn bam_block_size(component_sizes: &[usize]) -> Result<usize, BamRecordError> {
    component_sizes.iter().try_fold(0usize, |total, size| {
        total.checked_add(*size).ok_or_else(|| {
            BamRecordError::from(BamEncodingError::SizeOverflow {
                field:  "BAM block",
                target: NumberSizeTarget::MaxInclusive(usize::MAX),
            })
        })
    })
}
