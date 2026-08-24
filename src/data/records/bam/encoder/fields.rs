//! Encoding logic for individual BAM fields.

use crate::{
    DEFAULT_SIMD_LANES,
    alignment::AlignmentStates,
    data::{
        ByteIndexMap,
        bam::{
            encoder::MAX_CIGAR_INC,
            error::{BamEncodingError, BamRecordError},
        },
        cigar::CigarError,
        nucleotides::Nucleotides,
        phred::{QScoreInt, QualityScores, QualityScoresView},
        sam::{OptArray, SamOptField, SamOptRaw, SamOptValue, is_missing_sam_field},
        validation::CheckSequence,
        views::Len,
    },
    search::ByteSubstring,
};
use std::collections::{HashMap, hash_map::Entry};

/// Encodes a SAM `QNAME` as BAM's NUL-terminated read name field.
pub(super) fn encode_read_name(qname: &str) -> Result<Vec<u8>, BamRecordError> {
    let qname_bytes = qname.as_bytes();
    if qname_bytes.is_empty()
        || !qname_bytes.is_graphic_simd::<{ DEFAULT_SIMD_LANES }>()
        || qname_bytes.find_byte(b'@').is_some()
    {
        return Err(BamEncodingError::other("QNAME must match SAM/BAM [!-?A-~]{1,254}").into());
    }
    let total_len = qname.len().checked_add(1).ok_or(BamEncodingError::SizeOverflow {
        field:  "QNAME length",
        target: "usize",
    })?;
    if total_len > u8::MAX as usize {
        return Err(BamEncodingError::other(format!("QNAME is too long for BAM ({total_len} bytes including NUL)")).into());
    }

    let mut out = Vec::with_capacity(total_len);
    out.extend_from_slice(qname.as_bytes());
    out.push(0);
    Ok(out)
}

/// Encodes a [`Nucleotides`] sequence in BAM's packed 4-bit representation.
///
/// A missing sequence (`*`) or empty sequence is encoded as an empty byte
/// vector. Bases are packed two per byte, with the first base in the high
/// nibble and the second base in the low nibble. The BAM alphabet
/// `=ACMGRSVTWYHKDBN` is used; `U` is encoded as `T`, and any other byte is
/// encoded as `N`. Treating `U` as `T` is an intentional *Zoe* choice that
/// deviates from the SAM/BAM spec; the spec's BAM sequence alphabet does not
/// include a `U` code.
pub(super) fn encode_seq(seq: &Nucleotides) -> Vec<u8> {
    /// Sequence byte index map: `=ACMGRSVTWYHKDBN` are mapped to `[0,15]`. `N`
    /// is used as a catch-all for all other characters, and `U` is encoded as
    /// `T` as an intentional deviation from the BAM sequence alphabet. `=` is
    /// called a base code in the specs.
    const SEQ_MAP: ByteIndexMap<16> =
        ByteIndexMap::new_ignoring_case(*b"=ACMGRSVTWYHKDBN", b'N').add_synonym_ignore_case(b'U', b'T');

    if seq.as_bytes() == b"*" || seq.is_empty() {
        return Vec::new();
    }

    let mut out = Vec::with_capacity(seq.len().div_ceil(2));
    let mut iter = seq.iter().copied();
    while let Some(base_hi) = iter.next() {
        let hi = SEQ_MAP[base_hi] << 4;
        let lo = iter.next().map_or(0, |b| SEQ_MAP[b]);
        out.push(hi | lo);
    }
    out
}

/// Encodes [`QualityScores`] as BAM quality bytes.
///
/// BAM stores raw Phred scores, so the ASCII `+33` offset is removed. A SAM
/// quality value that is missing (`*` or empty) becomes `0xFF` repeated once
/// per sequence base. If `l_seq` is zero, this always returns an empty byte
/// vector.
pub(super) fn encode_qual(qual: &QualityScores, l_seq: usize) -> Result<Vec<u8>, BamRecordError> {
    if l_seq == 0 {
        return Ok(Vec::new());
    }
    if is_missing_sam_field(qual) {
        return Ok(vec![0xFF; l_seq]);
    }

    let qual_view = QualityScoresView::try_from(qual.as_bytes())
        .map_err(|source| BamEncodingError::other_with_source("Quality scores cannot be encoded as BAM", source))?;
    Ok(qual_view.iter().map(|&byte| QScoreInt::from(byte).as_u8()).collect())
}

/// Encodes parsed [`AlignmentStates`] as BAM CIGAR words.
///
/// Each word stores a 28-bit operation length and a 4-bit operation code.
pub(super) fn encode_cigar(ciglets: &AlignmentStates) -> Result<Vec<u32>, BamRecordError> {
    /// CIGAR byte index map: `MIDNSHP=X`→`012345678`. `?` is used as a
    /// catch-all for invalid CIGAR operations.
    const CIGAR_MAP: ByteIndexMap<10> = ByteIndexMap::new(*b"MIDNSHP=X?", b'?');

    let mut encoded = Vec::with_capacity(ciglets.len());
    for ciglet in ciglets {
        let cig_inc = u32::try_from(ciglet.inc).map_err(|_| BamEncodingError::SizeOverflow {
            field:  "CIGAR increment length",
            target: "u32",
        })?;
        let mapped_op = CIGAR_MAP[ciglet.op];
        if mapped_op == CIGAR_MAP[b'?'] {
            return Err(BamRecordError::InvalidCigar {
                source: CigarError::InvalidOperation,
            });
        }
        if cig_inc > MAX_CIGAR_INC {
            return Err(BamEncodingError::SizeOverflow {
                field:  "CIGAR increment length",
                target: "28-bit BAM CIGAR length field",
            }
            .into());
        }
        encoded.push(cig_inc << 4 | u32::from(mapped_op));
    }
    Ok(encoded)
}

/// Encodes SAM optional fields in BAM aux-field format.
///
/// When `cg_field` is present, the encoded output is extended with a `CG:B:I`
/// field containing the full CIGAR for records that overflow BAM's inline
/// `n_cigar_op` limit. The `CG` tag is reserved for long CIGARs.
pub(super) fn encode_aux_fields(aux_fields: &SamOptRaw, cg_field: Option<&[u32]>) -> Result<Vec<u8>, BamRecordError> {
    let mut out = Vec::new();
    let mut seen_tags: HashMap<[u8; 2], usize> = HashMap::with_capacity(aux_fields.len());
    let mut unique_fields: Vec<SamOptField> = Vec::with_capacity(aux_fields.len());
    for field in aux_fields.iter() {
        let field = field
            .map_err(|source| BamEncodingError::other_with_source("SAM optional field cannot be encoded as BAM", source))?;
        if &field.tag == b"CG" {
            return Err(
                BamEncodingError::other("Auxiliary tag CG is reserved for encoder-generated long-CIGAR data").into(),
            );
        }
        match seen_tags.entry(field.tag) {
            Entry::Occupied(entry) => {
                let field_idx = *entry.get();
                if unique_fields[field_idx].value != field.value {
                    return Err(BamEncodingError::other(format!(
                        "Conflicting duplicate auxiliary tag {}{}",
                        field.tag[0] as char, field.tag[1] as char
                    ))
                    .into());
                }
            }
            Entry::Vacant(entry) => {
                entry.insert(unique_fields.len());
                unique_fields.push(field);
            }
        }
    }

    for field in &unique_fields {
        encode_single_aux(field, &mut out)?;
    }

    if let Some(cigar) = cg_field {
        out.extend_from_slice(b"CG");
        out.push(b'B');
        out.push(b'I');
        out.extend_from_slice(
            &u32::try_from(cigar.len())
                .map_err(|_| BamEncodingError::SizeOverflow {
                    field:  "CG array",
                    target: "u32",
                })?
                .to_le_bytes(),
        );
        for cig in cigar {
            out.extend_from_slice(&cig.to_le_bytes());
        }
    }

    Ok(out)
}

/// Encodes one parsed SAM optional field in BAM aux-field format and appends it
/// to `out`.
fn encode_single_aux(field: &SamOptField, out: &mut Vec<u8>) -> Result<(), BamRecordError> {
    let tag = field.tag;
    out.extend_from_slice(&tag);

    match &field.value {
        SamOptValue::Char(c) => {
            out.push(b'A');
            out.push(*c);
        }
        SamOptValue::Int(i) => {
            encode_bam_aux_int(*i, out)?;
        }
        SamOptValue::Float(f) => {
            out.push(b'f');
            out.extend_from_slice(&f.to_le_bytes());
        }
        SamOptValue::String(s) => {
            out.push(b'Z');
            out.extend_from_slice(s.as_bytes());
            out.push(0);
        }
        SamOptValue::Hex(h) => {
            out.push(b'H');
            out.extend_from_slice(h.as_bytes());
            out.push(0);
        }
        SamOptValue::Array(array) => {
            encode_bam_aux_array(array, out)?;
        }
    }
    Ok(())
}

/// Encodes an integer auxiliary value using the narrowest BAM integer type
/// that can represent it.
///
/// Non-negative values are encoded as unsigned `C`, `S`, or `I` fields.
/// Negative values are encoded as signed `c`, `s`, or `i` fields.
fn encode_bam_aux_int(value: i64, out: &mut Vec<u8>) -> Result<(), BamRecordError> {
    if value < 0 {
        if let Ok(value) = i8::try_from(value) {
            out.push(b'c');
            out.push(value.cast_unsigned());
        } else if let Ok(value) = i16::try_from(value) {
            out.push(b's');
            out.extend_from_slice(&value.to_le_bytes());
        } else if let Ok(value) = i32::try_from(value) {
            out.push(b'i');
            out.extend_from_slice(&value.to_le_bytes());
        } else {
            return Err(BamEncodingError::other(format!(
                "Integer optional field value {value} is outside BAM's encodable range"
            ))
            .into());
        }
    } else {
        if let Ok(value) = u8::try_from(value) {
            out.push(b'C');
            out.push(value);
        } else if let Ok(value) = u16::try_from(value) {
            out.push(b'S');
            out.extend_from_slice(&value.to_le_bytes());
        } else if let Ok(value) = u32::try_from(value) {
            out.push(b'I');
            out.extend_from_slice(&value.to_le_bytes());
        } else {
            return Err(BamEncodingError::other(format!(
                "Integer optional field value {value} is outside BAM's encodable range"
            ))
            .into());
        }
    }

    Ok(())
}

/// Encodes an array of auxiliary values.
fn encode_bam_aux_array(array: &OptArray, out: &mut Vec<u8>) -> Result<(), BamRecordError> {
    out.push(b'B');
    match array {
        OptArray::I8(values) => write_opt_array(out, b'c', values, values.iter().map(|v| v.cast_unsigned()))?,
        OptArray::U8(values) => write_opt_array(out, b'C', values, values.iter().copied())?,
        OptArray::I16(values) => write_opt_array(out, b's', values, values.iter().flat_map(|v| v.to_le_bytes()))?,
        OptArray::U16(values) => write_opt_array(out, b'S', values, values.iter().flat_map(|v| v.to_le_bytes()))?,
        OptArray::I32(values) => write_opt_array(out, b'i', values, values.iter().flat_map(|v| v.to_le_bytes()))?,
        OptArray::U32(values) => write_opt_array(out, b'I', values, values.iter().flat_map(|v| v.to_le_bytes()))?,
        OptArray::F32(values) => write_opt_array(out, b'f', values, values.iter().flat_map(|v| v.to_le_bytes()))?,
    }

    Ok(())
}

/// Helper functions for encoding auxiliary arrays.
fn write_opt_array<V>(
    out: &mut Vec<u8>, subtype: u8, values: &[V], bytes: impl IntoIterator<Item = u8>,
) -> Result<(), BamRecordError> {
    out.push(subtype);
    out.extend_from_slice(
        &u32::try_from(values.len())
            .map_err(|_| BamEncodingError::SizeOverflow {
                field:  "'B' field array",
                target: "u32",
            })?
            .to_le_bytes(),
    );

    out.extend(bytes);

    Ok(())
}
