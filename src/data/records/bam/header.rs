//! Functions for parsing SAM headers as needed for BAM serialization.

use crate::data::{
    bam::error::{BamEncodingError, BamError, BamHeaderError, BamRecordError, NumberSizeTarget},
    err::ResultWithErrorContext,
};
use std::{
    collections::{HashMap, HashSet, hash_map::Entry},
    io::Write,
};

/// Parsed [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) header needed
/// to write a [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) header.
///
/// All lines are preserved are normalized for the BAM header block and `@SQ`
/// lines are parsed into reference metadata for record encoding. Each `@SQ`
/// line must provide a unique `SN` value and a positive signed-32-bit `LN`
/// value.
#[derive(Default)]
pub(super) struct Header {
    /// Stored header lines, preserved as SAM text for BAM serialization.
    raw_lines: Vec<String>,
    /// Parsed reference records in BAM reference-ID order.
    refs:      Vec<ReferenceInfo>,
    /// A map from each `@SQ SN` value to its BAM reference ID.
    ref_to_id: HashMap<String, i32>,
}

impl Header {
    /// Parses and stores one SAM header line.
    ///
    /// Trailing line terminators are removed and embedded NUL bytes or line
    /// terminators are rejected. Non-`@SQ` lines are otherwise stored without
    /// semantic validation. `@SQ` lines are parsed into the BAM reference
    /// dictionary and validated.
    pub(super) fn push_header_line(&mut self, header_line: &str) -> Result<(), BamHeaderError> {
        let header_line = header_line.trim_end_matches(['\n', '\r']);
        if header_line.bytes().any(|byte| matches!(byte, b'\0' | b'\r' | b'\n')) {
            return Err(BamEncodingError::other("Header line contains an embedded NUL or line terminator").into());
        }

        if let Some(suffix) = header_line.strip_prefix("@SQ") {
            if !suffix.is_empty() && !suffix.starts_with('\t') {
                return Err(
                    BamEncodingError::other(format!("Invalid @SQ record type in header line {header_line:?}")).into(),
                );
            }

            let ref_info = parse_sq_line(header_line)?;

            if self.refs.len() >= i32::MAX as usize {
                return Err(BamEncodingError::SizeOverflow {
                    field:  "Number of references",
                    target: NumberSizeTarget::MaxExclusive(1usize << 31),
                }
                .into());
            }
            let idx = i32::try_from(self.refs.len()).expect("size checked");

            match self.ref_to_id.entry(ref_info.name.clone()) {
                Entry::Occupied(_) => {
                    return Err(BamHeaderError::DuplicateReference { name: ref_info.name });
                }
                Entry::Vacant(entry) => entry.insert(idx),
            };
            self.refs.push(ref_info);
            self.raw_lines.push(String::from(header_line));
        } else {
            self.raw_lines.push(String::from(header_line));
        }

        Ok(())
    }

    /// Returns the BAM reference ID associated with a SAM reference name.
    ///
    /// The SAM sentinel `*` is translated to BAM's `-1` unmapped reference ID.
    pub(super) fn get_ref_id(&self, rname: &str) -> Result<i32, BamRecordError> {
        if rname == "*" {
            return Ok(-1);
        }
        self.ref_to_id
            .get(rname)
            .copied()
            .ok_or_else(|| BamRecordError::ReferenceNotFound { name: rname.to_string() })
    }

    /// Serializes the BAM header text and reference dictionary from a
    /// [`Header`].
    pub(super) fn write_to<W: Write>(&self, writer: &mut W) -> Result<(), BamError> {
        /// BAM magic string.
        const BAM_MAGIC: &[u8; 4] = b"BAM\x01";

        let header_text_len = self.raw_lines.iter().try_fold(0usize, |total, line| {
            let line_len = line.len().checked_add(1).ok_or_else(|| {
                BamHeaderError::from(BamEncodingError::SizeOverflow {
                    field:  "header text",
                    target: NumberSizeTarget::MaxInclusive(usize::MAX),
                })
            })?;
            total.checked_add(line_len).ok_or_else(|| {
                BamHeaderError::from(BamEncodingError::SizeOverflow {
                    field:  "header text",
                    target: NumberSizeTarget::MaxInclusive(usize::MAX),
                })
            })
        })?;
        if header_text_len >= (1usize << 31) {
            return Err(BamHeaderError::from(BamEncodingError::SizeOverflow {
                field:  "Header text length",
                target: NumberSizeTarget::MaxExclusive(1usize << 31),
            })
            .into());
        }
        let l_text = u32::try_from(header_text_len).expect("size checked");
        if self.refs.len() >= (1usize << 31) {
            return Err(BamHeaderError::from(BamEncodingError::SizeOverflow {
                field:  "Number of references",
                target: NumberSizeTarget::MaxExclusive(1usize << 31),
            })
            .into());
        }
        let n_ref = u32::try_from(self.refs.len()).expect("size checked");

        writer.write_all(BAM_MAGIC).with_context("Cannot write BAM magic string")?;
        writer
            .write_all(&l_text.to_le_bytes())
            .with_context("Cannot write BAM header text length")?;
        for line in &self.raw_lines {
            writer
                .write_all(line.as_bytes())
                .with_context("Cannot write BAM header text")?;
            writer
                .write_all(b"\n")
                .with_context("Cannot write BAM header line terminator")?;
        }
        writer
            .write_all(&n_ref.to_le_bytes())
            .with_context("Cannot write BAM number of references")?;

        for reference in &self.refs {
            let name_len_with_nul = reference.name.len().checked_add(1).ok_or_else(|| {
                BamHeaderError::from(BamEncodingError::SizeOverflow {
                    field:  "reference name length",
                    target: NumberSizeTarget::MaxInclusive(usize::MAX),
                })
            })?;
            let l_name = u32::try_from(name_len_with_nul).map_err(|_| {
                BamHeaderError::from(BamEncodingError::SizeOverflow {
                    field:  "reference name",
                    target: NumberSizeTarget::MaxInclusive(u32::MAX as usize),
                })
            })?;

            writer
                .write_all(&l_name.to_le_bytes())
                .with_context("Cannot write BAM reference name length")?;
            writer
                .write_all(reference.name.as_bytes())
                .with_context("Cannot write BAM reference name")?;
            writer
                .write_all(&[0])
                .with_context("Cannot write BAM reference name null terminator")?;
            writer
                .write_all(&reference.len.to_le_bytes())
                .with_context("Cannot write BAM reference sequence length")?;
        }

        Ok(())
    }
}

/// Parsed `@SQ` reference entry.
struct ReferenceInfo {
    /// Reference sequence name.
    name: String,
    /// Reference sequence length.
    len:  u32,
}

/// Parses a SAM `@SQ` header line, requiring valid, unique `SN` and `LN`
/// fields.
///
/// Other `@SQ` tags are preserved in raw header text.
fn parse_sq_line(line: &str) -> Result<ReferenceInfo, BamHeaderError> {
    let mut fields = line.split('\t');
    // may be a duplicate check but we need to iterate past anyway
    if fields.next() != Some("@SQ") {
        return Err(BamEncodingError::other(format!("Invalid @SQ header line {line:?}")).into());
    }

    let mut sn = None;
    let mut ln = None;
    let mut seen_tags = HashSet::new();

    for field in fields {
        let (tag, value) = field
            .split_once(':')
            .ok_or_else(|| BamEncodingError::other(format!("Invalid @SQ field {field:?}")))?;

        let bytes = tag.as_bytes();
        if bytes.len() != 2 || !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphanumeric() {
            return Err(BamEncodingError::other(format!("Invalid @SQ tag {tag:?}")).into());
        }
        // values should not be large in length, simd should not be necessary
        if value.is_empty() || value.bytes().any(|byte| matches!(byte, b'\0' | b'\r' | b'\n')) {
            return Err(BamEncodingError::other(format!("Invalid @SQ value for tag {tag:?}")).into());
        }
        if !seen_tags.insert(tag) {
            return Err(BamEncodingError::other(format!("Duplicate @SQ tag {tag:?}")).into());
        }

        match tag {
            "SN" => sn = Some(String::from(value)),
            "LN" => ln = Some(parse_ln(value)?),
            _ => {}
        }
    }

    let name = sn.ok_or_else(|| BamEncodingError::other("@SQ line is missing required SN tag"))?;
    if !name.as_bytes().iter().all(u8::is_ascii_graphic)
        || matches!(name.as_bytes().first(), None | Some(b'*' | b'='))
        || name.bytes().any(|byte| {
            matches!(
                byte,
                b'\\' | b',' | b'"' | b'\'' | b'(' | b')' | b'[' | b']' | b'{' | b'}' | b'<' | b'>' | b'`'
            )
        })
    {
        return Err(BamEncodingError::other(format!("Invalid @SQ SN reference name {name:?}")).into());
    }
    let len = ln.ok_or_else(|| BamEncodingError::other("@SQ line is missing required LN tag"))?;

    Ok(ReferenceInfo { name, len })
}

/// Parses an `@SQ:LN` value, requiring a positive 32-bit length in SAM's
/// allowed signed 32-bit range.
fn parse_ln(text: &str) -> Result<u32, BamHeaderError> {
    let value = text
        .parse::<i32>()
        .map_err(|source| BamEncodingError::other_with_source(format!("Invalid integer for @SQ:LN: '{text}'"), source))?;
    if value < 1 {
        return Err(BamEncodingError::other(format!("Value {value} for @SQ:LN is outside [1, {}]", i32::MAX)).into());
    }

    Ok(value.cast_unsigned())
}
