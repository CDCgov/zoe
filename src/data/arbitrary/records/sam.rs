//! Arbitrary implementations and specification structs for the [`SamData`]
//! record type and related optional SAM fields.
//!
//! The direct [`Arbitrary`] implementations preserve broad representational
//! coverage, while the specification structs provide optional SAM- and
//! BAM-oriented guarantees.

use crate::{
    data::{
        arbitrary::{
            ArbitrarySpecs, ByteSet, CigarSpecs, ClampAlignment, FloatSpecs, NucleotidesSpecs, QualityScoresSpecs,
            StringSpecs, VecSpecs,
        },
        cigar::{Cigar, LenInAlignment},
        sam::{OptArray, SamData, SamOptField, SamOptRaw, SamOptValue},
    },
    prelude::{Len, Nucleotides, QualityScores},
};
use arbitrary::{Arbitrary, Unstructured};
use std::num::NonZeroUsize;

impl<'a> Arbitrary<'a> for SamData {
    /// Generates an arbitrary [`SamData`] record from the given unstructured
    /// data.
    ///
    /// This ensures that `rnext` is `*`, `pnext` is 0, and `tlen` is 0.
    /// Optional fields are generated with [`SamOptRaw::arbitrary`] and may be
    /// malformed.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        let mut data = SamData::new(
            String::arbitrary(u)?,
            u16::arbitrary(u)?,
            String::arbitrary(u)?,
            usize::arbitrary(u)?,
            u8::arbitrary(u)?,
            Cigar::arbitrary(u)?,
            Nucleotides::arbitrary(u)?,
            QualityScores::arbitrary(u)?,
        );
        data.opt_fields = SamOptRaw::arbitrary(u)?;
        Ok(data)
    }
}

/// Specifications for generating arbitrary [`SamData`] records.
///
/// All generated records have `rnext` as `*`, `pnext` as 0, and `tlen` as 0.
/// Optional fields are generated according to [`SamDataSpecs::opt_fields`].
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct SamDataSpecs {
    /// The specifications for generating the `qname` field.
    pub qname: StringSpecs,

    /// The specifications for generating the `rname` field.
    pub rname: StringSpecs,

    /// The specifications for generating the `cigar` field.
    pub cigar: CigarSpecs,

    /// The specifications for generating the `seq` field.
    pub seq: NucleotidesSpecs,

    /// The specifications for generating the `qual` field.
    pub qual: QualityScoresSpecs,

    /// The specifications for generating the `opt_fields` field.
    pub opt_fields: SamOptRawSpecs,

    /// Whether to ensure that the `pos` field is nonzero (since it represents a
    /// 1-based position).
    pub nonzero_pos: bool,

    /// Ensures that the `pos` field is generated such that the end position of
    /// the alignment in the reference does not overflow.
    ///
    /// This assumes the match length of the CIGAR string also does not overflow
    /// (and that it does not equal [`usize::MAX`] when `nonzero_pos` is also
    /// specified).
    pub cap_end_pos: bool,

    /// Ensures that the length of the `seq` field agrees with the CIGAR string.
    ///
    /// As a side effect, this will cause the CIGAR string to only contain valid
    /// ciglets (the CIGAR string is iterated over and reconstructed with
    /// [`CigletIterator`], so see its documentation for more details).
    ///
    /// This may also cause the `seq` field to get shrunk. This may
    /// hypothetically invalidate arbitrary assumptions imposed by
    /// [`NucleotidesSpecs`].
    ///
    /// [`CigletIterator`]: crate::data::types::cigar::CigletIterator
    pub correct_seq_len: bool,

    /// Ensures that the length of the `qual` field agrees with the `seq` field.
    pub correct_qual_len: bool,
}

impl<'a> ArbitrarySpecs<'a> for SamDataSpecs {
    type Output = SamData;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        let mut cigar = self.cigar.make_arbitrary(u)?;
        let mut seq = self.seq.make_arbitrary(u)?;

        if self.correct_seq_len {
            cigar.clamp_query_len(seq.len());
            seq.shorten_to(cigar.query_len_in_alignment());
        }

        let qual = if self.correct_qual_len {
            self.qual.with_len(seq.len()).make_arbitrary(u)?
        } else {
            self.qual.make_arbitrary(u)?
        };

        let pos = match (self.nonzero_pos, self.cap_end_pos) {
            (false, false) => usize::arbitrary(u)?,
            (false, true) => {
                if let Some(match_len) = cigar.ref_len_in_alignment_checked() {
                    u.int_in_range(0..=(usize::MAX - match_len))?
                } else {
                    usize::arbitrary(u)?
                }
            }
            (true, false) => NonZeroUsize::arbitrary(u)?.into(),
            (true, true) => {
                if let Some(match_len) = cigar.ref_len_in_alignment_checked()
                    && match_len < usize::MAX
                {
                    u.int_in_range(1..=(usize::MAX - match_len))?
                } else {
                    NonZeroUsize::arbitrary(u)?.get()
                }
            }
        };

        let mut data = SamData::new(
            self.qname.make_arbitrary(u)?,
            u16::arbitrary(u)?,
            self.rname.make_arbitrary(u)?,
            pos,
            u8::arbitrary(u)?,
            cigar,
            seq,
            qual,
        );
        data.opt_fields = self.opt_fields.make_arbitrary(u)?;
        Ok(data)
    }
}

impl<'a> Arbitrary<'a> for SamOptRaw {
    /// Generates arbitrary raw optional SAM fields from the given unstructured
    /// data.
    ///
    /// Generated fields may have malformed tags or values, and tags may be
    /// duplicated.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        Ok(Vec::<String>::arbitrary(u)?.into_iter().collect())
    }
}

const SAM_OPT_TAG_FIRST_CHARS: &[u8; 52] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
const SAM_OPT_TAG_SECOND_CHARS: &[u8; 62] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

/// Specifications for generating arbitrary [`SamOptRaw`] optional fields.
///
/// Unlike [`Arbitrary`] for [`SamOptRaw`], these specifications construct each
/// field from a tag and a [`SamOptValue`]. Consequently, generated fields use
/// one of the supported optional-value type codes and are bounded in number by
/// [`SamOptRawSpecs::max_len`]. The defaults still allow arbitrary tags and
/// values, as well as duplicate tags, so they do not guarantee valid SAM
/// syntax.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SamOptRawSpecs {
    /// The maximum number of optional fields to generate.
    pub max_len: usize,

    /// The specifications for generating each optional field.
    pub field: SamOptFieldSpecs,

    /// Whether to ensure that generated tags are unique.
    ///
    /// When enabled, generated tags are retried until they do not duplicate an
    /// earlier tag.
    pub unique_tags: bool,
}

impl SamOptRawSpecs {
    /// Specifications for generating optional fields with valid, unique tags
    /// and BAM-compatible numeric and array-length bounds. Generated tags are
    /// valid and exclude `CG`.
    #[inline]
    #[must_use]
    pub fn bam_encodable(max_len: usize, max_value_len: usize) -> Self {
        Self {
            max_len,
            field: SamOptFieldSpecs {
                valid_tag: true,
                value:     SamOptValueSpecs {
                    valid_sam:     true,
                    bam_encodable: true,
                    max_len:       max_value_len,
                },
            },
            unique_tags: true,
        }
    }
}

impl Default for SamOptRawSpecs {
    fn default() -> Self {
        Self {
            max_len:     usize::MAX,
            field:       SamOptFieldSpecs::default(),
            unique_tags: false,
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for SamOptRawSpecs {
    type Output = SamOptRaw;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        let mut out = SamOptRaw::new();
        let mut tags = Vec::new();

        for _ in 0..self.max_len {
            if !bool::arbitrary(u).unwrap_or(false) {
                break;
            }

            let tag = loop {
                let tag = if self.field.valid_tag {
                    [*u.choose(SAM_OPT_TAG_FIRST_CHARS)?, *u.choose(SAM_OPT_TAG_SECOND_CHARS)?]
                } else {
                    Arbitrary::arbitrary(u)?
                };

                if (self.field.valid_tag && tag == *b"CG") || (self.unique_tags && tags.contains(&tag)) {
                    if u.is_empty() {
                        return Ok(out);
                    }
                    continue;
                }

                break tag;
            };

            if self.unique_tags {
                tags.push(tag);
            }

            let tag = String::from_utf8_lossy(&tag);
            let value = self.field.value.make_arbitrary(u)?;
            out.push(&tag, &value);
        }

        Ok(out)
    }
}

impl<'a> Arbitrary<'a> for SamOptField {
    /// Generates an arbitrary [`SamOptField`] from the given unstructured data.
    ///
    /// The tag and value are generated independently and are not required to
    /// satisfy SAM syntax.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        Ok(SamOptField {
            tag:   Arbitrary::arbitrary(u)?,
            value: SamOptValue::arbitrary(u)?,
        })
    }
}

/// Specifications for generating arbitrary parsed optional SAM fields.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct SamOptFieldSpecs {
    /// Whether to generate valid, non-`CG` tags.
    ///
    /// When false, tags are arbitrary two-byte values and may include `CG`.
    pub valid_tag: bool,

    /// The specifications for generating the optional field value.
    pub value: SamOptValueSpecs,
}

impl<'a> ArbitrarySpecs<'a> for SamOptFieldSpecs {
    type Output = SamOptField;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        let tag = if self.valid_tag {
            loop {
                let tag = [*u.choose(SAM_OPT_TAG_FIRST_CHARS)?, *u.choose(SAM_OPT_TAG_SECOND_CHARS)?];

                if tag != *b"CG" {
                    break tag;
                }
            }
        } else {
            Arbitrary::arbitrary(u)?
        };

        Ok(SamOptField {
            tag,
            value: self.value.make_arbitrary(u)?,
        })
    }
}

impl<'a> Arbitrary<'a> for SamOptValue {
    /// Generates an arbitrary [`SamOptValue`] from the given unstructured data.
    ///
    /// The generated value is not required to display as valid SAM text.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        match u.int_in_range(0_u8..=5)? {
            0 => Ok(SamOptValue::Char(u8::arbitrary(u)?)),
            1 => Ok(SamOptValue::Int(i64::arbitrary(u)?)),
            2 => Ok(SamOptValue::Float(f32::arbitrary(u)?)),
            3 => Ok(SamOptValue::String(String::arbitrary(u)?)),
            4 => Ok(SamOptValue::Hex(String::arbitrary(u)?)),
            5 => Ok(SamOptValue::Array(OptArray::arbitrary(u)?)),
            _ => unreachable!(),
        }
    }
}

/// Specifications for generating arbitrary optional SAM values.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct SamOptValueSpecs {
    /// Whether to restrict value characters and numeric values to
    /// SAM-compatible ranges. Hex values are limited to hexadecimal digits.
    pub valid_sam: bool,

    /// Whether to apply BAM-compatible restrictions to integer values and array
    /// lengths.
    ///
    /// Enabling this also enables `valid_sam`.
    pub bam_encodable: bool,

    /// The maximum text length for generated string and hex values, and the
    /// maximum element count for generated array values.
    pub max_len: usize,
}

impl Default for SamOptValueSpecs {
    fn default() -> Self {
        Self {
            valid_sam:     false,
            bam_encodable: false,
            max_len:       usize::MAX,
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for SamOptValueSpecs {
    type Output = SamOptValue;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        const HEX_VAL: &[u8; 16] = b"0123456789ABCDEF";
        let valid_sam = self.valid_sam || self.bam_encodable;

        match u.int_in_range(0_u8..=5)? {
            0 => {
                let val = if valid_sam {
                    u.int_in_range(b'!'..=b'~')?
                } else {
                    u8::arbitrary(u)?
                };

                Ok(SamOptValue::Char(val))
            }
            1 => {
                let val = if self.bam_encodable {
                    u.int_in_range(i64::from(i32::MIN)..=i64::from(u32::MAX))?
                } else {
                    i64::arbitrary(u)?
                };
                Ok(SamOptValue::Int(val))
            }
            2 => {
                let mut val = f32::arbitrary(u)?;
                if valid_sam && !val.is_finite() {
                    val = 0.0;
                }
                Ok(SamOptValue::Float(val))
            }
            3 => {
                let bytes = VecSpecs {
                    element_specs: ByteSet::AsciiGraphicOrSpace,
                    max_len: self.max_len,
                    ..Default::default()
                }
                .make_arbitrary(u)?;

                Ok(SamOptValue::String(
                    String::from_utf8(bytes).expect("generated bytes are ASCII"),
                ))
            }
            4 => {
                let hex = if self.bam_encodable {
                    let bytes = VecSpecs {
                        element_specs: ByteSet::Any,
                        max_len: self.max_len / 2,
                        ..Default::default()
                    }
                    .make_arbitrary(u)?;

                    let bytes = bytes
                        .into_iter()
                        .flat_map(|byte| [HEX_VAL[usize::from(byte >> 4)], HEX_VAL[usize::from(byte & 0x0f)]])
                        .collect();

                    String::from_utf8(bytes).expect("generated bytes are ASCII")
                } else {
                    StringSpecs {
                        set: ByteSet::Custom(HEX_VAL),
                        max_len: self.max_len,
                        ..Default::default()
                    }
                    .make_arbitrary(u)?
                };

                Ok(SamOptValue::Hex(hex))
            }
            5 => Ok(SamOptValue::Array(
                OptArraySpecs {
                    valid_sam,
                    bam_encodable: self.bam_encodable,
                    max_len: self.max_len,
                }
                .make_arbitrary(u)?,
            )),
            _ => unreachable!(),
        }
    }
}

impl<'a> Arbitrary<'a> for OptArray {
    /// Generates an arbitrary [`OptArray`] from the given unstructured data.
    ///
    /// The generated array subtype and elements are not required to satisfy
    /// SAM or BAM constraints.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> arbitrary::Result<Self> {
        match u.int_in_range(0_u8..=6)? {
            0 => Ok(OptArray::I8(Vec::arbitrary(u)?)),
            1 => Ok(OptArray::U8(Vec::arbitrary(u)?)),
            2 => Ok(OptArray::I16(Vec::arbitrary(u)?)),
            3 => Ok(OptArray::U16(Vec::arbitrary(u)?)),
            4 => Ok(OptArray::I32(Vec::arbitrary(u)?)),
            5 => Ok(OptArray::U32(Vec::arbitrary(u)?)),
            6 => Ok(OptArray::F32(Vec::arbitrary(u)?)),
            _ => unreachable!(),
        }
    }
}

/// Specifications for generating arbitrary `B` array optional values.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct OptArraySpecs {
    /// Whether to generate arrays that display as valid SAM text.
    pub valid_sam: bool,

    /// Whether to restrict the array element count to BAM's `u32` count field.
    pub bam_encodable: bool,

    /// The maximum number of array elements to generate. This is reduced to no
    /// more than `u32::MAX` if `bam_encodable` is enabled.
    pub max_len: usize,
}

impl Default for OptArraySpecs {
    fn default() -> Self {
        Self {
            valid_sam:     false,
            bam_encodable: false,
            max_len:       usize::MAX,
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for OptArraySpecs {
    type Output = OptArray;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> arbitrary::Result<Self::Output> {
        let max_len = if self.bam_encodable {
            self.max_len.min(usize::try_from(u32::MAX).unwrap_or(usize::MAX))
        } else {
            self.max_len
        };

        match u.int_in_range(0_u8..=6)? {
            0 => Ok(OptArray::I8(
                u.arbitrary_iter::<i8>()?.take(max_len).collect::<arbitrary::Result<_>>()?,
            )),
            1 => Ok(OptArray::U8(
                u.arbitrary_iter::<u8>()?.take(max_len).collect::<arbitrary::Result<_>>()?,
            )),
            2 => Ok(OptArray::I16(
                u.arbitrary_iter::<i16>()?.take(max_len).collect::<arbitrary::Result<_>>()?,
            )),
            3 => Ok(OptArray::U16(
                u.arbitrary_iter::<u16>()?.take(max_len).collect::<arbitrary::Result<_>>()?,
            )),
            4 => Ok(OptArray::I32(
                u.arbitrary_iter::<i32>()?.take(max_len).collect::<arbitrary::Result<_>>()?,
            )),
            5 => Ok(OptArray::U32(
                u.arbitrary_iter::<u32>()?.take(max_len).collect::<arbitrary::Result<_>>()?,
            )),
            6 => {
                let float_specs = FloatSpecs::<f32> {
                    include_nan: !self.valid_sam,
                    include_infinite: !self.valid_sam,
                    ..Default::default()
                };

                let vals = float_specs
                    .make_arbitrary_iter(u)
                    .take(max_len)
                    .collect::<arbitrary::Result<Vec<_>>>()?;

                Ok(OptArray::F32(vals))
            }
            _ => unreachable!(),
        }
    }
}
