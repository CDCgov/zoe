use crate::{
    alignment::{AlignOp, Alignment, AlignmentStates, MaybeAligned, NextCiglet},
    data::{
        cigar::{CigarError, LenInAlignment},
        err::impl_err_conv,
        nucleotides::Nucleotides,
        phred::QualityScores,
        sam::{Flag, SamData, SamOptRaw, is_missing_sam_field},
        views::Len,
    },
    math::AnyInt,
};
use std::{error::Error, fmt::Display};

impl SamData {
    /// Constructs a new [`SamData`] record from an [`Alignment`] struct as well
    /// as the other provided fields. The score is included as a field under the
    /// tag `AS`.
    ///
    /// For the opposite transformation, see [`SamData::to_alignment`].
    #[inline]
    #[must_use]
    pub fn from_alignment<T: AnyInt + Into<i64>>(
        alignment: &Alignment<T>, qname: String, flag: Flag, rname: String, mapq: u8, seq: Nucleotides, qual: QualityScores,
    ) -> Self {
        // Validate incoming arguments to ensure they are not missing
        let qname = if is_missing_sam_field(&qname) { String::new() } else { qname };
        let rname = if is_missing_sam_field(&rname) { String::new() } else { rname };
        let seq = if is_missing_sam_field(&seq) { Nucleotides::new() } else { seq };
        let qual = if is_missing_sam_field(&qual) {
            QualityScores::new()
        } else {
            qual
        };

        // Both SAM and Alignment exclude clipped bases when reporting
        // positions, so we just need to adjust to 1-based
        let pos = alignment.ref_range.start + 1;
        let cigar = alignment.states.to_cigar_unchecked();
        let opt_fields = SamOptRaw::new_with_score(alignment.score);

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
            opt_fields,
        }
    }

    /// Converts the [`SamData`] record into an [`Alignment`] struct (wrapped in
    /// [`MaybeAligned`]).
    ///
    /// If the entire CIGAR string is empty or consists solely of clipping, or
    /// if [`SamData::flag`] has the unmapped flag set, then
    /// [`MaybeAligned::Unmapped`] is returned.
    ///
    /// Any hard clipped bases in the query are not included in the resulting
    /// `query_range` or `query_len`. `query_range` represents the entire `SEQ`
    /// except for soft clipping, while `query_len` represents the length with
    /// soft clipping (which must also equal the length of `SEQ` if present).
    ///
    /// The `score` must be provided as an argument. A ZST like `()` can also be
    /// passed to avoid storing the score. If the score is present under the
    /// `AS` TAG, then it can be retrieved with:
    ///
    /// ```
    /// # use zoe::{
    /// #     data::{cigar::Cigar, sam::{Flag, SamOptValue, SamData}},
    /// #     prelude::{Nucleotides, QualityScores},
    /// # };
    /// #
    /// # let mut sam_data = SamData::new(
    /// #     String::new(),
    /// #     Flag::default(),
    /// #     String::new(),
    /// #     0,
    /// #     255,
    /// #     Cigar::new(),
    /// #     Nucleotides::new(),
    /// #     QualityScores::new(),
    /// # );
    /// #
    /// # sam_data.opt_fields.push("AS", &SamOptValue::Int(0));
    /// #
    /// let score = sam_data
    ///     .opt_fields
    ///     .get("AS")
    ///     .expect("The optional fields must be formatted properly")
    ///     .expect("The score TAG must be present")
    ///     .int()
    ///     .expect("The score VALUE should be an integer");
    /// ```
    ///
    /// For the opposite transformation, see [`SamData::from_alignment`].
    ///
    /// ## Errors
    ///
    /// Unmapped reads are those where the FLAG has the unmapped bit set, the
    /// CIGAR field is not populated, or the CIGAR field is solely clipping.
    /// These are checked in this order, and [`MaybeAligned::Unmapped`] is
    /// returned immediately when a condition is found true. During this
    /// process, an error is returned if the CIGAR field must be read and is
    /// malformed.
    ///
    /// For any mapped read, the following additional errors could occur:
    ///
    /// - The CIGAR field must have no clipping in the interior of the string
    /// - The number of reference bases consumed by the CIGAR string must not
    ///   overflow [`usize`], and similarly for the query
    /// - If the SEQ field is populated, then the sequence length must match the
    ///   number of query bases consumed by the CIGAR string
    /// - The POS field must be populated (i.e., not zero)
    /// - The POS and CIGAR field must yield an ending reference index that does
    ///   not exceed `ref_len`
    #[inline]
    pub fn to_alignment<T>(&self, score: T, ref_len: usize) -> Result<MaybeAligned<Alignment<T>>, ToAlignmentError> {
        // Handle FLAG
        if self.flag.is_unmapped() {
            return Ok(MaybeAligned::Unmapped);
        }

        // Parse CIGAR to AlignmentStates
        let states = AlignmentStates::try_from(&self.cigar)?;
        if states.is_empty() {
            return Ok(MaybeAligned::Unmapped);
        }

        // Compute the start of the aligned query by examining clipping in CIGAR
        let mut ciglets = states.as_slice();
        ciglets.next_ciglet_if_op(|op| op == AlignOp::H);
        let soft_clipping_front = ciglets.next_ciglet_if_op(|op| op == AlignOp::S).map_or(0, |c| c.inc);
        ciglets.next_ciglet_back_if_op(|op| op == AlignOp::H);
        let soft_clipping_back = ciglets.next_ciglet_back_if_op(|op| op == AlignOp::S).map_or(0, |c| c.inc);
        if ciglets.is_empty() {
            return Ok(MaybeAligned::Unmapped);
        }
        if ciglets.iter().any(|c| matches!(c.op, AlignOp::S | AlignOp::H)) {
            return Err(ToAlignmentError::InteriorClipping);
        }

        let query_range_start = soft_clipping_front;

        // Compute number of residues consumed in CIGAR
        let Some(aligned_ref_len) = states.ref_len_in_alignment_checked() else {
            return Err(ToAlignmentError::OverflowingCigarRefInc);
        };
        let Some(query_len) = states.query_len_in_alignment_checked() else {
            return Err(ToAlignmentError::OverflowingCigarQueryInc);
        };

        // Confirm CIGAR is compatible with SEQ, if SEQ is present
        if !is_missing_sam_field(&self.seq) && self.seq.len() != query_len {
            return Err(ToAlignmentError::IncorrectSeqLen {
                cigar_query_inc: query_len,
                seq_len:         self.seq.len(),
            });
        }

        // Parse POS
        let Some(ref_range_start) = self.pos.checked_sub(1) else {
            return Err(ToAlignmentError::MissingPosField);
        };

        // Compute and validate the index for the end of the reference
        let Some(ref_range_end) = aligned_ref_len.checked_add(ref_range_start).filter(|&end| end <= ref_len) else {
            return Err(ToAlignmentError::IncorrectRefLen { ref_len });
        };

        // Compute the end of the aligned query
        let query_range_end = query_len - soft_clipping_back;

        // Form ranges
        let ref_range = ref_range_start..ref_range_end;
        let query_range = query_range_start..query_range_end;

        Ok(MaybeAligned::Some(Alignment {
            score,
            ref_range,
            query_range,
            states,
            ref_len,
            query_len,
        }))
    }
}

#[derive(Debug)]
pub enum ToAlignmentError {
    InvalidCigar(CigarError),
    InteriorClipping,
    OverflowingCigarRefInc,
    OverflowingCigarQueryInc,
    IncorrectSeqLen { cigar_query_inc: usize, seq_len: usize },
    MissingPosField,
    IncorrectRefLen { ref_len: usize },
}

impl Display for ToAlignmentError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ToAlignmentError::InvalidCigar(_) => f.write_str("The record contained an invalid CIGAR string"),
            ToAlignmentError::InteriorClipping => f.write_str("Found clipping in the interior of the CIGAR string"),
            ToAlignmentError::OverflowingCigarRefInc => {
                f.write_str("The number of reference residues consumed by the CIGAR string overflowed usize")
            }
            ToAlignmentError::OverflowingCigarQueryInc => {
                f.write_str("The number of query residues consumed by the CIGAR string overflowed usize")
            }
            ToAlignmentError::IncorrectSeqLen {
                cigar_query_inc,
                seq_len,
            } => write!(
                f,
                "The CIGAR field consumed {cigar_query_inc} query residues, but SEQ contained {seq_len} residues"
            ),
            ToAlignmentError::MissingPosField => f.write_str("The POS field was not populated"),
            ToAlignmentError::IncorrectRefLen { ref_len } => {
                write!(f, "The alignment implied by POS and CIGAR extends beyond ref_len ({ref_len})")
            }
        }
    }
}

impl Error for ToAlignmentError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ToAlignmentError::InvalidCigar(err) => Some(err),
            _ => None,
        }
    }
}

impl_err_conv!(ToAlignmentError);

impl From<CigarError> for ToAlignmentError {
    fn from(value: CigarError) -> Self {
        Self::InvalidCigar(value)
    }
}
