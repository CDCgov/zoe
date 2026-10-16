use crate::data::cigar::InvalidCigarOpErr;
use std::fmt::Display;

/// An enum of the different alignment states that can exist in a CIGAR string,
/// as used in *Zoe*'s [`Cigar`] and [`Ciglet`].
///
/// The operations match those used by the
/// [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) specifications. The
/// integers used to represent each variant correspond to the ASCII code for the
/// character.
///
/// [`Cigar`]: crate::data::cigar::Cigar
/// [`Ciglet`]: crate::data::cigar::Ciglet
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
#[repr(u8)]
pub enum AlignOp {
    /// Alignment match (can be a sequence match or mismatch).
    ///
    /// This operation consumes positions from both the query and reference. It
    /// is represented by `b'M'`. This is prefered to [`AlignOp::Eq`] and
    /// [`AlignOp::X`] in *Zoe* since it allows more concise and readable
    /// alignment representations.
    M  = b'M',

    /// Insertion to the reference.
    ///
    /// This operation consumes positions from the query only. It is represented
    /// by `b'I'`.
    I  = b'I',

    /// Deletion from the reference.
    ///
    /// This operation consumes positions from the reference only. It is
    /// represented by `b'D'`.
    D  = b'D',

    /// Skipped region from the reference.
    ///
    /// This operation consumes positions from the reference only. It is
    /// represented by `b'N'`.
    N  = b'N',

    /// Soft clipping (clipped sequences present in SAM's `SEQ` field).
    ///
    /// This operation consumes positions from the query only. It is represented
    /// by `b'S'`.
    S  = b'S',

    /// Hard clipping (clipped sequences not present in SAM's `SEQ` field).
    ///
    /// This operation consumes positions from neither the query nor reference.
    /// It is represented by `b'H'`.
    H  = b'H',

    /// Padding (silent deletion from padded reference).
    ///
    /// This operation consumes positions from neither the query nor reference.
    /// It is represented by `b'P'`.
    P  = b'P',

    /// Sequence match.
    ///
    /// This operation consumes positions from both the query and reference. It
    /// is represented by `b'='`.
    Eq = b'=',

    /// Sequence mismatch.
    ///
    /// This operation consumes positions from both the query and reference. It
    /// is represented by `b'X'`.
    X  = b'X',
}

impl AlignOp {
    /// Converts the [`AlignOp`] to a `char` which can be displayed.
    #[must_use]
    pub fn as_char(self) -> char {
        (self as u8) as char
    }
}

impl TryFrom<u8> for AlignOp {
    type Error = InvalidCigarOpErr;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'M' => Ok(AlignOp::M),
            b'I' => Ok(AlignOp::I),
            b'D' => Ok(AlignOp::D),
            b'N' => Ok(AlignOp::N),
            b'S' => Ok(AlignOp::S),
            b'H' => Ok(AlignOp::H),
            b'P' => Ok(AlignOp::P),
            b'=' => Ok(AlignOp::Eq),
            b'X' => Ok(AlignOp::X),
            op => Err(InvalidCigarOpErr(op)),
        }
    }
}

impl Display for AlignOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

impl AlignOp {
    /// Returns true if the [`AlignOp`] consumes positions from the query.
    #[inline]
    #[must_use]
    pub fn consumes_query(self) -> bool {
        matches!(self, AlignOp::M | AlignOp::I | AlignOp::S | AlignOp::Eq | AlignOp::X)
    }

    /// Returns true if the [`AlignOp`] consumes positions from the reference.
    #[inline]
    #[must_use]
    pub fn consumes_ref(self) -> bool {
        matches!(self, AlignOp::M | AlignOp::D | AlignOp::N | AlignOp::Eq | AlignOp::X)
    }
}
