use crate::alignment::AlignmentStates;

#[cfg(test)]
mod test;

#[cfg(test)]
mod benches;

mod error;
mod iter;
mod views;

pub use error::*;
pub use iter::*;
pub use views::*;

/// A [CIGAR string] of increment-operation pairs used in sequence alignment.
///
/// [`Cigar`] internally stores a buffer of the bytes used when displaying a
/// CIGAR string. This makes it efficient for displaying, as well as reading
/// (since no parsing is done). Some algorithms in *Zoe* assume that a CIGAR
/// string meets certain assumptions, in which case these are documented.
/// Otherwise, a CIGAR string may contain arbitrary data, such as through
/// [`from_vec_unchecked`].
///
/// For an alternative data type which stores the increment-operation pairs
/// directly (and hence has more data quality guarantees), see
/// [`AlignmentStates`].
///
/// [CIGAR string]:
///     https://en.wikipedia.org/wiki/Sequence_alignment#CIGAR_Format
/// [`from_vec_unchecked`]: Cigar::from_vec_unchecked
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Cigar(pub(crate) Vec<u8>);

/// The number of bytes in a [`usize`]. This varies based on the target
/// architecture.
#[cfg(target_pointer_width = "16")]
const USIZE_WIDTH: usize = 5;
/// The number of bytes in a [`usize`]. This varies based on the target
/// architecture.
#[cfg(target_pointer_width = "32")]
const USIZE_WIDTH: usize = 10;
/// The number of bytes in a [`usize`]. This varies based on the target
/// architecture.
#[cfg(target_pointer_width = "64")]
const USIZE_WIDTH: usize = 20;

impl Cigar {
    /// Retrieves the underlying byte string from the CIGAR string.
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_slice()
    }

    /// Creates a CIGAR string from a Vec of bytes without checking for
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_vec_unchecked(v: Vec<u8>) -> Self {
        Cigar(v)
    }

    /// Creates a CIGAR string from a slice of bytes without checking for
    /// validity.
    #[inline]
    #[must_use]
    pub fn from_slice_unchecked<T: AsRef<[u8]>>(v: T) -> Self {
        Cigar(v.as_ref().to_vec())
    }

    /// Creates a CIGAR string from an iterator of ciglets without checking for
    /// validity.
    pub fn from_ciglets_unchecked<I>(ciglets: I) -> Self
    where
        I: IntoIterator<Item = Ciglet>, {
        let mut cigar = Vec::new();
        let mut format_buffer = itoa::Buffer::new();

        for ciglet in ciglets {
            let Ciglet { inc, op } = ciglet;
            cigar.extend_from_slice(format_buffer.format(inc).as_bytes());
            cigar.push(op);
        }

        Cigar(cigar)
    }

    /// Creates a new empty CIGAR string
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Cigar(Vec::new())
    }

    /// Expands a CIGAR to an opcode-only byte string, e.g. "3M" → "MMM".
    #[must_use]
    pub(crate) fn expand_cigar(&self) -> ExpandedCigar {
        let mut expanded = Vec::new();

        for Ciglet { inc, op } in self {
            expanded.extend(std::iter::repeat_n(op, inc));
        }

        ExpandedCigar(expanded)
    }

    /// Returns an iterator over the contained [`Ciglet`] values.
    #[inline]
    #[must_use]
    pub fn iter(&self) -> CigletIterator<'_> {
        CigletIterator::new(&self.0)
    }

    /// Returns an iterator over the opcodes in the CIGAR string, with each
    /// repeated the number of times specified in the CIGAR string.
    #[inline]
    pub fn expanded_cigar_iter(&self) -> impl Iterator<Item = u8> {
        self.iter().flat_map(|Ciglet { inc, op }| std::iter::repeat_n(op, inc))
    }

    /// Validates that a CIGAR only contains valid increment-operation pairs.
    ///
    /// Specifically, the assumptions that are checked are:
    ///
    /// - The CIGAR operations must be among `M, I, D, N, S, H, P, X, =`
    /// - Every operation in the CIGAR string must have a preceding increment
    /// - Every increment must be followed by an operation
    /// - The increment for each operation must be non-zero and less than
    ///   [`usize::MAX`]
    ///
    /// This does not confirm that adjacent operations are distinct.
    #[inline]
    #[must_use]
    pub fn is_valid(&self) -> bool {
        Self::check_for_err(&self.0).is_ok()
    }

    /// Checks for errors in the CIGAR string, returning the particular
    /// [`CigarError`] if one is present.
    ///
    /// See the documentation for [`CigletIteratorChecked`] for a list of the
    /// assumptions that are checked.
    fn check_for_err(bytes: &[u8]) -> Result<(), CigarError> {
        if bytes == b"*" || bytes.is_empty() {
            return Ok(());
        }
        for ciglet in CigletIteratorChecked::new(bytes) {
            ciglet?;
        }
        Ok(())
    }
}

impl Default for Cigar {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Display for Cigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.0.is_empty() || self.0 == b"*" {
            write!(f, "*")
        } else if self.0.is_ascii() {
            // SAFETY: we just checked it is ASCII and ASCII is valid UTF8.
            f.write_str(unsafe { std::str::from_utf8_unchecked(&self.0) })
        } else {
            f.write_str(&String::from_utf8_lossy(&self.0))
        }
    }
}

impl std::fmt::Debug for Cigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Cigar({})", String::from_utf8_lossy(&self.0))
    }
}

impl From<ExpandedCigar> for Cigar {
    #[inline]
    fn from(c: ExpandedCigar) -> Cigar {
        c.condense_to_cigar()
    }
}

impl TryFrom<String> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(s: String) -> Result<Self, Self::Error> {
        Cigar::try_from(s.into_bytes())
    }
}

impl TryFrom<Vec<u8>> for Cigar {
    type Error = CigarError;

    fn try_from(bytes: Vec<u8>) -> Result<Self, Self::Error> {
        Self::check_for_err(&bytes)?;
        Ok(Cigar(bytes))
    }
}

impl<const N: usize> TryFrom<[u8; N]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: [u8; N]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_vec())
    }
}

impl TryFrom<&[u8]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &[u8]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_owned())
    }
}

impl TryFrom<&mut [u8]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_vec())
    }
}

impl<const N: usize> TryFrom<&[u8; N]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &[u8; N]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_vec())
    }
}

impl<const N: usize> TryFrom<&mut [u8; N]> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8; N]) -> Result<Self, Self::Error> {
        Cigar::try_from(v.to_vec())
    }
}

impl TryFrom<&str> for Cigar {
    type Error = CigarError;

    #[inline]
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        Cigar::try_from(s.as_bytes())
    }
}

/// A single increment-opcode pair.
///
/// This is yielded when iterating over a [`Cigar`] string. One can also collect
/// an iterator of [`Ciglet`]s into a `Result<Cigar, CigarError>`.
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Ciglet {
    /// Increment or repetition count.
    pub inc: usize,
    /// Alignment operation code.
    pub op:  u8,
}

impl PartialEq<AlignmentStates> for Cigar {
    #[inline]
    fn eq(&self, other: &AlignmentStates) -> bool {
        other.eq(self)
    }
}

impl std::fmt::Display for Ciglet {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}{}", self.inc, self.op as char)
    }
}

impl std::fmt::Debug for Ciglet {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({},{})", self.inc, self.op as char)
    }
}

impl FromIterator<Ciglet> for Result<Cigar, CigarError> {
    /// Converts an iterator of [`Ciglet`]s into a [`Cigar`] string, and checks
    /// that the result is valid.
    #[inline]
    fn from_iter<T: IntoIterator<Item = Ciglet>>(data: T) -> Self {
        let mut byte_string = Vec::new();
        let mut format_buffer = itoa::Buffer::new();

        for Ciglet { inc, op } in data {
            if inc == 0 {
                return Err(CigarError::IncZero);
            } else if !is_valid_op(op) {
                return Err(CigarError::InvalidOperation);
            }

            byte_string.extend_from_slice(format_buffer.format(inc).as_bytes());
            byte_string.push(op);
        }

        Ok(Cigar(byte_string))
    }
}

/// A CIGAR string as an "expanded" opcode-only byte string, e.g. `3M` → `MMM`.
#[derive(Clone, PartialEq)]
pub(crate) struct ExpandedCigar(Vec<u8>);

impl ExpandedCigar {
    /// Condenses the [`ExpandedCigar`] back to its standard form,  e.g. `MMM` →
    /// `3M`.
    #[inline]
    #[must_use]
    pub(crate) fn condense_to_cigar(self) -> Cigar {
        let mut condensed: Vec<u8> = Vec::new();
        let mut format_buffer = itoa::Buffer::new();

        let mut cigars = self.0.iter().copied().filter(|op| is_valid_op(*op));

        let Some(mut previous) = cigars.next() else {
            return Cigar(condensed);
        };

        let mut count: usize = 1;

        for op in cigars {
            if previous == op {
                count += 1;
            } else {
                condensed.extend_from_slice(format_buffer.format(count).as_bytes());
                condensed.push(previous);
                previous = op;
                count = 1;
            }
        }

        condensed.extend_from_slice(format_buffer.format(count).as_bytes());
        condensed.push(previous);

        Cigar(condensed)
    }
}

impl std::fmt::Display for ExpandedCigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.0))
    }
}

impl std::fmt::Debug for ExpandedCigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.0))
    }
}

impl From<Cigar> for ExpandedCigar {
    #[inline]
    fn from(c: Cigar) -> ExpandedCigar {
        c.expand_cigar()
    }
}

impl From<Vec<u8>> for ExpandedCigar {
    #[inline]
    fn from(vec: Vec<u8>) -> Self {
        ExpandedCigar(vec)
    }
}

impl From<&str> for ExpandedCigar {
    #[inline]
    fn from(s: &str) -> Self {
        ExpandedCigar(s.as_bytes().to_owned())
    }
}

impl From<&[u8]> for ExpandedCigar {
    #[inline]
    fn from(v: &[u8]) -> Self {
        ExpandedCigar(v.to_vec())
    }
}

impl<const N: usize> From<&[u8; N]> for ExpandedCigar {
    #[inline]
    fn from(v: &[u8; N]) -> Self {
        ExpandedCigar(v.to_vec())
    }
}

/// Returns whether a byte represents a valid CIGAR operation.
///
/// The following operations are valid: `MIDNSHPX=`.
#[inline]
pub(crate) const fn is_valid_op(op: u8) -> bool {
    matches!(op, b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'=')
}

/// A trait allowing the length of the query or reference in an alignment to be
/// re-computed for CIGAR-like data.
///
/// This is implemented for anything that can be converted into an iterator of
/// [`Ciglet`] values via the trait [`ToCigletIterator`].
pub trait LenInAlignment {
    /// Sums lengths for operations `M`, `D`, `N`, `X`, and `=`, which is the
    /// length of the reference that is included in the alignment.
    #[must_use]
    fn ref_len_in_alignment(&self) -> usize;

    /// Sums lengths for operations `M`, `I`, `S`, `X`, and `=`, which is the
    /// length of the query that is included in the alignment.
    #[must_use]
    fn query_len_in_alignment(&self) -> usize;

    /// A checked version of [`ref_len_in_alignment`], returning `None` on
    /// overflow.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *fuzzing* feature in your `Cargo.toml` to use this
    /// method.
    ///
    /// </div>
    ///
    /// [`ref_len_in_alignment`]: LenInAlignment::ref_len_in_alignment
    #[must_use]
    #[cfg(feature = "fuzzing")]
    fn ref_len_in_alignment_checked(&self) -> Option<usize>;

    /// A checked version of [`query_len_in_alignment`], returning `None` on
    /// overflow.
    ///
    /// <div class="warning note">
    ///
    /// **Note**
    ///
    /// You must enable the *fuzzing* feature in your `Cargo.toml` to use this
    /// method.
    ///
    /// </div>
    ///
    /// [`query_len_in_alignment`]: LenInAlignment::query_len_in_alignment
    #[must_use]
    #[cfg(feature = "fuzzing")]
    fn query_len_in_alignment_checked(&self) -> Option<usize>;
}

impl<T> LenInAlignment for T
where
    T: ToCigletIterator,
{
    #[inline]
    fn ref_len_in_alignment(&self) -> usize {
        self.to_ciglet_iterator()
            .filter_map(|Ciglet { inc, op }| matches!(op, b'M' | b'D' | b'N' | b'=' | b'X').then_some(inc))
            .sum()
    }

    #[inline]
    fn query_len_in_alignment(&self) -> usize {
        self.to_ciglet_iterator()
            .filter_map(|Ciglet { inc, op }| matches!(op, b'M' | b'I' | b'S' | b'=' | b'X').then_some(inc))
            .sum()
    }

    #[inline]
    #[cfg(feature = "fuzzing")]
    fn ref_len_in_alignment_checked(&self) -> Option<usize> {
        self.to_ciglet_iterator()
            .filter_map(|Ciglet { inc, op }| matches!(op, b'M' | b'D' | b'N' | b'=' | b'X').then_some(inc))
            .try_fold(0, usize::checked_add)
    }

    #[inline]
    #[cfg(feature = "fuzzing")]
    fn query_len_in_alignment_checked(&self) -> Option<usize> {
        self.to_ciglet_iterator()
            .filter_map(|Ciglet { inc, op }| matches!(op, b'M' | b'I' | b'S' | b'=' | b'X').then_some(inc))
            .try_fold(0, usize::checked_add)
    }
}

/// A trait enabling inspection of a CIGAR string (or similar struct) for
/// certain conditions.
pub trait InspectCigar {
    /// Returns whether the alignment contains any insertions (`I`) or deletions
    /// (`D`).
    #[must_use]
    fn has_indels(&self) -> bool;
}

impl<T> InspectCigar for T
where
    T: ToCigletIterator,
{
    #[inline]
    fn has_indels(&self) -> bool {
        self.to_ciglet_iterator().any(|c| matches!(c.op, b'I' | b'D'))
    }
}
