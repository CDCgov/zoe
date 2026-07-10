use crate::{
    alignment::{Alignment, MaybeAligned, types::AlignmentStates},
    data::{
        cigar::{CigarView, is_valid_op},
        types::cigar::{Cigar, CigarError, Ciglet, CigletIteratorChecked},
    },
};
use std::{cmp::Ordering, fmt::Write};

impl AsRef<[Ciglet]> for AlignmentStates {
    #[inline]
    fn as_ref(&self) -> &[Ciglet] {
        &self.0
    }
}

impl PartialEq<Cigar> for AlignmentStates {
    #[inline]
    fn eq(&self, other: &Cigar) -> bool {
        let mut o = other.iter();
        let matches = self.iter().copied().eq(o.by_ref());
        matches && o.valid()
    }
}

impl<'a> IntoIterator for &'a AlignmentStates {
    type Item = Ciglet;
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, Ciglet>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter().copied()
    }
}

impl IntoIterator for AlignmentStates {
    type Item = Ciglet;
    type IntoIter = <Vec<Ciglet> as IntoIterator>::IntoIter;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl FromIterator<Ciglet> for Result<AlignmentStates, CigarError> {
    #[inline]
    fn from_iter<T: IntoIterator<Item = Ciglet>>(iter: T) -> Self {
        let mut last_op = 0;
        let mut repeated_op = false;

        // We do not expect repeated operations, so use collect and its
        // specialized implementation relying on size_hint and/or in-place
        // iteration
        let mut out = iter
            .into_iter()
            .map(|ciglet| {
                if ciglet.inc == 0 {
                    Err(CigarError::IncZero)
                } else if !is_valid_op(ciglet.op) {
                    Err(CigarError::InvalidOperation)
                } else if ciglet.op == last_op {
                    repeated_op = true;
                    Ok(ciglet)
                } else {
                    last_op = ciglet.op;
                    Ok(ciglet)
                }
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Check whether a second pass for merging repeated operations is needed
        if repeated_op {
            // The original index in the output being read from
            let mut original_index = 1;
            // The current index being accumulated into
            let mut merged_index = 0;
            while let Some(ciglet) = out.get(original_index).copied() {
                let merged_ciglet = &mut out[merged_index];

                if ciglet.op == merged_ciglet.op {
                    merged_ciglet.inc = merged_ciglet
                        .inc
                        .checked_add(ciglet.inc)
                        .ok_or(CigarError::MergeIncOverflow)?;
                } else {
                    // Advance to the next available position to accumulate into
                    merged_index += 1;
                    out[merged_index] = out[original_index];
                }
                original_index += 1;
            }
            // We have written to merged_index, so we must keep merged_index+1
            // length
            out.truncate(merged_index + 1);
        }

        Ok(AlignmentStates(out))
    }
}

impl std::fmt::Debug for AlignmentStates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.0.iter().map(|c| (c.inc, c.op as char))).finish()
    }
}

impl std::fmt::Display for AlignmentStates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buff = core::fmt::NumBuffer::new();
        for Ciglet { inc, op } in self {
            f.write_str(inc.format_into(&mut buff))?;
            f.write_char(op as char)?;
        }
        Ok(())
    }
}

impl TryFrom<&[u8]> for AlignmentStates {
    type Error = CigarError;

    fn try_from(v: &[u8]) -> Result<Self, CigarError> {
        let mut states = AlignmentStates::new();

        if v == b"*" {
            Ok(states)
        } else {
            let mut iter = CigletIteratorChecked::new(v);
            let Some(prev_ciglet) = iter.next() else {
                return Ok(states);
            };
            let mut prev_ciglet = prev_ciglet?;

            // More efficient version of add_ciglet: we don't need to check
            // whether a previous ciglet exists, nor do we have to check for
            // non-zero increments
            for ciglet in iter {
                let ciglet = ciglet?;
                if prev_ciglet.op == ciglet.op {
                    prev_ciglet.inc = prev_ciglet.inc.checked_add(ciglet.inc).ok_or(CigarError::MergeIncOverflow)?;
                } else {
                    // Validity: CigletIteratorChecked errors on zero
                    // increments, and we have ensured the operation is
                    // different from the previous
                    states.0.push(prev_ciglet);
                    prev_ciglet = ciglet;
                }
            }

            // Push the last ciglet
            states.0.push(prev_ciglet);
            Ok(states)
        }
    }
}

impl TryFrom<String> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(s: String) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(s.as_bytes())
    }
}

impl TryFrom<Vec<u8>> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(bytes: Vec<u8>) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(bytes.as_slice())
    }
}

impl<const N: usize> TryFrom<[u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: [u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl TryFrom<&mut [u8]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(&*v)
    }
}

impl<const N: usize> TryFrom<&[u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &[u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl<const N: usize> TryFrom<&mut [u8; N]> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(v: &mut [u8; N]) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(v.as_slice())
    }
}

impl TryFrom<&str> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(s.as_bytes())
    }
}

impl TryFrom<Cigar> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(value: Cigar) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(value.as_bytes())
    }
}

impl TryFrom<&Cigar> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(value: &Cigar) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(value.as_bytes())
    }
}

impl<'a> TryFrom<CigarView<'a>> for AlignmentStates {
    type Error = CigarError;

    #[inline]
    fn try_from(value: CigarView<'a>) -> Result<Self, Self::Error> {
        AlignmentStates::try_from(value.as_bytes())
    }
}

impl<T> PartialOrd for Alignment<T>
where
    T: PartialOrd,
{
    /// This method returns an ordering between `self` and `other` values if one
    /// exists.
    ///
    /// The ordering is based on the `score` field, and equivalent scores either
    /// produce `None` or [`Ordering::Equal`] (if all other fields are
    /// identical).
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self.score.partial_cmp(&other.score) {
            Some(Ordering::Equal) => {
                if self == other {
                    Some(Ordering::Equal)
                } else {
                    None
                }
            }
            ord => ord,
        }
    }
}

impl<T> PartialOrd for MaybeAligned<T>
where
    T: PartialOrd,
{
    /// This method returns an ordering between `self` and `other` values if one
    /// exists.
    ///
    /// The [`Unmapped`] variant is considered the smallest, while the
    /// [`Overflowed`] variant is considered the largest. Comparisons involving
    /// two [`MaybeAligned::Some`] values defer to the implementation of
    /// [`PartialOrd`] on `T`.
    ///
    /// [`Unmapped`]: MaybeAligned::Unmapped
    /// [`Overflowed`]: MaybeAligned::Overflowed
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self, other) {
            // Equality of state-less variants
            (MaybeAligned::Overflowed, MaybeAligned::Overflowed) | (MaybeAligned::Unmapped, MaybeAligned::Unmapped) => {
                Some(Ordering::Equal)
            }

            // Comparing state-less variants
            (MaybeAligned::Overflowed, _) | (_, MaybeAligned::Unmapped) => Some(Ordering::Greater),
            (_, MaybeAligned::Overflowed) | (MaybeAligned::Unmapped, _) => Some(Ordering::Less),

            // Comparing two alignments
            (MaybeAligned::Some(left), MaybeAligned::Some(right)) => left.partial_cmp(right),
        }
    }
}
