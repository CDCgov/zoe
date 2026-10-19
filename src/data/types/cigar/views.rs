//! View types for CIGAR strings.

use crate::data::{
    cigar::{Cigar, CigarError, CigletIterator, CigletIteratorChecked},
    views::{AsView, AssocOwnedType, AssocViewType, DataView, Len, ToOwnedData},
};

/// A view of a [`Cigar`] string.
///
/// See [Views](crate::data#views) for more details. This struct is still in
/// development; it currently only supports parsing and displaying, and not all
/// the other operations that [`Cigar`] supports.
#[derive(Copy, Clone, Eq, PartialEq, Hash)]
pub struct CigarView<'a>(&'a [u8]);

impl<'a> CigarView<'a> {
    /// Retrieves the underlying byte string from the CIGAR string.
    #[inline]
    #[must_use]
    pub fn as_bytes(self) -> &'a [u8] {
        self.0
    }

    /// Creates a new [`CigarView`] from `v` without checking for validity.
    ///
    /// If `*` is passed, this is replaced with [`CigarView::new`].
    #[inline]
    #[must_use]
    pub fn from_slice_unchecked(v: &'a [u8]) -> Self {
        if v == b"*" { CigarView::new() } else { CigarView(v) }
    }

    /// Creates a new empty [`CigarView`].
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self(b"")
    }

    /// Returns an iterator of the contained [`Ciglet`] values.
    ///
    /// [`Ciglet`]: crate::data::types::cigar::Ciglet
    #[inline]
    #[must_use]
    pub fn iter(self) -> CigletIteratorChecked<'a> {
        CigletIteratorChecked::new(self.0)
    }

    /// Returns an iterator of the contained [`Ciglet`] values, truncating on
    /// the first error.
    ///
    /// [`Ciglet`]: crate::data::types::cigar::Ciglet
    #[inline]
    #[must_use]
    pub fn iter_unchecked(self) -> CigletIterator<'a> {
        CigletIterator::new(self.0)
    }
}

impl Default for CigarView<'_> {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Display for CigarView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.0.is_empty() || self.0 == b"*" {
            write!(f, "*")
        } else if self.0.is_ascii() {
            // SAFETY: we just checked it is ASCII and ASCII is valid UTF8.
            f.write_str(unsafe { std::str::from_utf8_unchecked(self.0) })
        } else {
            f.write_str(&String::from_utf8_lossy(self.0))
        }
    }
}

impl std::fmt::Debug for CigarView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CigarView({})", String::from_utf8_lossy(self.0))
    }
}

impl<'a> TryFrom<&'a [u8]> for CigarView<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a [u8]) -> Result<Self, Self::Error> {
        Cigar::check_for_err(bytes)?;
        Ok(CigarView::from_slice_unchecked(bytes))
    }
}

impl<'a> TryFrom<&'a mut [u8]> for CigarView<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a mut [u8]) -> Result<Self, Self::Error> {
        CigarView::try_from(&*bytes)
    }
}

impl<'a, const N: usize> TryFrom<&'a [u8; N]> for CigarView<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a [u8; N]) -> Result<Self, Self::Error> {
        CigarView::try_from(bytes.as_slice())
    }
}

impl<'a, const N: usize> TryFrom<&'a mut [u8; N]> for CigarView<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a mut [u8; N]) -> Result<Self, Self::Error> {
        CigarView::try_from(bytes.as_slice())
    }
}

impl Len for Cigar {
    #[inline]
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl Len for CigarView<'_> {
    #[inline]
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl AssocOwnedType for Cigar {
    type Owned = Cigar;
}

impl AssocViewType for Cigar {
    type View<'a> = CigarView<'a>;
}

impl AssocOwnedType for CigarView<'_> {
    type Owned = Cigar;
}

impl AssocViewType for CigarView<'_> {
    type View<'a> = CigarView<'a>;
}

impl ToOwnedData for CigarView<'_> {
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        Cigar(self.0.into())
    }
}

impl AsView for Cigar {
    fn as_view(&self) -> CigarView<'_> {
        CigarView(&self.0)
    }
}

impl<'a> DataView<'a> for CigarView<'a> {
    fn reborrow_view<'b>(&'b self) -> CigarView<'b>
    where
        'a: 'b, {
        CigarView(self.0)
    }
}
