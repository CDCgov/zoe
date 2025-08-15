use crate::data::{
    cigar::{Cigar, CigarError},
    views::impl_views_for_wrapper,
};

/// A view of a [`Cigar`] string.
///
/// See [Views](crate::data#views) for more details. This struct is still in
/// development; it currently only supports parsing and displaying, and not all
/// the other operations that [`Cigar`] supports.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct CigarView<'a>(&'a [u8]);

/// A mutable view of a [`Cigar`] string.
///
/// See [Views](crate::data#views) for more details. This struct is still in
/// development; it currently only supports parsing and displaying, and not all
/// the other operations that [`Cigar`] supports.
pub struct CigarViewMut<'a>(&'a mut [u8]);

impl<'a> CigarView<'a> {
    /// Creates a new [`CigarView`] from `v` without checking for validity.
    #[inline]
    #[must_use]
    pub fn from_slice_unchecked(v: &'a [u8]) -> Self {
        Self(v)
    }

    /// Creates a new empty [`CigarView`].
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self(b"")
    }
}

impl<'a> CigarViewMut<'a> {
    /// Creates a new [`CigarViewMut`] from `v` without checking for validity.
    #[inline]
    #[must_use]
    pub fn from_slice_unchecked(v: &'a mut [u8]) -> Self {
        Self(v)
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
        Ok(CigarView(bytes))
    }
}

impl<'a> TryFrom<&'a mut [u8]> for CigarView<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a mut [u8]) -> Result<Self, Self::Error> {
        Cigar::check_for_err(bytes)?;
        Ok(CigarView(bytes))
    }
}

impl<'a> TryFrom<&'a mut [u8]> for CigarViewMut<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a mut [u8]) -> Result<Self, Self::Error> {
        Cigar::check_for_err(bytes)?;
        Ok(CigarViewMut(bytes))
    }
}

impl<'a, const N: usize> TryFrom<&'a [u8; N]> for CigarView<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a [u8; N]) -> Result<Self, Self::Error> {
        Cigar::check_for_err(bytes)?;
        Ok(CigarView(bytes))
    }
}

impl<'a, const N: usize> TryFrom<&'a mut [u8; N]> for CigarView<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a mut [u8; N]) -> Result<Self, Self::Error> {
        Cigar::check_for_err(bytes)?;
        Ok(CigarView(bytes))
    }
}

impl<'a, const N: usize> TryFrom<&'a mut [u8; N]> for CigarViewMut<'a> {
    type Error = CigarError;

    fn try_from(bytes: &'a mut [u8; N]) -> Result<Self, Self::Error> {
        Cigar::check_for_err(bytes)?;
        Ok(CigarViewMut(bytes))
    }
}

impl_views_for_wrapper!(Cigar, CigarView, CigarViewMut);
