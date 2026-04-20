use crate::{
    data::views::impl_view_assoc_types,
    prelude::{DataOwned, DataView, DataViewMut},
};

/// A taxon for a record.
#[repr(transparent)]
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct Taxon(pub(super) String);

/// The corresponding immutable view type for [`Taxon`].
///
/// See [Views](crate::data#views) for more details.
#[repr(transparent)]
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct TaxonView<'a>(pub(super) &'a str);

/// The corresponding mutable view type for [`Taxon`].
///
/// See [Views](crate::data#views) for more details.
#[repr(transparent)]
#[derive(Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct TaxonViewMut<'a>(pub(super) &'a mut String);

impl Taxon {
    /// Returns the taxon as a string slice.
    #[inline]
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl<'a> TaxonView<'a> {
    /// Returns the taxon as a string slice.
    #[inline]
    #[must_use]
    pub fn as_str(&self) -> &'a str {
        self.0
    }
}

impl TaxonViewMut<'_> {
    /// Returns the taxon as a string slice.
    #[inline]
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.0
    }
}

impl From<String> for Taxon {
    #[inline]
    fn from(value: String) -> Self {
        Self(value)
    }
}

impl From<&str> for Taxon {
    #[inline]
    fn from(value: &str) -> Self {
        Self(String::from(value))
    }
}

impl<'a> From<TaxonView<'a>> for Taxon {
    #[inline]
    fn from(value: TaxonView<'a>) -> Self {
        Self::from(value.as_str())
    }
}

impl<'a> From<TaxonViewMut<'a>> for Taxon {
    #[inline]
    fn from(value: TaxonViewMut<'a>) -> Self {
        Self::from(value.as_str())
    }
}

impl<'a> From<&'a str> for TaxonView<'a> {
    #[inline]
    fn from(value: &'a str) -> Self {
        Self(value)
    }
}

impl<'a> From<&'a mut String> for TaxonViewMut<'a> {
    #[inline]
    fn from(value: &'a mut String) -> Self {
        Self(value)
    }
}

impl std::fmt::Display for Taxon {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::fmt::Display for TaxonView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::fmt::Display for TaxonViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl_view_assoc_types!(Taxon, TaxonView, TaxonViewMut);

impl DataOwned for Taxon {
    #[inline]
    fn as_view(&self) -> TaxonView<'_> {
        TaxonView(&self.0)
    }

    #[inline]
    fn as_view_mut(&mut self) -> TaxonViewMut<'_> {
        TaxonViewMut(&mut self.0)
    }
}

impl<'a> DataView<'a> for TaxonView<'a> {
    #[inline]
    fn to_owned_data(&self) -> Taxon {
        Taxon(self.0.into())
    }

    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        TaxonView(self.0)
    }
}

impl<'a> DataViewMut<'a> for TaxonViewMut<'a> {
    #[inline]
    fn as_view(&self) -> TaxonView<'_> {
        TaxonView(self.0)
    }

    #[inline]
    fn to_view(self) -> TaxonView<'a> {
        TaxonView(self.0)
    }

    #[inline]
    fn to_owned_data(&self) -> Taxon {
        Taxon(self.0.clone())
    }

    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        TaxonViewMut(self.0)
    }
}
