/// A trait for types which have a length. The benefit of including this in a
/// trait is that it can be used in a trait bound, such as [`Slice`].
pub trait Len {
    /// Return whether the data type is empty or not.
    fn is_empty(&self) -> bool;

    /// Get the length of the data.
    fn len(&self) -> usize;
}

/// A trait for data which is owned, and from which views and mutable views can
/// be created.
pub trait DataOwned: Sealed {
    type View<'a>
    where
        Self: 'a;

    type ViewMut<'a>
    where
        Self: 'a;

    /// Create an immutable view of the data.
    fn as_view(&self) -> Self::View<'_>;

    /// Create a mutable view of the data.
    fn as_view_mut(&mut self) -> Self::ViewMut<'_>;
}

/// A trait for data which is an immutable view, and from which owned data can
/// be created (via cloning).
pub trait DataView: Sealed {
    type Owned;

    /// Create an owned copy of the data via cloning.
    fn to_owned_data(&self) -> Self::Owned;
}

/// A trait for data which is a mutable view, and from which an immutable view
/// or owned data can be created (the latter requiring cloning).
pub trait DataViewMut<'b>: Sealed {
    type View<'a>
    where
        Self: 'a;
    type Owned;

    /// Create an immutable view of the data.
    fn as_view(&self) -> Self::View<'_>;

    /// Creates an immutable view by consuming the mutable one.
    fn to_view(self) -> Self::View<'b>;

    /// Create an owned copy of the data via cloning.
    fn to_owned_data(&self) -> Self::Owned;
}

/// Provides the ability to resize a view in-place.
pub trait Restrict {
    /// Re-slice the view in-place, changing the view to hold only a subslice of
    /// its original data.
    fn restrict<R: SliceRange>(&mut self, range: R);

    /// Clear the view so that it is empty (this does not affect the underlying
    /// data).
    fn clear(&mut self);
}

/// Provides the ability to obtain an immutable view of a range of the data.
pub trait Slice: Len + Sealed {
    type View<'a>
    where
        Self: 'a;

    /// Create an immutable view of the data contained in `range`.
    fn slice<R: SliceRange>(&self, range: R) -> Self::View<'_>;
    /// Create an immutable view of the data contained in `range`, or return
    /// `None` if it is out of bounds.
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<Self::View<'_>>;
}

/// Provides the ability to obtain a mutable view of a range of the data.
pub trait SliceMut: Slice {
    type ViewMut<'a>
    where
        Self: 'a;

    /// Create a mutable view of the data contained in `range`.
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> Self::ViewMut<'_>;
    /// Create a mutable view of the data contained in `range`, or return `None`
    /// if it is out of bounds.
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<Self::ViewMut<'_>>;
}

use crate::{data::views::SliceRange, private::Sealed};
