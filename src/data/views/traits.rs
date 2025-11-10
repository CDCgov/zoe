use crate::{data::views::SliceRange, private::Sealed};

/// A trait for types which have a length. The benefit of including this in a
/// trait is that it can be used in a trait bound, such as [`Slice`].
pub trait Len {
    /// Return whether the data type is empty or not.
    fn is_empty(&self) -> bool;

    /// Get the length of the data.
    fn len(&self) -> usize;
}

/// A trait for associating a type with its corresponding views and owned type.
///
/// The trait bounds on the associated types require that each type implements
/// the corresponding trait ([`DataOwned`], [`DataView`], or [`DataViewMut`]).
/// Furthermore, each must also implement [`ViewAssocTypes`] with the same
/// associated types.
pub trait ViewAssocTypes: Sealed {
    /// The associated owned type.
    type Owned: DataOwned
        + for<'a> ViewAssocTypes<Owned = Self::Owned, View<'a> = Self::View<'a>, ViewMut<'a> = Self::ViewMut<'a>>;

    /// The associated view type.
    type View<'a>: DataView<'a>
        + for<'b> ViewAssocTypes<Owned = Self::Owned, View<'b> = Self::View<'b>, ViewMut<'b> = Self::ViewMut<'b>>;

    /// The associated mutable view type.
    type ViewMut<'a>: DataViewMut<'a>
        + for<'b> ViewAssocTypes<Owned = Self::Owned, View<'b> = Self::View<'b>, ViewMut<'b> = Self::ViewMut<'b>>;
}

/// A trait for data which is owned, and from which views and mutable views can
/// be created.
pub trait DataOwned: ViewAssocTypes<Owned = Self> {
    /// Creates an immutable view of the data.
    #[must_use]
    fn as_view(&self) -> Self::View<'_>;

    /// Creates a mutable view of the data.
    #[must_use]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_>;
}

/// A trait for data which is an immutable view, and from which owned data can
/// be created (via cloning).
pub trait DataView<'a>: ViewAssocTypes {
    /// Creates an owned copy of the data via cloning.
    #[must_use]
    fn to_owned_data(&self) -> Self::Owned;

    /// Reborrows the view, reducing the lifetime to `'b`.
    #[must_use]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b;
}

/// A trait for data which is a mutable view, and from which an immutable view
/// or owned data can be created (the latter requiring cloning).
///
/// ## Parameters
///
/// The lifetime parameter is the lifetime of the stored data inside the view,
/// required to facilitate [`to_view`].
///
/// [`to_view`]: DataViewMut::to_view
pub trait DataViewMut<'a>: ViewAssocTypes {
    /// Creates an immutable view of the data.
    fn as_view(&self) -> Self::View<'_>;

    /// Creates an immutable view by consuming the mutable one.
    fn to_view(self) -> Self::View<'a>;

    /// Creates an owned copy of the data via cloning.
    fn to_owned_data(&self) -> Self::Owned;

    /// Reborrows the mutable view, reducing the lifetime to `'b`.
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b;
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
pub trait Slice: ViewAssocTypes + Len {
    /// Create an immutable view of the data contained in `range`.
    fn slice<R: SliceRange>(&self, range: R) -> Self::View<'_>;
    /// Create an immutable view of the data contained in `range`, or return
    /// `None` if it is out of bounds.
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<Self::View<'_>>;
}

/// Provides the ability to obtain a mutable view of a range of the data.
pub trait SliceMut: Slice {
    /// Create a mutable view of the data contained in `range`.
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> Self::ViewMut<'_>;

    /// Create a mutable view of the data contained in `range`, or return `None`
    /// if it is out of bounds.
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<Self::ViewMut<'_>>;
}

impl Len for Vec<u8> {
    #[inline]
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.len()
    }
}

impl Len for &[u8] {
    #[inline]
    fn is_empty(&self) -> bool {
        (*self).is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        (*self).len()
    }
}

impl Len for &mut [u8] {
    #[inline]
    fn is_empty(&self) -> bool {
        self.as_ref().is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.as_ref().len()
    }
}

impl ViewAssocTypes for Vec<u8> {
    type Owned = Vec<u8>;
    type View<'a> = &'a [u8];
    type ViewMut<'a> = &'a mut [u8];
}

impl ViewAssocTypes for &[u8] {
    type Owned = Vec<u8>;
    type View<'a> = &'a [u8];
    type ViewMut<'a> = &'a mut [u8];
}

impl ViewAssocTypes for &mut [u8] {
    type Owned = Vec<u8>;
    type View<'a> = &'a [u8];
    type ViewMut<'a> = &'a mut [u8];
}

impl DataOwned for Vec<u8> {
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        self
    }

    #[inline]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_> {
        self
    }
}

impl<'a> DataView<'a> for &'a [u8] {
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        self.to_vec()
    }

    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        self
    }
}

impl<'a> DataViewMut<'a> for &'a mut [u8] {
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        self
    }

    #[inline]
    fn to_view(self) -> Self::View<'a> {
        self
    }

    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        self.to_vec()
    }

    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        self
    }
}

impl Restrict for &[u8] {
    #[inline]
    fn restrict<R: SliceRange>(&mut self, range: R) {
        *self = &self[range];
    }

    #[inline]
    fn clear(&mut self) {
        *self = &[];
    }
}

impl Restrict for &mut [u8] {
    #[inline]
    fn restrict<R: SliceRange>(&mut self, range: R) {
        let data = std::mem::take(self);
        *self = &mut data[range];
    }

    #[inline]
    fn clear(&mut self) {
        *self = &mut [];
    }
}

impl Slice for Vec<u8> {
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> Self::View<'_> {
        &self[range]
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<Self::View<'_>> {
        self.get(range)
    }
}

impl Slice for &[u8] {
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> Self::View<'_> {
        &self[range]
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<Self::View<'_>> {
        self.get(range)
    }
}

impl Slice for &mut [u8] {
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> Self::View<'_> {
        &self[range]
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<Self::View<'_>> {
        self.get(range)
    }
}

impl SliceMut for Vec<u8> {
    #[inline]
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> &mut [u8] {
        &mut self[range]
    }

    #[inline]
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<&mut [u8]> {
        self.get_mut(range)
    }
}

impl SliceMut for &mut [u8] {
    #[inline]
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> &mut [u8] {
        &mut self[range]
    }

    #[inline]
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<&mut [u8]> {
        self.get_mut(range)
    }
}
