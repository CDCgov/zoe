use crate::{data::views::SliceRange, private::Sealed};

/// A trait for types which have a length. The benefit of including this in a
/// trait is that it can be used in a trait bound or in generic code.
pub trait Len {
    /// Return whether the data type is empty or not.
    fn is_empty(&self) -> bool;

    /// Get the length of the data.
    fn len(&self) -> usize;
}

/// Associates an owned data type with `Self`. This can be used as a trait
/// bound, such as for [`ToOwnedData`].
///
/// The view API expects this to also be implemented for owned types, with the
/// associated type equal to `Self`.
pub trait AssocOwnedType: Sealed {
    /// The owned data type associated with `Self`.
    type Owned: AssocOwnedType<Owned = Self::Owned>;
}

/// Associates a view data type with `Self`. This can be used as a trait bound,
/// such as for [`AsView`], [`ToView`], [`DataView`], or [`Slice`].
///
/// The view API expects this to also be implemented for view types, with the
/// associated type equal to `Self` but with the generic lifetime.
pub trait AssocViewType: Sealed {
    /// The view data type associated with `Self`.
    type View<'a>: DataView<'a> + for<'b> AssocViewType<View<'b> = Self::View<'b>>;
}

/// Associates a mutable view data type with `Self`. This can be used as a trait
/// bound, such as for [`AsViewMut`], [`DataViewMut`], or [`SliceMut`].
///
/// The view API expects this to also be implemented for mutable view types,
/// with the associated type equal to `Self` but with the generic lifetime.
pub trait AssocViewMutType: Sealed {
    /// The mutable view data type associated with `Self`.
    type ViewMut<'a>: DataViewMut<'a> + for<'b> AssocViewMutType<ViewMut<'b> = Self::ViewMut<'b>>;
}

/// A trait for converting immutable or mutable view types into owned types via
/// cloning.
pub trait ToOwnedData: AssocOwnedType {
    /// Creates an owned copy of the data via cloning.
    #[must_use]
    fn to_owned_data(&self) -> Self::Owned;
}

/// A trait for forming an immutable view from an owned type (or a mutable
/// view).
pub trait AsView: AssocViewType {
    /// Creates an immutable view of the data.
    #[must_use]
    fn as_view(&self) -> Self::View<'_>;
}

/// A trait for forming a mutable view from an owned type.
pub trait AsViewMut: AssocViewMutType {
    /// Creates a mutable view of the data.
    #[must_use]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_>;
}

/// A trait for converting a mutable view into an immutable view.
pub trait ToView<'a>: AssocViewType {
    /// Creates an immutable view by consuming the mutable one.
    fn to_view(self) -> Self::View<'a>;
}

/// A trait for other functions specific to immutable view types.
pub trait DataView<'a>: AssocViewType {
    /// Reborrows the view, reducing the lifetime to `'b`.
    #[must_use]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b;
}

/// A trait for other functions specific to mutable view types.
pub trait DataViewMut<'a>: AssocViewMutType {
    /// Reborrows the mutable view, reducing the lifetime to `'b`.
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b;
}

/// Provides the ability to resize a view in-place.
pub trait Restrict {
    /// Re-slices the view in-place, changing the view to hold only a subslice
    /// of its original data.
    fn restrict<R: SliceRange>(&mut self, range: R);

    /// Clears the view so that it is empty (this does not affect the underlying
    /// data).
    fn clear(&mut self);
}

/// Provides the ability to obtain an immutable view of a range of the data.
pub trait Slice: AssocViewType {
    /// Creates an immutable view of the data contained in `range`.
    fn slice<R: SliceRange>(&self, range: R) -> Self::View<'_>;

    /// Creates an immutable view of the data contained in `range`, or return
    /// `None` if it is out of bounds.
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<Self::View<'_>>;
}

/// Similar to [`Slice`], but designed to be implemented on immutable view types
/// that implement [`Copy`]. The functionality is the same, but the lifetimes of
/// the returned views are not tied to the input view.
pub trait SliceCopy: Copy {
    /// Creates an immutable view of the data contained in `range`.
    #[must_use]
    fn slice<R: SliceRange>(self, range: R) -> Self;

    /// Creates an immutable view of the data contained in `range`, or return
    /// `None` if it is out of bounds.
    #[must_use]
    fn get_slice<R: SliceRange>(self, range: R) -> Option<Self>;
}

/// Provides the ability to obtain a mutable view of a range of the data.
pub trait SliceMut: AssocViewMutType {
    /// Creates a mutable view of the data contained in `range`.
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> Self::ViewMut<'_>;

    /// Creates a mutable view of the data contained in `range`, or return
    /// `None` if it is out of bounds.
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<Self::ViewMut<'_>>;
}

impl Len for Vec<u8> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }
}

impl Len for &[u8] {
    fn is_empty(&self) -> bool {
        (*self).is_empty()
    }

    fn len(&self) -> usize {
        (*self).len()
    }
}

impl Len for &mut [u8] {
    fn is_empty(&self) -> bool {
        (**self).is_empty()
    }

    fn len(&self) -> usize {
        (**self).len()
    }
}

impl AssocOwnedType for Vec<u8> {
    type Owned = Vec<u8>;
}

impl AssocViewType for Vec<u8> {
    type View<'a> = &'a [u8];
}

impl AssocViewMutType for Vec<u8> {
    type ViewMut<'a> = &'a mut [u8];
}

impl AssocOwnedType for &[u8] {
    type Owned = Vec<u8>;
}

impl AssocViewType for &[u8] {
    type View<'a> = &'a [u8];
}

impl AssocViewMutType for &[u8] {
    type ViewMut<'a> = &'a mut [u8];
}

impl AssocOwnedType for &mut [u8] {
    type Owned = Vec<u8>;
}

impl AssocViewType for &mut [u8] {
    type View<'a> = &'a [u8];
}

impl AssocViewMutType for &mut [u8] {
    type ViewMut<'a> = &'a mut [u8];
}

impl AsView for Vec<u8> {
    fn as_view(&self) -> Self::View<'_> {
        self
    }
}

impl AsView for &mut [u8] {
    fn as_view(&self) -> Self::View<'_> {
        self
    }
}

impl AsViewMut for Vec<u8> {
    fn as_view_mut(&mut self) -> Self::ViewMut<'_> {
        self
    }
}

impl ToOwnedData for &[u8] {
    fn to_owned_data(&self) -> Self::Owned {
        self.to_vec()
    }
}

impl ToOwnedData for &mut [u8] {
    fn to_owned_data(&self) -> Self::Owned {
        self.to_vec()
    }
}

impl<'a> ToView<'a> for &'a mut [u8] {
    fn to_view(self) -> Self::View<'a> {
        self
    }
}

impl<'a> DataView<'a> for &'a [u8] {
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        self
    }
}

impl<'a> DataViewMut<'a> for &'a mut [u8] {
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

impl SliceCopy for &[u8] {
    #[inline]
    fn slice<R: SliceRange>(self, range: R) -> Self {
        &self[range]
    }

    #[inline]
    fn get_slice<R: SliceRange>(self, range: R) -> Option<Self> {
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
