use std::slice::SliceIndex;

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
pub trait DataOwned {
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
pub trait DataView {
    type Owned;

    /// Create an owned copy of the data via cloning.
    fn to_owned_data(&self) -> Self::Owned;
    // TODO: Bike shed on naming (inspiration from other crates)
    /// Re-slice the view in-place.
    fn restrict<R: SliceRange>(&mut self, range: R);
}

/// A trait for data which is a mutable view, and from which an immutable view
/// or owned data can be created (the latter requiring cloning).
pub trait DataViewMut {
    type View<'a>
    where
        Self: 'a;
    type Owned;

    /// Create an immutable view of the data.
    fn as_view(&self) -> Self::View<'_>;
    /// Create an owned copy of the data via cloning.
    fn to_owned_data(&self) -> Self::Owned;
    /// Re-slice the view in-place.
    fn restrict<R: SliceRange>(&mut self, range: R);
}

// TODO: Maybe remove, since it seems there isn't a way to make this pub(crate)
/// Alias for the types that can be used to index into an array and get an array
/// out (i.e., ranges)
pub trait SliceRange: SliceIndex<[u8], Output = [u8]> + Clone {}
impl<T: SliceIndex<[u8], Output = [u8]> + Clone> SliceRange for T {}

/// Provides the ability to obtain an immutable view of a range of the data.
pub trait Slice: Len {
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
    fn get_mut_slice<R: SliceRange>(&mut self, range: R) -> Option<Self::ViewMut<'_>>;
}

/// A macro for implementing [`Len`], given the owning type, the view type, the
/// mutable view type, and the name of the field which determines the length.
macro_rules! impl_len {
    ($owned:ident, $view:ident, $viewmut:ident, $lenfield:tt) => {
        impl Len for $owned {
            #[inline]
            fn is_empty(&self) -> bool {
                self.$lenfield.is_empty()
            }

            #[inline]
            fn len(&self) -> usize {
                self.$lenfield.len()
            }
        }

        impl<'a> Len for $view<'a> {
            #[inline]
            fn is_empty(&self) -> bool {
                self.$lenfield.is_empty()
            }

            #[inline]
            fn len(&self) -> usize {
                self.$lenfield.len()
            }
        }

        impl<'a> Len for $viewmut<'a> {
            #[inline]
            fn is_empty(&self) -> bool {
                self.$lenfield.is_empty()
            }

            #[inline]
            fn len(&self) -> usize {
                self.$lenfield.len()
            }
        }
    };
}

/// A macro for implementing the conversions between the owned types and views,
/// given the owning type, the view type, and the mutable view type. Also
/// implements [`Slice`] and [`SliceMut`]. It is assumed
/// that the types are all wrappers around a single field.
macro_rules! impl_views_for_wrapper {
    ($owned:ident, $view:ident, $viewmut:ident) => {
        impl DataOwned for $owned {
            type View<'a> = $view<'a>;
            type ViewMut<'a> = $viewmut<'a>;

            #[inline]
            fn as_view(&self) -> $view<'_> {
                $view(&self.0)
            }

            #[inline]
            fn as_view_mut(&mut self) -> $viewmut<'_> {
                $viewmut(&mut self.0)
            }
        }

        impl<'a> DataView for $view<'a> {
            type Owned = $owned;

            #[inline]
            fn to_owned_data(&self) -> $owned {
                $owned(self.0.into())
            }

            #[inline]
            fn restrict<R: SliceRange>(&mut self, range: R) {
                self.0 = &self.0[range]
            }
        }

        impl DataViewMut for $viewmut<'_> {
            type View<'a>
                = $view<'a>
            where
                Self: 'a;

            type Owned = $owned;

            #[inline]
            fn as_view(&self) -> $view<'_> {
                $view(&self.0)
            }

            #[inline]
            fn to_owned_data(&self) -> $owned {
                $owned(self.0.to_vec())
            }

            #[inline]
            fn restrict<R: SliceRange>(&mut self, range: R) {
                let data = std::mem::take(&mut self.0);
                self.0 = &mut data[range];
            }
        }

        impl Slice for $owned {
            type View<'a> = $view<'a>;

            #[inline]
            fn slice<R: SliceRange>(&self, range: R) -> $view {
                $view(&self.0[range])
            }

            #[inline]
            fn get_slice<R: SliceRange>(&self, range: R) -> Option<$view> {
                Some($view(self.0.get(range)?))
            }
        }

        impl Slice for $view<'_> {
            type View<'a>
                = $view<'a>
            where
                Self: 'a;

            #[inline]
            fn slice<R: SliceRange>(&self, range: R) -> $view {
                $view(&self.0[range])
            }

            #[inline]
            fn get_slice<R: SliceRange>(&self, range: R) -> Option<$view> {
                Some($view(self.0.get(range)?))
            }
        }

        impl Slice for $viewmut<'_> {
            type View<'a>
                = $view<'a>
            where
                Self: 'a;

            #[inline]
            fn slice<R: SliceRange>(&self, range: R) -> $view {
                $view(&self.0[range])
            }

            #[inline]
            fn get_slice<R: SliceRange>(&self, range: R) -> Option<$view> {
                Some($view(self.0.get(range)?))
            }
        }

        impl SliceMut for $owned {
            type ViewMut<'a> = $viewmut<'a>;

            #[inline]
            fn slice_mut<R: SliceRange>(&mut self, range: R) -> $viewmut {
                $viewmut(&mut self.0[range])
            }

            #[inline]
            fn get_mut_slice<R: SliceRange>(&mut self, range: R) -> Option<$viewmut> {
                Some($viewmut(self.0.get_mut(range)?))
            }
        }

        impl SliceMut for $viewmut<'_> {
            type ViewMut<'a>
                = $viewmut<'a>
            where
                Self: 'a;

            #[inline]
            fn slice_mut<R: SliceRange>(&mut self, range: R) -> $viewmut {
                $viewmut(&mut self.0[range])
            }

            #[inline]
            fn get_mut_slice<R: SliceRange>(&mut self, range: R) -> Option<$viewmut> {
                Some($viewmut(self.0.get_mut(range)?))
            }
        }
    };
}

pub(crate) use impl_len;
pub(crate) use impl_views_for_wrapper;
