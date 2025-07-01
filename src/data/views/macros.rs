/// A macro for implementing [`Len`], given the owning type, the view type, the
/// mutable view type, and the name of the field which determines the length.
///
/// [`Len`]: crate::data::views::Len
macro_rules! impl_len {
    ($owned:ident, $view:ident, $viewmut:ident, $lenfield:tt) => {
        impl $crate::data::views::Len for $owned {
            #[inline]
            fn is_empty(&self) -> bool {
                self.$lenfield.is_empty()
            }

            #[inline]
            fn len(&self) -> usize {
                self.$lenfield.len()
            }
        }

        impl<'a> $crate::data::views::Len for $view<'a> {
            #[inline]
            fn is_empty(&self) -> bool {
                self.$lenfield.is_empty()
            }

            #[inline]
            fn len(&self) -> usize {
                self.$lenfield.len()
            }
        }

        impl<'a> $crate::data::views::Len for $viewmut<'a> {
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
/// implements [`Slice`] and [`SliceMut`]. It is assumed that the types are all
/// wrappers around a single field.
///
/// [`Slice`]: crate::data::views::traits::Slice
/// [`SliceMut`]: crate::data::views::traits::SliceMut
macro_rules! impl_views_for_wrapper {
    ($owned:ident, $view:ident, $viewmut:ident) => {
        impl $crate::data::views::DataOwned for $owned {
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

        impl<'a> $crate::data::views::DataView for $view<'a> {
            type Owned = $owned;

            #[inline]
            fn to_owned_data(&self) -> $owned {
                $owned(self.0.into())
            }
        }

        impl<'b> $crate::data::views::DataViewMut<'b> for $viewmut<'b> {
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
            fn to_view(self) -> $view<'b> {
                $view(self.0)
            }

            #[inline]
            fn to_owned_data(&self) -> $owned {
                $owned(self.0.to_vec())
            }
        }

        impl $crate::data::views::Slice for $owned {
            type View<'a> = $view<'a>;

            #[inline]
            fn slice<R: $crate::data::views::SliceRange>(&self, range: R) -> $view<'_> {
                $view(&self.0[range])
            }

            #[inline]
            fn get_slice<R: $crate::data::views::SliceRange>(&self, range: R) -> Option<$view<'_>> {
                Some($view(self.0.get(range)?))
            }
        }

        impl $crate::data::views::Slice for $view<'_> {
            type View<'a>
                = $view<'a>
            where
                Self: 'a;

            #[inline]
            fn slice<R: $crate::data::views::SliceRange>(&self, range: R) -> $view<'_> {
                $view(&self.0[range])
            }

            #[inline]
            fn get_slice<R: $crate::data::views::SliceRange>(&self, range: R) -> Option<$view<'_>> {
                Some($view(self.0.get(range)?))
            }
        }

        impl $crate::data::views::Slice for $viewmut<'_> {
            type View<'a>
                = $view<'a>
            where
                Self: 'a;

            #[inline]
            fn slice<R: $crate::data::views::SliceRange>(&self, range: R) -> $view<'_> {
                $view(&self.0[range])
            }

            #[inline]
            fn get_slice<R: $crate::data::views::SliceRange>(&self, range: R) -> Option<$view<'_>> {
                Some($view(self.0.get(range)?))
            }
        }

        impl $crate::data::views::SliceMut for $owned {
            type ViewMut<'a> = $viewmut<'a>;

            #[inline]
            fn slice_mut<R: $crate::data::views::SliceRange>(&mut self, range: R) -> $viewmut<'_> {
                $viewmut(&mut self.0[range])
            }

            #[inline]
            fn get_slice_mut<R: $crate::data::views::SliceRange>(&mut self, range: R) -> Option<$viewmut<'_>> {
                Some($viewmut(self.0.get_mut(range)?))
            }
        }

        impl $crate::data::views::SliceMut for $viewmut<'_> {
            type ViewMut<'a>
                = $viewmut<'a>
            where
                Self: 'a;

            #[inline]
            fn slice_mut<R: $crate::data::views::SliceRange>(&mut self, range: R) -> $viewmut<'_> {
                $viewmut(&mut self.0[range])
            }

            #[inline]
            fn get_slice_mut<R: $crate::data::views::SliceRange>(&mut self, range: R) -> Option<$viewmut<'_>> {
                Some($viewmut(self.0.get_mut(range)?))
            }
        }
    };
}

macro_rules! impl_restrict_for_wrapper {
    ($owned:ident, $view:ident, $viewmut:ident) => {
        impl<'a> $crate::data::views::Restrict for $view<'a> {
            #[inline]
            fn restrict<R: $crate::data::views::SliceRange>(&mut self, range: R) {
                self.0 = &self.0[range];
            }

            #[inline]
            fn clear(&mut self) {
                self.0 = &[];
            }
        }

        impl<'a> $crate::data::views::Restrict for $viewmut<'a> {
            #[inline]
            fn restrict<R: $crate::data::views::SliceRange>(&mut self, range: R) {
                let data = std::mem::take(&mut self.0);
                self.0 = &mut data[range];
            }

            #[inline]
            fn clear(&mut self) {
                self.0 = &mut [];
            }
        }
    };
}

pub(crate) use impl_len;
pub(crate) use impl_restrict_for_wrapper;
pub(crate) use impl_views_for_wrapper;
