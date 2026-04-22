use crate::{
    data::views::{
        AsView, AsViewMut, Restrict, SliceRange, ToOwnedData, ToView, impl_len_for_views_generic,
        impl_view_assoc_types_generic,
    },
    prelude::*,
};

impl_len_for_views_generic!(FastQ, FastQView, FastQViewMut, sequence);
impl_view_assoc_types_generic!(FastQ, FastQView, FastQViewMut);

impl ToOwnedData for FastQView<'_> {
    #[inline]
    fn to_owned_data(&self) -> FastQ {
        FastQ {
            header:   self.header.to_string(),
            sequence: self.sequence.to_owned_data(),
            quality:  self.quality.to_owned_data(),
        }
    }
}

impl ToOwnedData for FastQViewMut<'_> {
    #[inline]
    fn to_owned_data(&self) -> FastQ {
        FastQ {
            header:   (*self.header).clone(),
            sequence: self.sequence.to_owned_data(),
            quality:  self.quality.to_owned_data(),
        }
    }
}

impl AsView for FastQ {
    #[inline]
    fn as_view(&self) -> FastQView<'_> {
        FastQView {
            header:   &self.header,
            sequence: self.sequence.as_view(),
            quality:  self.quality.as_view(),
        }
    }
}

impl AsView for FastQViewMut<'_> {
    #[inline]
    fn as_view(&self) -> FastQView<'_> {
        FastQView {
            header:   self.header,
            sequence: self.sequence.as_view(),
            quality:  self.quality.as_view(),
        }
    }
}

impl AsViewMut for FastQ {
    #[inline]
    fn as_view_mut(&mut self) -> FastQViewMut<'_> {
        FastQViewMut {
            header:   &mut self.header,
            sequence: self.sequence.as_view_mut(),
            quality:  self.quality.as_view_mut(),
        }
    }
}

impl<'a> ToView<'a> for FastQViewMut<'a> {
    #[inline]
    fn to_view(self) -> FastQView<'a> {
        FastQView {
            header:   self.header,
            sequence: self.sequence.to_view(),
            quality:  self.quality.to_view(),
        }
    }
}

impl<'a> DataView<'a> for FastQView<'a> {
    #[inline]
    fn reborrow_view<'b>(&'b self) -> FastQView<'b>
    where
        'a: 'b, {
        FastQView {
            header:   self.header,
            sequence: self.sequence.reborrow_view(),
            quality:  self.quality.reborrow_view(),
        }
    }
}

impl<'a> DataViewMut<'a> for FastQViewMut<'a> {
    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        FastQViewMut {
            header:   self.header,
            sequence: self.sequence.reborrow_view_mut(),
            quality:  self.quality.reborrow_view_mut(),
        }
    }
}

impl Restrict for FastQView<'_> {
    #[inline]
    fn restrict<R: SliceRange>(&mut self, range: R) {
        self.sequence.restrict(range.clone());
        self.quality.restrict(range);
    }

    #[inline]
    fn clear(&mut self) {
        self.sequence.clear();
        self.quality.clear();
    }
}

impl Restrict for FastQViewMut<'_> {
    #[inline]
    fn restrict<R: SliceRange>(&mut self, range: R) {
        self.sequence.restrict(range.clone());
        self.quality.restrict(range);
    }

    #[inline]
    fn clear(&mut self) {
        self.sequence.clear();
        self.quality.clear();
    }
}

impl Slice for FastQ {
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastQView<'_> {
        FastQView {
            header:   &self.header,
            sequence: self.sequence.slice(range.clone()),
            quality:  self.quality.slice(range),
        }
    }

    // TODO: Maybe optimize to only perform one get or check
    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastQView<'_>> {
        Some(FastQView {
            header:   &self.header,
            sequence: self.sequence.get_slice(range.clone())?,
            quality:  self.quality.get_slice(range)?,
        })
    }
}

impl SliceMut for FastQ {
    #[inline]
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> FastQViewMut<'_> {
        FastQViewMut {
            header:   &mut self.header,
            sequence: self.sequence.slice_mut(range.clone()),
            quality:  self.quality.slice_mut(range),
        }
    }

    #[inline]
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<FastQViewMut<'_>> {
        Some(FastQViewMut {
            header:   &mut self.header,
            sequence: self.sequence.get_slice_mut(range.clone())?,
            quality:  self.quality.get_slice_mut(range)?,
        })
    }
}

impl Slice for FastQView<'_> {
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastQView<'_> {
        FastQView {
            header:   self.header,
            sequence: self.sequence.slice(range.clone()),
            quality:  self.quality.slice(range),
        }
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastQView<'_>> {
        Some(FastQView {
            header:   self.header,
            sequence: self.sequence.get_slice(range.clone())?,
            quality:  self.quality.get_slice(range)?,
        })
    }
}

impl Slice for FastQViewMut<'_> {
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastQView<'_> {
        FastQView {
            header:   self.header,
            sequence: self.sequence.slice(range.clone()),
            quality:  self.quality.slice(range),
        }
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastQView<'_>> {
        Some(FastQView {
            header:   self.header,
            sequence: self.sequence.get_slice(range.clone())?,
            quality:  self.quality.get_slice(range)?,
        })
    }
}

impl SliceMut for FastQViewMut<'_> {
    #[inline]
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> FastQViewMut<'_> {
        FastQViewMut {
            header:   self.header,
            sequence: self.sequence.slice_mut(range.clone()),
            quality:  self.quality.slice_mut(range),
        }
    }

    #[inline]
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<FastQViewMut<'_>> {
        Some(FastQViewMut {
            header:   self.header,
            sequence: self.sequence.get_slice_mut(range.clone())?,
            quality:  self.quality.get_slice_mut(range)?,
        })
    }
}
