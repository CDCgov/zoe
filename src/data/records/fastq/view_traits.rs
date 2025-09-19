use crate::{
    data::views::{Restrict, SliceRange, impl_len_for_views},
    prelude::*,
};

impl_len_for_views!(FastQ, FastQView, FastQViewMut, sequence);

impl DataOwned for FastQ {
    type View<'a> = FastQView<'a>;
    type ViewMut<'a> = FastQViewMut<'a>;

    #[inline]
    fn as_view(&self) -> FastQView<'_> {
        FastQView {
            header:   &self.header,
            sequence: self.sequence.as_view(),
            quality:  self.quality.as_view(),
        }
    }

    #[inline]
    fn as_view_mut(&mut self) -> FastQViewMut<'_> {
        FastQViewMut {
            header:   &mut self.header,
            sequence: self.sequence.as_view_mut(),
            quality:  self.quality.as_view_mut(),
        }
    }
}

impl DataView for FastQView<'_> {
    type Owned = FastQ;

    #[inline]
    fn to_owned_data(&self) -> FastQ {
        FastQ {
            header:   self.header.to_string(),
            sequence: self.sequence.to_owned_data(),
            quality:  self.quality.to_owned_data(),
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

impl<'b> DataViewMut<'b> for FastQViewMut<'b> {
    type View<'a>
        = FastQView<'a>
    where
        Self: 'a;

    type Owned = FastQ;

    #[inline]
    fn as_view(&self) -> FastQView<'_> {
        FastQView {
            header:   self.header,
            sequence: self.sequence.as_view(),
            quality:  self.quality.as_view(),
        }
    }

    #[inline]
    fn to_view(self) -> FastQView<'b> {
        FastQView {
            header:   self.header,
            sequence: self.sequence.to_view(),
            quality:  self.quality.to_view(),
        }
    }

    #[inline]
    fn to_owned_data(&self) -> FastQ {
        FastQ {
            header:   (*self.header).clone(),
            sequence: self.sequence.to_owned_data(),
            quality:  self.quality.to_owned_data(),
        }
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
    type View<'a> = FastQView<'a>;

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
    type ViewMut<'a> = FastQViewMut<'a>;

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
    type View<'a>
        = FastQView<'a>
    where
        Self: 'a;

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
    type View<'a>
        = FastQView<'a>
    where
        Self: 'a;

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
    type ViewMut<'a>
        = FastQViewMut<'a>
    where
        Self: 'a;

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
