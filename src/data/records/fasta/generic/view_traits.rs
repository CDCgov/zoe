use crate::{
    data::{
        fasta::generic::{Fasta, FastaView, FastaViewMut},
        views::{AsView, AsViewMut, AssocOwnedType, AssocViewMutType, AssocViewType, SliceRange, ToOwnedData, ToView},
    },
    prelude::{DataView, DataViewMut, Len, Restrict, Slice, SliceMut},
};

impl<S: Len> Len for Fasta<S> {
    #[inline]
    fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.sequence.len()
    }
}

impl<'a, S> Len for FastaView<'a, S>
where
    S: AssocViewType<View<'a>: Len>,
{
    #[inline]
    fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.sequence.len()
    }
}

impl<'a, S> Len for FastaViewMut<'a, S>
where
    S: AssocViewMutType<ViewMut<'a>: Len>,
{
    #[inline]
    fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.sequence.len()
    }
}

impl<S> AssocOwnedType for Fasta<S> {
    type Owned = Fasta<S>;
}

impl<S> AssocOwnedType for FastaView<'_, S>
where
    S: AssocViewType,
{
    type Owned = Fasta<S>;
}

impl<S> AssocOwnedType for FastaViewMut<'_, S>
where
    S: AssocViewMutType,
{
    type Owned = Fasta<S>;
}

impl<S> AssocViewType for Fasta<S>
where
    S: AssocViewType,
{
    type View<'a> = FastaView<'a, S>;
}

impl<S> AssocViewType for FastaView<'_, S>
where
    S: AssocViewType,
{
    type View<'a> = FastaView<'a, S>;
}

impl<S> AssocViewType for FastaViewMut<'_, S>
where
    S: AssocViewType + AssocViewMutType,
{
    type View<'a> = FastaView<'a, S>;
}

impl<S> AssocViewMutType for Fasta<S>
where
    S: AssocViewMutType,
{
    type ViewMut<'a> = FastaViewMut<'a, S>;
}

impl<S> AssocViewMutType for FastaView<'_, S>
where
    S: AssocViewType + AssocViewMutType,
{
    type ViewMut<'a> = FastaViewMut<'a, S>;
}

impl<S> AssocViewMutType for FastaViewMut<'_, S>
where
    S: AssocViewMutType,
{
    type ViewMut<'a> = FastaViewMut<'a, S>;
}

impl<'a, S> ToOwnedData for FastaView<'a, S>
where
    S: AssocViewType<View<'a>: ToOwnedData<Owned = S>>,
{
    fn to_owned_data(&self) -> Self::Owned {
        Fasta {
            header:   self.header.to_string(),
            sequence: self.sequence.to_owned_data(),
        }
    }
}

impl<'a, S> ToOwnedData for FastaViewMut<'a, S>
where
    S: AssocViewMutType<ViewMut<'a>: ToOwnedData<Owned = S>>,
{
    #[inline]
    fn to_owned_data(&self) -> Fasta<S> {
        Fasta {
            header:   (*self.header).clone(),
            sequence: self.sequence.to_owned_data(),
        }
    }
}

impl<S> AsView for Fasta<S>
where
    S: AssocViewType + AsView,
{
    #[inline]
    fn as_view(&self) -> FastaView<'_, S> {
        FastaView {
            header:   &self.header,
            sequence: self.sequence.as_view(),
        }
    }
}

impl<'a, S> AsView for FastaViewMut<'a, S>
where
    S: AssocViewType + for<'b> AssocViewMutType<ViewMut<'a>: AsView<View<'b> = S::View<'b>>>,
{
    #[inline]
    fn as_view(&self) -> FastaView<'_, S> {
        FastaView {
            header:   self.header,
            sequence: self.sequence.as_view(),
        }
    }
}

impl<S> AsViewMut for Fasta<S>
where
    S: AsViewMut,
{
    #[inline]
    fn as_view_mut(&mut self) -> FastaViewMut<'_, S> {
        FastaViewMut {
            header:   &mut self.header,
            sequence: self.sequence.as_view_mut(),
        }
    }
}

impl<'a, S> ToView<'a> for FastaViewMut<'a, S>
where
    S: AssocViewType + AssocViewMutType<ViewMut<'a>: ToView<'a, View<'a> = S::View<'a>>>,
{
    #[inline]
    fn to_view(self) -> FastaView<'a, S> {
        FastaView {
            header:   self.header,
            sequence: self.sequence.to_view(),
        }
    }
}

impl<'a, S> DataView<'a> for FastaView<'a, S>
where
    S: AssocViewType,
{
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        FastaView {
            header:   self.header,
            sequence: self.sequence.reborrow_view(),
        }
    }
}

impl<'a, S> DataViewMut<'a> for FastaViewMut<'a, S>
where
    S: AssocViewMutType,
{
    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        FastaViewMut {
            header:   self.header,
            sequence: self.sequence.reborrow_view_mut(),
        }
    }
}

impl<'a, S> Restrict for FastaView<'a, S>
where
    S: AssocViewType<View<'a>: Restrict>,
{
    #[inline]
    fn restrict<R: SliceRange>(&mut self, range: R) {
        self.sequence.restrict(range);
    }

    #[inline]
    fn clear(&mut self) {
        self.sequence.clear();
    }
}

impl<'a, S> Restrict for FastaViewMut<'a, S>
where
    S: AssocViewMutType<ViewMut<'a>: Restrict>,
{
    #[inline]
    fn restrict<R: SliceRange>(&mut self, range: R) {
        self.sequence.restrict(range);
    }

    #[inline]
    fn clear(&mut self) {
        self.sequence.clear();
    }
}

impl<S> Slice for Fasta<S>
where
    S: AssocViewType + Slice,
{
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastaView<'_, S> {
        FastaView {
            header:   &self.header,
            sequence: self.sequence.slice(range),
        }
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastaView<'_, S>> {
        Some(FastaView {
            header:   &self.header,
            sequence: self.sequence.get_slice(range)?,
        })
    }
}

impl<S> SliceMut for Fasta<S>
where
    S: AssocViewMutType + SliceMut,
{
    #[inline]
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> FastaViewMut<'_, S> {
        FastaViewMut {
            header:   &mut self.header,
            sequence: self.sequence.slice_mut(range),
        }
    }

    #[inline]
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<FastaViewMut<'_, S>> {
        Some(FastaViewMut {
            header:   &mut self.header,
            sequence: self.sequence.get_slice_mut(range)?,
        })
    }
}

impl<'a, S> Slice for FastaView<'a, S>
where
    S: AssocViewType<View<'a>: Slice> + 'a,
{
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastaView<'_, S> {
        FastaView {
            header:   self.header,
            sequence: self.sequence.slice(range),
        }
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastaView<'_, S>> {
        Some(FastaView {
            header:   self.header,
            sequence: self.sequence.get_slice(range)?,
        })
    }
}

impl<'a, S> Slice for FastaViewMut<'a, S>
where
    S: AssocViewType + AssocViewMutType<ViewMut<'a>: for<'b> Slice<View<'b> = S::View<'b>>> + 'a,
{
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastaView<'_, S> {
        FastaView {
            header:   self.header,
            sequence: self.sequence.slice(range),
        }
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastaView<'_, S>> {
        Some(FastaView {
            header:   self.header,
            sequence: self.sequence.get_slice(range)?,
        })
    }
}

impl<'a, S> SliceMut for FastaViewMut<'a, S>
where
    S: AssocViewMutType<ViewMut<'a>: SliceMut> + 'a,
{
    #[inline]
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> FastaViewMut<'_, S> {
        FastaViewMut {
            header:   self.header,
            sequence: self.sequence.slice_mut(range),
        }
    }

    #[inline]
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<FastaViewMut<'_, S>> {
        Some(FastaViewMut {
            header:   self.header,
            sequence: self.sequence.get_slice_mut(range)?,
        })
    }
}
