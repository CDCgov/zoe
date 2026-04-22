use crate::{
    data::{
        fasta::generic::{Fasta, FastaAnnot, FastaAnnotView, FastaAnnotViewMut, FastaView, FastaViewMut},
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

impl<M, S: Len> Len for FastaAnnot<M, S> {
    #[inline]
    fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.sequence.len()
    }
}

impl<'a, M, S> Len for FastaAnnotView<'a, M, S>
where
    M: AssocViewType,
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

impl<'a, M, S> Len for FastaAnnotViewMut<'a, M, S>
where
    M: AssocViewMutType,
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

impl<M, S> AssocOwnedType for FastaAnnot<M, S> {
    type Owned = FastaAnnot<M, S>;
}

impl<M, S> AssocOwnedType for FastaAnnotView<'_, M, S>
where
    M: AssocViewType,
    S: AssocViewType,
{
    type Owned = FastaAnnot<M, S>;
}

impl<M, S> AssocOwnedType for FastaAnnotViewMut<'_, M, S>
where
    M: AssocViewMutType,
    S: AssocViewMutType,
{
    type Owned = FastaAnnot<M, S>;
}

impl<M, S> AssocViewType for FastaAnnot<M, S>
where
    M: AssocViewType,
    S: AssocViewType,
{
    type View<'a> = FastaAnnotView<'a, M, S>;
}

impl<M, S> AssocViewType for FastaAnnotView<'_, M, S>
where
    M: AssocViewType,
    S: AssocViewType,
{
    type View<'a> = FastaAnnotView<'a, M, S>;
}

impl<M, S> AssocViewType for FastaAnnotViewMut<'_, M, S>
where
    M: AssocViewType + AssocViewMutType,
    S: AssocViewType + AssocViewMutType,
{
    type View<'a> = FastaAnnotView<'a, M, S>;
}

impl<M, S> AssocViewMutType for FastaAnnot<M, S>
where
    M: AssocViewMutType,
    S: AssocViewMutType,
{
    type ViewMut<'a> = FastaAnnotViewMut<'a, M, S>;
}

impl<M, S> AssocViewMutType for FastaAnnotView<'_, M, S>
where
    M: AssocViewType + AssocViewMutType,
    S: AssocViewType + AssocViewMutType,
{
    type ViewMut<'a> = FastaAnnotViewMut<'a, M, S>;
}

impl<M, S> AssocViewMutType for FastaAnnotViewMut<'_, M, S>
where
    M: AssocViewMutType,
    S: AssocViewMutType,
{
    type ViewMut<'a> = FastaAnnotViewMut<'a, M, S>;
}

impl<'a, M, S> ToOwnedData for FastaAnnotView<'a, M, S>
where
    M: AssocViewType<View<'a>: ToOwnedData<Owned = M>>,
    S: AssocViewType<View<'a>: ToOwnedData<Owned = S>>,
{
    fn to_owned_data(&self) -> FastaAnnot<M, S> {
        FastaAnnot {
            header:   self.header.to_string(),
            sequence: self.sequence.to_owned_data(),
            annot:    self.annot.to_owned_data(),
        }
    }
}

impl<'a, M, S> ToOwnedData for FastaAnnotViewMut<'a, M, S>
where
    M: AssocViewMutType<ViewMut<'a>: ToOwnedData<Owned = M>>,
    S: AssocViewMutType<ViewMut<'a>: ToOwnedData<Owned = S>>,
{
    #[inline]
    fn to_owned_data(&self) -> FastaAnnot<M, S> {
        FastaAnnot {
            header:   (*self.header).clone(),
            sequence: self.sequence.to_owned_data(),
            annot:    self.annot.to_owned_data(),
        }
    }
}

impl<M, S> AsView for FastaAnnot<M, S>
where
    M: AsView,
    S: AsView,
{
    #[inline]
    fn as_view(&self) -> FastaAnnotView<'_, M, S> {
        FastaAnnotView {
            header:   &self.header,
            sequence: self.sequence.as_view(),
            annot:    self.annot.as_view(),
        }
    }
}

impl<'a, M, S> AsView for FastaAnnotViewMut<'a, M, S>
where
    M: AssocViewType + for<'b> AssocViewMutType<ViewMut<'a>: AsView<View<'b> = M::View<'b>>>,
    S: AssocViewType + for<'b> AssocViewMutType<ViewMut<'a>: AsView<View<'b> = S::View<'b>>>,
{
    #[inline]
    fn as_view(&self) -> FastaAnnotView<'_, M, S> {
        FastaAnnotView {
            header:   self.header,
            sequence: self.sequence.as_view(),
            annot:    self.annot.as_view(),
        }
    }
}

impl<M, S> AsViewMut for FastaAnnot<M, S>
where
    M: AsViewMut,
    S: AsViewMut,
{
    #[inline]
    fn as_view_mut(&mut self) -> FastaAnnotViewMut<'_, M, S> {
        FastaAnnotViewMut {
            header:   &mut self.header,
            sequence: self.sequence.as_view_mut(),
            annot:    self.annot.as_view_mut(),
        }
    }
}

impl<'a, M, S> DataView<'a> for FastaAnnotView<'a, M, S>
where
    M: AssocViewType<View<'a>: DataView<'a>>,
    S: AssocViewType<View<'a>: DataView<'a>>,
{
    #[inline]
    fn reborrow_view<'b>(&'b self) -> FastaAnnotView<'b, M, S>
    where
        'a: 'b, {
        FastaAnnotView {
            header:   self.header,
            sequence: self.sequence.reborrow_view(),
            annot:    self.annot.reborrow_view(),
        }
    }
}

impl<'a, M, S> DataViewMut<'a> for FastaAnnotViewMut<'a, M, S>
where
    M: AssocViewMutType,
    S: AssocViewMutType,
{
    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> FastaAnnotViewMut<'b, M, S>
    where
        'a: 'b, {
        FastaAnnotViewMut {
            header:   self.header,
            sequence: self.sequence.reborrow_view_mut(),
            annot:    self.annot.reborrow_view_mut(),
        }
    }
}

impl<'a, M, S> Restrict for FastaAnnotView<'a, M, S>
where
    M: AssocViewType,
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

impl<'a, M, S> Restrict for FastaAnnotViewMut<'a, M, S>
where
    M: AssocViewMutType,
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

impl<M, S> Slice for FastaAnnot<M, S>
where
    M: AsView,
    S: Slice,
{
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastaAnnotView<'_, M, S> {
        FastaAnnotView {
            header:   &self.header,
            sequence: self.sequence.slice(range),
            annot:    self.annot.as_view(),
        }
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastaAnnotView<'_, M, S>> {
        Some(FastaAnnotView {
            header:   &self.header,
            sequence: self.sequence.get_slice(range)?,
            annot:    self.annot.as_view(),
        })
    }
}

impl<M, S> SliceMut for FastaAnnot<M, S>
where
    M: AsViewMut,
    S: SliceMut,
{
    #[inline]
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> FastaAnnotViewMut<'_, M, S> {
        FastaAnnotViewMut {
            header:   &mut self.header,
            sequence: self.sequence.slice_mut(range),
            annot:    self.annot.as_view_mut(),
        }
    }

    #[inline]
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<FastaAnnotViewMut<'_, M, S>> {
        Some(FastaAnnotViewMut {
            header:   &mut self.header,
            sequence: self.sequence.get_slice_mut(range)?,
            annot:    self.annot.as_view_mut(),
        })
    }
}

impl<'a, M, S> Slice for FastaAnnotView<'a, M, S>
where
    M: AssocViewType,
    S: AssocViewType<View<'a>: Slice>,
{
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastaAnnotView<'_, M, S> {
        FastaAnnotView {
            header:   self.header,
            sequence: self.sequence.slice(range),
            annot:    self.annot.reborrow_view(),
        }
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastaAnnotView<'_, M, S>> {
        Some(FastaAnnotView {
            header:   self.header,
            sequence: self.sequence.get_slice(range)?,
            annot:    self.annot.reborrow_view(),
        })
    }
}

impl<'a, M, S> Slice for FastaAnnotViewMut<'a, M, S>
where
    M: AssocViewType + for<'b> AssocViewMutType<ViewMut<'a>: AsView<View<'b> = M::View<'b>>>,
    S: AssocViewType + for<'b> AssocViewMutType<ViewMut<'a>: Slice<View<'b> = S::View<'b>>>,
{
    #[inline]
    fn slice<R: SliceRange>(&self, range: R) -> FastaAnnotView<'_, M, S> {
        FastaAnnotView {
            header:   self.header,
            sequence: self.sequence.slice(range),
            annot:    self.annot.as_view(),
        }
    }

    #[inline]
    fn get_slice<R: SliceRange>(&self, range: R) -> Option<FastaAnnotView<'_, M, S>> {
        Some(FastaAnnotView {
            header:   self.header,
            sequence: self.sequence.get_slice(range)?,
            annot:    self.annot.as_view(),
        })
    }
}

impl<'a, M, S> SliceMut for FastaAnnotViewMut<'a, M, S>
where
    M: AssocViewMutType,
    S: AssocViewMutType<ViewMut<'a>: SliceMut>,
{
    #[inline]
    fn slice_mut<R: SliceRange>(&mut self, range: R) -> FastaAnnotViewMut<'_, M, S> {
        FastaAnnotViewMut {
            header:   self.header,
            sequence: self.sequence.slice_mut(range),
            annot:    self.annot.reborrow_view_mut(),
        }
    }

    #[inline]
    fn get_slice_mut<R: SliceRange>(&mut self, range: R) -> Option<FastaAnnotViewMut<'_, M, S>> {
        Some(FastaAnnotViewMut {
            header:   self.header,
            sequence: self.sequence.get_slice_mut(range)?,
            annot:    self.annot.reborrow_view_mut(),
        })
    }
}
