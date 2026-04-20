use crate::{
    data::{
        fasta::generic::{Fasta, FastaAnnot, FastaAnnotView, FastaAnnotViewMut, FastaView, FastaViewMut},
        views::{SliceRange, ViewAssocTypes},
    },
    prelude::{DataOwned, DataView, DataViewMut, Len, Restrict, Slice, SliceMut},
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
    S: ViewAssocTypes<View<'a>: Len>,
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
    S: ViewAssocTypes<ViewMut<'a>: Len>,
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

impl<S> ViewAssocTypes for Fasta<S>
where
    S: DataOwned,
{
    type Owned = Fasta<S>;
    type View<'a> = FastaView<'a, S>;
    type ViewMut<'a> = FastaViewMut<'a, S>;
}

impl<S> ViewAssocTypes for FastaView<'_, S>
where
    S: DataOwned,
{
    type Owned = Fasta<S>;
    type View<'a> = FastaView<'a, S>;
    type ViewMut<'a> = FastaViewMut<'a, S>;
}

impl<S> ViewAssocTypes for FastaViewMut<'_, S>
where
    S: DataOwned,
{
    type Owned = Fasta<S>;

    type View<'a> = FastaView<'a, S>;

    type ViewMut<'a> = FastaViewMut<'a, S>;
}

impl<S: DataOwned> DataOwned for Fasta<S> {
    #[inline]
    fn as_view(&self) -> FastaView<'_, S> {
        FastaView {
            header:   &self.header,
            sequence: self.sequence.as_view(),
        }
    }

    #[inline]
    fn as_view_mut(&mut self) -> FastaViewMut<'_, S> {
        FastaViewMut {
            header:   &mut self.header,
            sequence: self.sequence.as_view_mut(),
        }
    }
}

impl<'a, S> DataView<'a> for FastaView<'a, S>
where
    S: DataOwned,
{
    #[inline]
    fn to_owned_data(&self) -> Fasta<S> {
        Fasta {
            header:   self.header.to_string(),
            sequence: self.sequence.to_owned_data(),
        }
    }

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
    S: DataOwned,
{
    #[inline]
    fn as_view(&self) -> FastaView<'_, S> {
        FastaView {
            header:   self.header,
            sequence: self.sequence.as_view(),
        }
    }

    #[inline]
    fn to_view(self) -> FastaView<'a, S> {
        FastaView {
            header:   self.header,
            sequence: self.sequence.to_view(),
        }
    }

    #[inline]
    fn to_owned_data(&self) -> Fasta<S> {
        Fasta {
            header:   (*self.header).clone(),
            sequence: self.sequence.to_owned_data(),
        }
    }

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
    S: ViewAssocTypes<View<'a>: Restrict>,
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
    S: ViewAssocTypes<ViewMut<'a>: Restrict>,
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
    S: DataOwned + Slice,
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
    S: DataOwned + SliceMut,
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
    S: DataOwned<View<'a>: Slice> + 'a,
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
    S: DataOwned<ViewMut<'a>: Slice> + 'a,
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
    S: DataOwned<ViewMut<'a>: SliceMut> + 'a,
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
    M: ViewAssocTypes,
    S: ViewAssocTypes<View<'a>: Len>,
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
    M: ViewAssocTypes,
    S: ViewAssocTypes<ViewMut<'a>: Len>,
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

impl<M, S> ViewAssocTypes for FastaAnnot<M, S>
where
    M: DataOwned,
    S: DataOwned,
{
    type Owned = FastaAnnot<M, S>;
    type View<'a> = FastaAnnotView<'a, M, S>;
    type ViewMut<'a> = FastaAnnotViewMut<'a, M, S>;
}

impl<M, S> ViewAssocTypes for FastaAnnotView<'_, M, S>
where
    M: DataOwned,
    S: DataOwned,
{
    type Owned = FastaAnnot<M, S>;
    type View<'a> = FastaAnnotView<'a, M, S>;
    type ViewMut<'a> = FastaAnnotViewMut<'a, M, S>;
}

impl<M, S> ViewAssocTypes for FastaAnnotViewMut<'_, M, S>
where
    M: DataOwned,
    S: DataOwned,
{
    type Owned = FastaAnnot<M, S>;
    type View<'a> = FastaAnnotView<'a, M, S>;
    type ViewMut<'a> = FastaAnnotViewMut<'a, M, S>;
}

impl<M, S> DataOwned for FastaAnnot<M, S>
where
    M: DataOwned,
    S: DataOwned,
{
    #[inline]
    fn as_view(&self) -> FastaAnnotView<'_, M, S> {
        FastaAnnotView {
            header:   &self.header,
            sequence: self.sequence.as_view(),
            annot:    self.annot.as_view(),
        }
    }

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
    M: DataOwned<View<'a>: DataView<'a, Owned = M>>,
    S: DataOwned<View<'a>: DataView<'a, Owned = S>>,
{
    #[inline]
    fn to_owned_data(&self) -> FastaAnnot<M, S> {
        FastaAnnot {
            header:   self.header.to_string(),
            sequence: self.sequence.to_owned_data(),
            annot:    self.annot.to_owned_data(),
        }
    }

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

impl<'a, M, S> Restrict for FastaAnnotView<'a, M, S>
where
    M: DataOwned,
    S: DataOwned<View<'a>: Restrict>,
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

impl<'a, M, S> DataViewMut<'a> for FastaAnnotViewMut<'a, M, S>
where
    M: DataOwned,
    S: DataOwned,
{
    #[inline]
    fn as_view(&self) -> FastaAnnotView<'_, M, S> {
        FastaAnnotView {
            header:   self.header,
            sequence: self.sequence.as_view(),
            annot:    self.annot.as_view(),
        }
    }

    #[inline]
    fn to_view(self) -> FastaAnnotView<'a, M, S> {
        FastaAnnotView {
            header:   self.header,
            sequence: self.sequence.to_view(),
            annot:    self.annot.to_view(),
        }
    }

    #[inline]
    fn to_owned_data(&self) -> FastaAnnot<M, S> {
        FastaAnnot {
            header:   (*self.header).clone(),
            sequence: self.sequence.to_owned_data(),
            annot:    self.annot.to_owned_data(),
        }
    }

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

impl<'a, M, S> Restrict for FastaAnnotViewMut<'a, M, S>
where
    M: DataOwned,
    S: DataOwned<ViewMut<'a>: Restrict>,
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
    M: DataOwned,
    S: DataOwned + Slice,
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
    M: DataOwned,
    S: DataOwned + SliceMut,
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
    M: DataOwned,
    S: DataOwned<View<'a>: Slice> + 'a,
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
    M: DataOwned,
    S: DataOwned<ViewMut<'a>: SliceMut>,
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
    M: DataOwned,
    S: DataOwned<ViewMut<'a>: SliceMut>,
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
