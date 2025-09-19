use super::*;
use crate::{
    data::{sam::SamData, views::impl_len_for_views},
    prelude::{DataOwned, DataView, DataViewMut},
};

impl_len_for_views!(SamData, SamDataView, SamDataViewMut, seq);

impl DataOwned for SamData {
    type View<'a>
        = SamDataView<'a>
    where
        Self: 'a;

    type ViewMut<'a>
        = SamDataViewMut<'a>
    where
        Self: 'a;

    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        SamDataView::new(
            &self.qname,
            self.flag,
            &self.rname,
            self.pos,
            self.mapq,
            self.cigar.as_view(),
            self.seq.as_view(),
            self.qual.as_view(),
        )
    }

    #[inline]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_> {
        SamDataViewMut::new(
            &mut self.qname,
            self.flag,
            &mut self.rname,
            self.pos,
            self.mapq,
            self.cigar.as_view_mut(),
            self.seq.as_view_mut(),
            self.qual.as_view_mut(),
        )
    }
}

impl DataView for SamDataView<'_> {
    type Owned = SamData;

    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        SamData::new(
            self.qname.to_string(),
            self.flag,
            self.rname.to_string(),
            self.pos,
            self.mapq,
            self.cigar.to_owned_data(),
            self.seq.to_owned_data(),
            self.qual.to_owned_data(),
        )
    }
}

impl<'b> DataViewMut<'b> for SamDataViewMut<'b> {
    type View<'a>
        = SamDataView<'a>
    where
        Self: 'a;

    type Owned = SamData;

    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        SamDataView::new(
            self.qname,
            self.flag,
            self.rname,
            self.pos,
            self.mapq,
            self.cigar.as_view(),
            self.seq.as_view(),
            self.qual.as_view(),
        )
    }

    #[inline]
    fn to_view(self) -> Self::View<'b> {
        SamDataView::new(
            self.qname,
            self.flag,
            self.rname,
            self.pos,
            self.mapq,
            self.cigar.to_view(),
            self.seq.to_view(),
            self.qual.to_view(),
        )
    }

    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        SamData::new(
            (*self.qname).clone(),
            self.flag,
            (*self.rname).clone(),
            self.pos,
            self.mapq,
            self.cigar.to_owned_data(),
            self.seq.to_owned_data(),
            self.qual.to_owned_data(),
        )
    }
}
