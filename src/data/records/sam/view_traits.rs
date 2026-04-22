use super::*;
use crate::{
    data::{
        sam::SamData,
        views::{AsView, AsViewMut, ToOwnedData, ToView, impl_len_for_views_generic, impl_view_assoc_types_generic},
    },
    prelude::{DataView, DataViewMut},
};

impl_len_for_views_generic!(SamData, SamDataView, SamDataViewMut, seq);
impl_view_assoc_types_generic!(SamData, SamDataView, SamDataViewMut);

impl ToOwnedData for SamDataView<'_> {
    #[inline]
    fn to_owned_data(&self) -> SamData {
        SamData {
            qname:      self.qname.to_string(),
            flag:       self.flag,
            rname:      self.rname.to_string(),
            pos:        self.pos,
            mapq:       self.mapq,
            cigar:      self.cigar.to_owned_data(),
            rnext:      self.rnext,
            pnext:      self.pnext,
            tlen:       self.tlen,
            seq:        self.seq.to_owned_data(),
            qual:       self.qual.to_owned_data(),
            opt_fields: SamOptRaw::new(),
        }
    }
}

impl ToOwnedData for SamDataViewMut<'_> {
    #[inline]
    fn to_owned_data(&self) -> SamData {
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

impl AsView for SamData {
    #[inline]
    fn as_view(&self) -> SamDataView<'_> {
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
}

impl AsView for SamDataViewMut<'_> {
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
}

impl AsViewMut for SamData {
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

impl<'a> ToView<'a> for SamDataViewMut<'a> {
    #[inline]
    fn to_view(self) -> SamDataView<'a> {
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
}

impl<'a> DataView<'a> for SamDataView<'a> {
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        SamDataView {
            qname: self.qname,
            flag:  self.flag,
            rname: self.rname,
            pos:   self.pos,
            mapq:  self.mapq,
            cigar: self.cigar.reborrow_view(),
            rnext: self.rnext,
            pnext: self.pnext,
            tlen:  self.tlen,
            seq:   self.seq.reborrow_view(),
            qual:  self.qual.reborrow_view(),
        }
    }
}

impl<'a> DataViewMut<'a> for SamDataViewMut<'a> {
    #[inline]
    fn reborrow_view_mut<'b>(&'b mut self) -> Self::ViewMut<'b>
    where
        'a: 'b, {
        SamDataViewMut {
            qname: self.qname,
            flag:  self.flag,
            rname: self.rname,
            pos:   self.pos,
            mapq:  self.mapq,
            cigar: self.cigar.reborrow_view_mut(),
            rnext: self.rnext,
            pnext: self.pnext,
            tlen:  self.tlen,
            seq:   self.seq.reborrow_view_mut(),
            qual:  self.qual.reborrow_view_mut(),
        }
    }
}
