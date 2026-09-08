use super::*;
use crate::{
    data::{
        sam::SamData,
        views::{
            AsView, AsViewMut, AssocOwnedType, AssocViewType, ToOwnedData, ToView, impl_len_for_views_generic,
            impl_view_assoc_types_generic,
        },
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
            opt_fields: self.opt_fields.to_owned_data(),
        }
    }
}

impl ToOwnedData for SamDataViewMut<'_> {
    #[inline]
    fn to_owned_data(&self) -> SamData {
        SamData {
            qname:      (*self.qname).clone(),
            flag:       self.flag,
            rname:      (*self.rname).clone(),
            pos:        self.pos,
            mapq:       self.mapq,
            cigar:      self.cigar.to_owned_data(),
            rnext:      self.rnext,
            pnext:      self.pnext,
            tlen:       self.tlen,
            seq:        self.seq.to_owned_data(),
            qual:       self.qual.to_owned_data(),
            // TODO: SamDataViewMut doesn't currently contain tags
            opt_fields: SamOptRaw::new(),
        }
    }
}

impl AsView for SamData {
    #[inline]
    fn as_view(&self) -> SamDataView<'_> {
        SamDataView {
            qname:      &self.qname,
            flag:       self.flag,
            rname:      &self.rname,
            pos:        self.pos,
            mapq:       self.mapq,
            cigar:      self.cigar.as_view(),
            rnext:      self.rnext,
            pnext:      self.pnext,
            tlen:       self.tlen,
            seq:        self.seq.as_view(),
            qual:       self.qual.as_view(),
            opt_fields: self.opt_fields.as_view(),
        }
    }
}

impl AsView for SamDataViewMut<'_> {
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        SamDataView {
            qname:      self.qname,
            flag:       self.flag,
            rname:      self.rname,
            pos:        self.pos,
            mapq:       self.mapq,
            cigar:      self.cigar.as_view(),
            rnext:      self.rnext,
            pnext:      self.pnext,
            tlen:       self.tlen,
            seq:        self.seq.as_view(),
            qual:       self.qual.as_view(),
            // TODO: SamDataViewMut doesn't currently contain tags
            opt_fields: SamOptRawView::new(),
        }
    }
}

impl AsViewMut for SamData {
    #[inline]
    fn as_view_mut(&mut self) -> Self::ViewMut<'_> {
        SamDataViewMut {
            qname: &mut self.qname,
            flag:  self.flag,
            rname: &mut self.rname,
            pos:   self.pos,
            mapq:  self.mapq,
            cigar: self.cigar.as_view_mut(),
            rnext: self.rnext,
            pnext: self.pnext,
            tlen:  self.tlen,
            seq:   self.seq.as_view_mut(),
            qual:  self.qual.as_view_mut(),
        }
    }
}

impl<'a> ToView<'a> for SamDataViewMut<'a> {
    #[inline]
    fn to_view(self) -> SamDataView<'a> {
        SamDataView {
            qname:      self.qname,
            flag:       self.flag,
            rname:      self.rname,
            pos:        self.pos,
            mapq:       self.mapq,
            cigar:      self.cigar.to_view(),
            rnext:      self.rnext,
            pnext:      self.pnext,
            tlen:       self.tlen,
            seq:        self.seq.to_view(),
            qual:       self.qual.to_view(),
            // TODO: SamDataViewMut doesn't currently contain tags
            opt_fields: SamOptRawView::new(),
        }
    }
}

impl<'a> DataView<'a> for SamDataView<'a> {
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        SamDataView {
            qname:      self.qname,
            flag:       self.flag,
            rname:      self.rname,
            pos:        self.pos,
            mapq:       self.mapq,
            cigar:      self.cigar.reborrow_view(),
            rnext:      self.rnext,
            pnext:      self.pnext,
            tlen:       self.tlen,
            seq:        self.seq.reborrow_view(),
            qual:       self.qual.reborrow_view(),
            opt_fields: self.opt_fields.reborrow_view(),
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

impl AssocViewType for SamOptRaw {
    type View<'a> = SamOptRawView<'a>;
}

impl AssocViewType for SamOptRawView<'_> {
    type View<'a> = SamOptRawView<'a>;
}

impl AssocOwnedType for SamOptRaw {
    type Owned = SamOptRaw;
}

impl AssocOwnedType for SamOptRawView<'_> {
    type Owned = SamOptRaw;
}

impl<'a> DataView<'a> for SamOptRawView<'a> {
    #[inline]
    fn reborrow_view<'b>(&'b self) -> Self::View<'b>
    where
        'a: 'b, {
        SamOptRawView(self.0)
    }
}

impl AsView for SamOptRaw {
    #[inline]
    fn as_view(&self) -> Self::View<'_> {
        SamOptRawView(&self.0)
    }
}

impl ToOwnedData for SamOptRawView<'_> {
    #[inline]
    fn to_owned_data(&self) -> Self::Owned {
        SamOptRaw(self.0.to_vec())
    }
}
