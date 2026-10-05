use super::*;
use crate::{
    data::{
        sam::SamData,
        views::{AsView, AssocOwnedType, AssocViewType, ToOwnedData, impl_len_for_wrapper},
    },
    prelude::DataView,
};

impl_len_for_wrapper!(SamData, seq);
impl_len_for_wrapper!(SamDataView<'_>, seq);

impl AssocOwnedType for SamData {
    type Owned = SamData;
}

impl AssocViewType for SamData {
    type View<'a> = SamDataView<'a>;
}

impl AssocOwnedType for SamDataView<'_> {
    type Owned = SamData;
}

impl AssocViewType for SamDataView<'_> {
    type View<'a> = SamDataView<'a>;
}

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
        SamOptRaw(self.0.to_string())
    }
}
