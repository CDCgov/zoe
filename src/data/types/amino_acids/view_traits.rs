use crate::{
    data::views::{impl_len_for_views, impl_restrict_for_wrapper, impl_slice_for_wrapper, impl_views_for_wrapper},
    prelude::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
};

impl_len_for_views!(AminoAcids, AminoAcidsView, AminoAcidsViewMut, 0);
impl_views_for_wrapper!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
impl_slice_for_wrapper!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
impl_restrict_for_wrapper!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
