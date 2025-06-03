use crate::{
    data::views::{impl_len, impl_restrict_for_wrapper, impl_views_for_wrapper},
    prelude::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
};

impl_len!(AminoAcids, AminoAcidsView, AminoAcidsViewMut, 0);
impl_views_for_wrapper!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
impl_restrict_for_wrapper!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
