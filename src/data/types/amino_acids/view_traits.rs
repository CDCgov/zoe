use crate::{
    data::view_traits::{SliceRange, impl_len, impl_views_for_wrapper},
    prelude::*,
};

impl_len!(AminoAcids, AminoAcidsView, AminoAcidsViewMut, 0);
impl_views_for_wrapper!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
