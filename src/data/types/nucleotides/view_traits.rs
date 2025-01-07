use crate::{
    data::view_traits::{SliceRange, impl_len, impl_views_for_wrapper},
    prelude::*,
};

impl_len!(Nucleotides, NucleotidesView, NucleotidesViewMut, 0);
impl_views_for_wrapper!(Nucleotides, NucleotidesView, NucleotidesViewMut);
