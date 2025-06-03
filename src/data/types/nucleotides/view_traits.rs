use crate::{
    data::views::{impl_len, impl_restrict_for_wrapper, impl_views_for_wrapper},
    prelude::*,
};

impl_len!(Nucleotides, NucleotidesView, NucleotidesViewMut, 0);
impl_views_for_wrapper!(Nucleotides, NucleotidesView, NucleotidesViewMut);
impl_restrict_for_wrapper!(Nucleotides, NucleotidesView, NucleotidesViewMut);
