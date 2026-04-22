use crate::{
    data::views::{
        impl_len_for_views_generic, impl_restrict_for_wrapper, impl_slice_for_wrapper, impl_view_assoc_types_generic,
        impl_view_conversion_generic,
    },
    prelude::*,
};

impl_len_for_views_generic!(Nucleotides, NucleotidesView, NucleotidesViewMut, 0);
impl_view_assoc_types_generic!(Nucleotides, NucleotidesView, NucleotidesViewMut);
impl_view_conversion_generic!(Nucleotides, NucleotidesView, NucleotidesViewMut);
impl_slice_for_wrapper!(Nucleotides, NucleotidesView, NucleotidesViewMut);
impl_restrict_for_wrapper!(Nucleotides, NucleotidesView, NucleotidesViewMut);
