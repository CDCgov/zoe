use crate::{
    data::views::{
        impl_len_for_views_generic, impl_restrict_for_wrapper, impl_slice_for_wrapper, impl_view_assoc_types_generic,
        impl_view_conversion_generic,
    },
    prelude::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
};

impl_len_for_views_generic!(AminoAcids, AminoAcidsView, AminoAcidsViewMut, 0);
impl_view_assoc_types_generic!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
impl_view_conversion_generic!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
impl_slice_for_wrapper!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
impl_restrict_for_wrapper!(AminoAcids, AminoAcidsView, AminoAcidsViewMut);
