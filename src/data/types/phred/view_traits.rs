use crate::{
    data::views::{
        impl_len_for_views_generic, impl_restrict_for_wrapper, impl_slice_for_wrapper, impl_view_assoc_types_generic,
        impl_view_conversion_generic,
    },
    prelude::{QualityScores, QualityScoresView, QualityScoresViewMut},
};

impl_len_for_views_generic!(QualityScores, QualityScoresView, QualityScoresViewMut, 0);
impl_view_assoc_types_generic!(QualityScores, QualityScoresView, QualityScoresViewMut);
impl_view_conversion_generic!(QualityScores, QualityScoresView, QualityScoresViewMut);
impl_slice_for_wrapper!(QualityScores, QualityScoresView, QualityScoresViewMut);
impl_restrict_for_wrapper!(QualityScores, QualityScoresView, QualityScoresViewMut);
