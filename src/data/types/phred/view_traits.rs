use crate::{
    data::views::{impl_len_for_views, impl_restrict_for_wrapper, impl_slice_for_wrapper, impl_views_for_wrapper},
    prelude::{QualityScores, QualityScoresView, QualityScoresViewMut},
};

impl_len_for_views!(QualityScores, QualityScoresView, QualityScoresViewMut, 0);
impl_views_for_wrapper!(QualityScores, QualityScoresView, QualityScoresViewMut);
impl_slice_for_wrapper!(QualityScores, QualityScoresView, QualityScoresViewMut);
impl_restrict_for_wrapper!(QualityScores, QualityScoresView, QualityScoresViewMut);
