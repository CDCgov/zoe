use crate::{
    data::view_traits::{SliceRange, impl_len, impl_views_for_wrapper},
    prelude::*,
};

impl_len!(QualityScores, QualityScoresView, QualityScoresViewMut, 0);
impl_views_for_wrapper!(QualityScores, QualityScoresView, QualityScoresViewMut);
