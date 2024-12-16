use super::{LocalProfile, SharedProfile};
use std::{
    hash::{Hash, Hasher},
    simd::{LaneCount, SupportedLaneCount},
};

impl<const N: usize, const S: usize> PartialEq for LocalProfile<'_, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn eq(&self, other: &Self) -> bool {
        (self.query == other.query)
            && (self.matrix == other.matrix)
            && (self.gap_open == other.gap_open)
            && (self.gap_extend == other.gap_extend)
    }
}

impl<const N: usize, const S: usize> Eq for LocalProfile<'_, N, S> where LaneCount<N>: SupportedLaneCount {}

impl<const N: usize, const S: usize> Hash for LocalProfile<'_, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.query.hash(state);
        self.matrix.hash(state);
        self.gap_open.hash(state);
        self.gap_extend.hash(state);
    }
}

impl<const N: usize, const S: usize> PartialEq for SharedProfile<'_, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn eq(&self, other: &Self) -> bool {
        (self.query == other.query)
            && (self.matrix == other.matrix)
            && (self.gap_open == other.gap_open)
            && (self.gap_extend == other.gap_extend)
    }
}

impl<const N: usize, const S: usize> Eq for SharedProfile<'_, N, S> where LaneCount<N>: SupportedLaneCount {}

impl<const N: usize, const S: usize> Hash for SharedProfile<'_, N, S>
where
    LaneCount<N>: SupportedLaneCount,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.query.hash(state);
        self.matrix.hash(state);
        self.gap_open.hash(state);
        self.gap_extend.hash(state);
    }
}
