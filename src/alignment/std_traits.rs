use crate::alignment::{LocalProfiles, SharedProfiles};
use std::{
    hash::{Hash, Hasher},
    simd::{LaneCount, SupportedLaneCount},
};

impl<const M: usize, const N: usize, const O: usize, const S: usize> PartialEq for LocalProfiles<'_, M, N, O, S>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        (self.seq == other.seq)
            && (self.matrix == other.matrix)
            && (self.gap_open == other.gap_open)
            && (self.gap_extend == other.gap_extend)
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> Eq for LocalProfiles<'_, M, N, O, S>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
{
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> Hash for LocalProfiles<'_, M, N, O, S>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
{
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seq.hash(state);
        self.matrix.hash(state);
        self.gap_open.hash(state);
        self.gap_extend.hash(state);
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> PartialEq for SharedProfiles<'_, M, N, O, S>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        (self.seq == other.seq)
            && (self.matrix == other.matrix)
            && (self.gap_open == other.gap_open)
            && (self.gap_extend == other.gap_extend)
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> Eq for SharedProfiles<'_, M, N, O, S>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
{
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> Hash for SharedProfiles<'_, M, N, O, S>
where
    LaneCount<M>: SupportedLaneCount,
    LaneCount<N>: SupportedLaneCount,
    LaneCount<O>: SupportedLaneCount,
{
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seq.hash(state);
        self.matrix.hash(state);
        self.gap_open.hash(state);
        self.gap_extend.hash(state);
    }
}
