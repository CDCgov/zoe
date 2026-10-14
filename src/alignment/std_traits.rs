use crate::alignment::{LocalProfiles, SharedProfiles};
use std::hash::{Hash, Hasher};

impl<const M: usize, const N: usize, const O: usize, const S: usize> PartialEq for LocalProfiles<'_, M, N, O, S> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        (self.seq == other.seq) && (self.params == other.params)
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> Eq for LocalProfiles<'_, M, N, O, S> {}

impl<const M: usize, const N: usize, const O: usize, const S: usize> Hash for LocalProfiles<'_, M, N, O, S> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seq.hash(state);
        self.params.hash(state);
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> PartialEq for SharedProfiles<'_, M, N, O, S> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        (self.seq == other.seq) && (self.params == other.params)
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> Eq for SharedProfiles<'_, M, N, O, S> {}

impl<const M: usize, const N: usize, const O: usize, const S: usize> Hash for SharedProfiles<'_, M, N, O, S> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seq.hash(state);
        self.params.hash(state);
    }
}
