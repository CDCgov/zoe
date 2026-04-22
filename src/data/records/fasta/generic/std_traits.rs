use crate::data::{
    fasta::generic::{Fasta, FastaAnnotView, FastaView, FastaViewMut},
    views::{AssocViewMutType, AssocViewType},
};
use std::fmt::Display;

impl<S: Display> Display for Fasta<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.header, self.sequence)
    }
}

impl<'a, S> Display for FastaView<'a, S>
where
    S: AssocViewType<View<'a>: Display>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.header, self.sequence)
    }
}

impl<'a, S> Display for FastaViewMut<'a, S>
where
    S: AssocViewMutType<ViewMut<'a>: Display>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.header, self.sequence)
    }
}

// A non-derived impl is needed since S may not implement clone (but S::View
// does)
impl<'a, S> Clone for FastaView<'a, S>
where
    S: AssocViewType<View<'a>: Clone>,
{
    #[inline]
    fn clone(&self) -> Self {
        Self {
            header:   self.header,
            sequence: self.sequence.clone(),
        }
    }
}

impl<'a, S> Copy for FastaView<'a, S> where S: AssocViewType<View<'a>: Copy> {}

// A non-derived impl is needed since M and S may not implement clone (but
// S::View and M::View do)
impl<'a, M, S> Clone for FastaAnnotView<'a, M, S>
where
    M: AssocViewType<View<'a>: Clone>,
    S: AssocViewType<View<'a>: Clone>,
{
    #[inline]
    fn clone(&self) -> Self {
        Self {
            header:   self.header,
            sequence: self.sequence.clone(),
            annot:    self.annot.clone(),
        }
    }
}

impl<'a, M, S> Copy for FastaAnnotView<'a, M, S>
where
    M: AssocViewType<View<'a>: Copy>,
    S: AssocViewType<View<'a>: Copy>,
{
}
