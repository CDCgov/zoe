use crate::data::{
    fasta::generic::{Fasta, FastaAnnotView, FastaView},
    views::ViewAssocTypes,
};
use std::fmt::Display;

impl<S: Display> Display for Fasta<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.header, self.sequence)
    }
}

// A non-derived impl is needed since S may not implement clone (but S::View
// does)
impl<'a, S> Clone for FastaView<'a, S>
where
    S: ViewAssocTypes<View<'a>: Clone>,
{
    #[inline]
    fn clone(&self) -> Self {
        Self {
            header:   self.header,
            sequence: self.sequence.clone(),
        }
    }
}

impl<'a, S> Copy for FastaView<'a, S> where S: ViewAssocTypes<View<'a>: Copy> {}

// A non-derived impl is needed since M and S may not implement clone (but
// S::View and M::View do)
impl<'a, M, S> Clone for FastaAnnotView<'a, M, S>
where
    M: ViewAssocTypes<View<'a>: Clone>,
    S: ViewAssocTypes<View<'a>: Clone>,
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
    M: ViewAssocTypes<View<'a>: Copy>,
    S: ViewAssocTypes<View<'a>: Copy>,
{
}
