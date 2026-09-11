//! Structs and traits to enable more readable/correct indexing into
//! pHMM-related data structures.

use crate::alignment::phmm::{
    DomainPhmm, GlobalPhmm, LocalPhmm, SemiLocalPhmm,
    components::CorePhmm,
    indexing::{GetCore, GetLayer},
    modules::{PrecomputedDomainModule, PrecomputedLocalModule, SemiLocalModule},
};

/// A trait for structures that can be indexed via a [`AlnIndex`], such as byte
/// sequences, pHMMs, and modules.
///
/// [`AlnIndex`]: crate::alignment::phmm::indexing::AlnIndex
pub trait AlnIndexable {
    /// Returns the length of the query or reference sequence corresponding to
    /// the structure.
    #[must_use]
    fn seq_len(&self) -> usize;
}

impl<P: AlnIndexable> AlnIndexable for &P {
    #[inline]
    fn seq_len(&self) -> usize {
        P::seq_len(self)
    }
}

impl<P: AlnIndexable> AlnIndexable for &mut P {
    #[inline]
    fn seq_len(&self) -> usize {
        P::seq_len(self)
    }
}

impl<T, const S: usize> AlnIndexable for CorePhmm<T, S> {
    #[inline]
    fn seq_len(&self) -> usize {
        // The END state does not have an index in the CorePhmm, so we subtract
        // one just for the BEGIN state
        self.layers().len() - 1
    }
}

impl<T, const S: usize> AlnIndexable for GlobalPhmm<T, S> {
    #[inline]
    fn seq_len(&self) -> usize {
        self.core().seq_len()
    }
}

impl<T, const S: usize> AlnIndexable for LocalPhmm<T, S> {
    #[inline]
    fn seq_len(&self) -> usize {
        self.core().seq_len()
    }
}

impl<T, const S: usize> AlnIndexable for SemiLocalPhmm<T, S> {
    #[inline]
    fn seq_len(&self) -> usize {
        self.core().seq_len()
    }
}

impl<T, const S: usize> AlnIndexable for DomainPhmm<T, S> {
    #[inline]
    fn seq_len(&self) -> usize {
        self.core().seq_len()
    }
}

impl<T> AlnIndexable for SemiLocalModule<T> {
    #[inline]
    fn seq_len(&self) -> usize {
        self.0.len() - 2
    }
}

impl<T, const S: usize> AlnIndexable for PrecomputedLocalModule<'_, T, S> {
    #[inline]
    fn seq_len(&self) -> usize {
        self.semilocal_params.seq_len()
    }
}

impl<T, const S: usize> AlnIndexable for PrecomputedDomainModule<T, S> {
    #[inline]
    fn seq_len(&self) -> usize {
        self.0.len() - 1
    }
}

impl AlnIndexable for [u8] {
    #[inline]
    fn seq_len(&self) -> usize {
        self.len()
    }
}

impl AlnIndexable for &[u8] {
    #[inline]
    fn seq_len(&self) -> usize {
        self.len()
    }
}

/// A trait providing an extension of [`AlnIndexable`] specifically for pHMMs.
pub trait PhmmLen: AlnIndexable {
    /// Returns the number of pseudo-match states in the pHMM, which includes
    /// both match states with emissions as well as the BEGIN and END states.
    fn num_pseudomatch(&self) -> usize {
        self.seq_len() + 2
    }
}

impl<T, const S: usize> PhmmLen for CorePhmm<T, S> {}
impl<T, const S: usize> PhmmLen for GlobalPhmm<T, S> {}
impl<T, const S: usize> PhmmLen for LocalPhmm<T, S> {}
impl<T, const S: usize> PhmmLen for SemiLocalPhmm<T, S> {}
impl<T, const S: usize> PhmmLen for DomainPhmm<T, S> {}
impl<T> PhmmLen for SemiLocalModule<T> {}
impl<T, const S: usize> PhmmLen for PrecomputedLocalModule<'_, T, S> {}
impl<T, const S: usize> PhmmLen for PrecomputedDomainModule<T, S> {}
