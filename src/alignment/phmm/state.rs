use std::{
    marker::PhantomData,
    ops::{Index, IndexMut},
};

use crate::alignment::phmm::PhmmError;

/// An enum representing the three states within each layer of a pHMM. This is
/// used for readability when indexing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhmmState {
    Delete = 0,
    Match  = 1,
    Insert = 2,
}

/// An enum representing the three states within each layer of a pHMM, in
/// addition to `Enter`. This is useful for local pHMMs.
#[derive(Clone, Copy, Debug)]
pub enum PhmmStateOrEnter {
    Delete = 0,
    Match  = 1,
    Insert = 2,
    Enter  = 3,
}

impl From<usize> for PhmmState {
    #[inline]
    fn from(value: usize) -> Self {
        // WARNING: enum order must be maintained
        [Self::Delete, Self::Match, Self::Insert][value]
    }
}

impl From<PhmmState> for usize {
    #[inline]
    fn from(value: PhmmState) -> Self {
        value as usize
    }
}

impl From<usize> for PhmmStateOrEnter {
    fn from(value: usize) -> Self {
        // WARNING: enum order must be maintained
        [Self::Delete, Self::Match, Self::Insert, Self::Enter][value]
    }
}

impl From<PhmmStateOrEnter> for usize {
    #[inline]
    fn from(value: PhmmStateOrEnter) -> Self {
        value as usize
    }
}

impl From<PhmmState> for PhmmStateOrEnter {
    fn from(value: PhmmState) -> Self {
        match value {
            PhmmState::Delete => PhmmStateOrEnter::Delete,
            PhmmState::Match => PhmmStateOrEnter::Match,
            PhmmState::Insert => PhmmStateOrEnter::Insert,
        }
    }
}

impl PhmmState {
    #[inline]
    pub(crate) fn get_from(value: PhmmStateOrEnter) -> Option<PhmmState> {
        match value {
            PhmmStateOrEnter::Delete => Some(PhmmState::Delete),
            PhmmStateOrEnter::Match => Some(PhmmState::Match),
            PhmmStateOrEnter::Insert => Some(PhmmState::Insert),
            PhmmStateOrEnter::Enter => None,
        }
    }

    #[inline]
    pub(crate) fn to_op(self) -> u8 {
        match self {
            PhmmState::Delete => b'D',
            PhmmState::Match => b'M',
            PhmmState::Insert => b'I',
        }
    }

    #[inline]
    pub(crate) fn from_op(op: u8) -> Result<Self, PhmmError> {
        match op {
            b'D' => Ok(PhmmState::Delete),
            b'M' | b'=' | b'X' => Ok(PhmmState::Match),
            b'I' => Ok(PhmmState::Insert),
            _ => Err(PhmmError::InvalidCigarOp),
        }
    }
}

/// An array of type `[T; N]` indexed by an enum `E`. This is useful for
/// readability.
#[derive(Clone, Copy)]
pub(crate) struct EnumArray<T, E, const N: usize> {
    pub(crate) inner: [T; N],
    phantom:          PhantomData<E>,
}

impl<T, E, const N: usize> EnumArray<T, E, N> {
    /// Create a new array of type `[T; N]` indexed by enum `E`.
    #[inline]
    #[must_use]
    pub(crate) fn new(inner: [T; N]) -> Self {
        Self {
            inner,
            phantom: PhantomData,
        }
    }
}

impl<T: PartialOrd + Copy, E: From<usize>, const N: usize> EnumArray<T, E, N> {
    /// Locate the minimum value in the array and the corresponding enum variant
    pub(crate) fn locate_min(&self) -> (E, T) {
        let mut argmin = 0;
        let mut min = self.inner[0];

        for i in 1..N {
            if self.inner[i] < min {
                argmin = i;
                min = self.inner[i];
            }
        }

        (E::from(argmin), min)
    }
}

impl<T: Default + Copy, E, const N: usize> Default for EnumArray<T, E, N> {
    #[inline]
    fn default() -> Self {
        Self {
            inner:   [T::default(); N],
            phantom: PhantomData,
        }
    }
}

impl<T, E: Into<usize>, const N: usize> Index<E> for EnumArray<T, E, N> {
    type Output = T;

    #[inline]
    fn index(&self, index: E) -> &Self::Output {
        &self.inner[index.into()]
    }
}

impl<T, E: Into<usize>, const N: usize> IndexMut<E> for EnumArray<T, E, N> {
    #[inline]
    fn index_mut(&mut self, index: E) -> &mut Self::Output {
        &mut self.inner[index.into()]
    }
}

/// An array holding values of type `T`, indexed by the variants in
/// [`PhmmState`]
pub(crate) type PhmmStateArray<T> = EnumArray<T, PhmmState, 3>;

/// An array holding values of type `T`, indexed by the variants in
/// [`PhmmStateOrEnter`]
pub(crate) type PhmmStateOrEnterArray<T> = EnumArray<T, PhmmStateOrEnter, 4>;

impl<T: Copy> PhmmStateArray<T> {
    pub(crate) fn with_enter(self, enter: T) -> PhmmStateOrEnterArray<T> {
        PhmmStateOrEnterArray::new([self.inner[0], self.inner[1], self.inner[2], enter])
    }
}

/// A trait for enums of which [`PhmmState`] is a subset (and which can be
/// stored as usize)
pub(crate) trait PhmmStateEnum: From<PhmmState> + Into<usize> + From<usize> {}
impl<T: From<PhmmState> + Into<usize> + From<usize>> PhmmStateEnum for T {}
