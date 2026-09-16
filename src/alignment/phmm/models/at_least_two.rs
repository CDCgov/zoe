//! A module providing a [`Vec`] wrapper guaranteed to be of length at least 2.
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]

use std::{
    borrow::{Borrow, BorrowMut},
    collections::TryReserveError,
    error::Error,
    fmt::Display,
    ops::{Deref, DerefMut, Index, IndexMut, RangeBounds},
    slice::{Iter, IterMut},
    vec::IntoIter,
};

/// A wrapper around a [`Vec`] guaranteeing that it has at least two elements.
///
/// This offers infallible methods for [`first`], [`last`], [`split_first`], and
/// so on.
///
/// [`first`]: VecAtLeast2::first
/// [`last`]: VecAtLeast2::last
/// [`split_first`]: VecAtLeast2::split_first
#[repr(transparent)]
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct VecAtLeast2<T>(Vec<T>);

impl<T> VecAtLeast2<T> {
    #[inline]
    pub fn push(&mut self, value: T) {
        self.0.push(value);
    }

    #[inline]
    pub fn push_mut(&mut self, value: T) -> &mut T {
        self.0.push_mut(value)
    }

    #[inline]
    #[must_use]
    pub const fn capacity(&self) -> usize {
        self.0.capacity()
    }

    #[inline]
    pub fn reserve(&mut self, additional: usize) {
        self.0.reserve(additional);
    }

    #[inline]
    pub fn reserve_exact(&mut self, additional: usize) {
        self.0.reserve_exact(additional);
    }

    #[inline]
    pub fn try_reserve(&mut self, additional: usize) -> Result<(), TryReserveError> {
        self.0.try_reserve(additional)
    }

    #[inline]
    pub fn try_reserve_exact(&mut self, additional: usize) -> Result<(), TryReserveError> {
        self.0.try_reserve_exact(additional)
    }

    #[inline]
    pub fn shrink_to_fit(&mut self) {
        self.0.shrink_to_fit();
    }

    #[inline]
    pub fn shrink_to(&mut self, min_capacity: usize) {
        self.0.shrink_to(min_capacity);
    }

    #[inline]
    #[must_use]
    pub fn into_boxed_slice(self) -> Box<[T]> {
        self.0.into_boxed_slice()
    }

    #[inline]
    #[must_use]
    pub const fn as_slice(&self) -> &[T] {
        self.0.as_slice()
    }

    #[inline]
    #[must_use]
    pub const fn as_mut_slice(&mut self) -> &mut [T] {
        self.0.as_mut_slice()
    }

    #[inline]
    pub fn insert(&mut self, index: usize, element: T) {
        self.0.insert(index, element);
    }

    #[inline]
    pub fn insert_mut(&mut self, index: usize, element: T) -> &mut T {
        self.0.insert_mut(index, element)
    }

    #[inline]
    pub fn append(&mut self, other: &mut Vec<T>) {
        self.0.append(other);
    }

    #[inline]
    #[must_use]
    pub fn first(&self) -> &T {
        self.0.first().expect("at least two elements are present")
    }

    #[inline]
    #[must_use]
    pub fn first_mut(&mut self) -> &mut T {
        self.0.first_mut().expect("at least two elements are present")
    }

    #[inline]
    #[must_use]
    pub fn second(&self) -> &T {
        &self.0[1]
    }

    #[inline]
    #[must_use]
    pub fn second_mut(&mut self) -> &mut T {
        &mut self.0[1]
    }

    #[inline]
    #[must_use]
    pub fn split_first(&self) -> (&T, &[T]) {
        self.0.split_first().expect("at least two elements are present")
    }

    #[inline]
    #[must_use]
    pub fn split_first_mut(&mut self) -> (&mut T, &mut [T]) {
        self.0.split_first_mut().expect("at least two elements are present")
    }

    #[inline]
    #[must_use]
    pub fn split_last(&self) -> (&T, &[T]) {
        self.0.split_last().expect("at least two elements are present")
    }

    #[inline]
    #[must_use]
    pub fn split_last_mut(&mut self) -> (&mut T, &mut [T]) {
        self.0.split_last_mut().expect("at least two elements are present")
    }

    #[inline]
    #[must_use]
    pub fn last(&self) -> &T {
        self.0.last().expect("at least two elements are present")
    }

    #[inline]
    #[must_use]
    pub fn last_mut(&mut self) -> &mut T {
        self.0.last_mut().expect("at least two elements are present")
    }
}

impl<T: Clone> VecAtLeast2<T> {
    #[inline]
    pub fn extend_from_slice(&mut self, other: &[T]) {
        self.0.extend_from_slice(other);
    }

    #[inline]
    pub fn extend_from_within<R>(&mut self, src: R)
    where
        R: RangeBounds<usize>, {
        self.0.extend_from_within(src);
    }
}

impl<T> AsMut<[T]> for VecAtLeast2<T> {
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.0.as_mut()
    }
}

impl<T> AsRef<Vec<T>> for VecAtLeast2<T> {
    #[inline]
    fn as_ref(&self) -> &Vec<T> {
        self.0.as_ref()
    }
}

impl<T> AsRef<[T]> for VecAtLeast2<T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0.as_ref()
    }
}

impl<T> Borrow<Vec<T>> for VecAtLeast2<T> {
    #[inline]
    fn borrow(&self) -> &Vec<T> {
        &self.0
    }
}

impl<T> Borrow<[T]> for VecAtLeast2<T> {
    #[inline]
    fn borrow(&self) -> &[T] {
        self.0.borrow()
    }
}

impl<T> BorrowMut<[T]> for VecAtLeast2<T> {
    #[inline]
    fn borrow_mut(&mut self) -> &mut [T] {
        self.0.borrow_mut()
    }
}

impl<T> Deref for VecAtLeast2<T> {
    type Target = [T];

    #[inline]
    fn deref(&self) -> &Self::Target {
        self.0.deref()
    }
}

impl<T> DerefMut for VecAtLeast2<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.0.deref_mut()
    }
}

impl<T, M> Extend<M> for VecAtLeast2<T>
where
    Vec<T>: Extend<M>,
{
    #[inline]
    fn extend<I: IntoIterator<Item = M>>(&mut self, iter: I) {
        self.0.extend(iter);
    }
}

impl<T, I> Index<I> for VecAtLeast2<T>
where
    Vec<T>: Index<I>,
{
    type Output = <Vec<T> as Index<I>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        self.0.index(index)
    }
}

impl<T, I> IndexMut<I> for VecAtLeast2<T>
where
    Vec<T>: IndexMut<I>,
{
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        self.0.index_mut(index)
    }
}

impl<T> IntoIterator for VecAtLeast2<T> {
    type Item = T;
    type IntoIter = IntoIter<T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a VecAtLeast2<T> {
    type Item = &'a T;
    type IntoIter = Iter<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut VecAtLeast2<T> {
    type Item = &'a mut T;
    type IntoIter = IterMut<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<T> TryFrom<Vec<T>> for VecAtLeast2<T> {
    type Error = TooShort;

    #[inline]
    fn try_from(value: Vec<T>) -> Result<Self, Self::Error> {
        if value.is_empty() {
            Err(TooShort::Empty)
        } else if value.len() == 1 {
            Err(TooShort::Len1)
        } else {
            Ok(Self(value))
        }
    }
}

/// An error arising when attempting to construct a [`VecAtLeast2`] with
/// insufficient elements.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum TooShort {
    Empty,
    Len1,
}

impl Display for TooShort {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TooShort::Empty => write!(f, "Expected at least two elements, but none were found"),
            TooShort::Len1 => write!(f, "Expected at least two elements, but one was found"),
        }
    }
}

impl Error for TooShort {}
