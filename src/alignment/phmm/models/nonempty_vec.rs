//! A module providing a [`Vec`] wrapper with a non-empty guarantee.
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]

use crate::data::err::GetCode;
use std::{
    borrow::{Borrow, BorrowMut},
    collections::TryReserveError,
    error::Error,
    fmt::Display,
    ops::{Deref, DerefMut, Index, IndexMut, RangeBounds},
    slice::{Iter, IterMut},
    vec::IntoIter,
};

/// A wrapper around a [`Vec`] guaranteeing that it is non-empty.
///
/// This offers infallible methods for [`first`], [`last`], [`split_first`], and
/// so on.
///
/// [`first`]: NonEmptyVec::first
/// [`last`]: NonEmptyVec::last
/// [`split_first`]: NonEmptyVec::split_first
#[repr(transparent)]
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct NonEmptyVec<T>(Vec<T>);

impl<T> NonEmptyVec<T> {
    pub fn push(&mut self, value: T) {
        self.0.push(value);
    }

    pub fn push_mut(&mut self, value: T) -> &mut T {
        self.0.push_mut(value)
    }

    #[must_use]
    pub const fn capacity(&self) -> usize {
        self.0.capacity()
    }

    pub fn reserve(&mut self, additional: usize) {
        self.0.reserve(additional);
    }

    pub fn reserve_exact(&mut self, additional: usize) {
        self.0.reserve_exact(additional);
    }

    pub fn try_reserve(&mut self, additional: usize) -> Result<(), TryReserveError> {
        self.0.try_reserve(additional)
    }

    pub fn try_reserve_exact(&mut self, additional: usize) -> Result<(), TryReserveError> {
        self.0.try_reserve_exact(additional)
    }

    pub fn shrink_to_fit(&mut self) {
        self.0.shrink_to_fit();
    }

    pub fn shrink_to(&mut self, min_capacity: usize) {
        self.0.shrink_to(min_capacity);
    }

    #[must_use]
    pub fn into_boxed_slice(self) -> Box<[T]> {
        self.0.into_boxed_slice()
    }

    #[must_use]
    pub const fn as_slice(&self) -> &[T] {
        self.0.as_slice()
    }

    #[must_use]
    pub const fn as_mut_slice(&mut self) -> &mut [T] {
        self.0.as_mut_slice()
    }

    pub fn insert(&mut self, index: usize, element: T) {
        self.0.insert(index, element);
    }

    pub fn insert_mut(&mut self, index: usize, element: T) -> &mut T {
        self.0.insert_mut(index, element)
    }

    pub fn append(&mut self, other: &mut Vec<T>) {
        self.0.append(other);
    }

    #[must_use]
    pub fn first(&self) -> &T {
        self.0.first().expect("at least one element is present")
    }

    #[must_use]
    pub fn first_mut(&mut self) -> &mut T {
        self.0.first_mut().expect("at least one element is present")
    }

    #[must_use]
    pub fn split_first(&self) -> (&T, &[T]) {
        self.0.split_first().expect("at least one element is present")
    }

    #[must_use]
    pub fn split_first_mut(&mut self) -> (&mut T, &mut [T]) {
        self.0.split_first_mut().expect("at least one element is present")
    }

    #[must_use]
    pub fn split_last(&self) -> (&T, &[T]) {
        self.0.split_last().expect("at least one element is present")
    }

    #[must_use]
    pub fn split_last_mut(&mut self) -> (&mut T, &mut [T]) {
        self.0.split_last_mut().expect("at least one element is present")
    }

    #[must_use]
    pub fn last(&self) -> &T {
        self.0.last().expect("at least one element is present")
    }

    #[must_use]
    pub fn last_mut(&mut self) -> &mut T {
        self.0.last_mut().expect("at least one element is present")
    }
}

impl<T: Clone> NonEmptyVec<T> {
    pub fn extend_from_slice(&mut self, other: &[T]) {
        self.0.extend_from_slice(other);
    }

    pub fn extend_from_within<R>(&mut self, src: R)
    where
        R: RangeBounds<usize>, {
        self.0.extend_from_within(src);
    }
}

impl<T> AsMut<[T]> for NonEmptyVec<T> {
    fn as_mut(&mut self) -> &mut [T] {
        self.0.as_mut()
    }
}

impl<T> AsRef<Vec<T>> for NonEmptyVec<T> {
    fn as_ref(&self) -> &Vec<T> {
        self.0.as_ref()
    }
}

impl<T> AsRef<[T]> for NonEmptyVec<T> {
    fn as_ref(&self) -> &[T] {
        self.0.as_ref()
    }
}

impl<T> Borrow<Vec<T>> for NonEmptyVec<T> {
    fn borrow(&self) -> &Vec<T> {
        &self.0
    }
}

impl<T> Borrow<[T]> for NonEmptyVec<T> {
    fn borrow(&self) -> &[T] {
        self.0.borrow()
    }
}

impl<T> BorrowMut<[T]> for NonEmptyVec<T> {
    fn borrow_mut(&mut self) -> &mut [T] {
        self.0.borrow_mut()
    }
}

impl<T> Deref for NonEmptyVec<T> {
    type Target = [T];

    fn deref(&self) -> &Self::Target {
        self.0.deref()
    }
}

impl<T> DerefMut for NonEmptyVec<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.0.deref_mut()
    }
}

impl<T, M> Extend<M> for NonEmptyVec<T>
where
    Vec<T>: Extend<M>,
{
    fn extend<I: IntoIterator<Item = M>>(&mut self, iter: I) {
        self.0.extend(iter);
    }
}

impl<T, I> Index<I> for NonEmptyVec<T>
where
    Vec<T>: Index<I>,
{
    type Output = <Vec<T> as Index<I>>::Output;

    fn index(&self, index: I) -> &Self::Output {
        self.0.index(index)
    }
}

impl<T, I> IndexMut<I> for NonEmptyVec<T>
where
    Vec<T>: IndexMut<I>,
{
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        self.0.index_mut(index)
    }
}

impl<T> IntoIterator for NonEmptyVec<T> {
    type Item = T;
    type IntoIter = IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}
impl<'a, T> IntoIterator for &'a NonEmptyVec<T> {
    type Item = &'a T;
    type IntoIter = Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut NonEmptyVec<T> {
    type Item = &'a mut T;
    type IntoIter = IterMut<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<T> TryFrom<Vec<T>> for NonEmptyVec<T> {
    type Error = EmptyElements;

    fn try_from(value: Vec<T>) -> Result<Self, Self::Error> {
        if value.is_empty() {
            Err(EmptyElements)
        } else {
            Ok(Self(value))
        }
    }
}

/// An error arising when attempting to construct a [`NonEmptyVec`] with no
/// elements.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct EmptyElements;

impl Display for EmptyElements {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Expected at least one element, but none were found")
    }
}

impl Error for EmptyElements {}
impl GetCode for EmptyElements {}
