//! ## Iterator Utilities
//!
//! This module provides miscellaneous iterators and tools for working with
//! iterators.

use std::{cell::OnceCell, iter::FusedIterator, ops::ControlFlow};

#[cfg(feature = "rand")]
pub mod sampling;

#[cfg(not(feature = "fuzzing"))]
#[doc(auto_cfg(hide(feature, values("fuzzing"))))]
mod stepped_windows;
#[cfg(not(feature = "fuzzing"))]
#[doc(auto_cfg(hide(feature, values("fuzzing"))))]
pub(crate) use stepped_windows::*;

#[cfg(feature = "fuzzing")]
mod stepped_windows;
#[cfg(feature = "fuzzing")]
pub use stepped_windows::*;

/// An iterator that extracts the `Ok` variants from an input iterator of
/// results, ending upon the first encountered `Err` variant and storing that
/// error.
///
/// ## Acknowledgements
///
/// This is inspired by similar iterators in
/// [Itertools](https://docs.rs/itertools/latest/itertools/) and
/// [`iterr`](https://docs.rs/iterr/latest/iterr/).
#[derive(Debug)]
pub struct ProcessResults<'a, I, E: 'a> {
    /// The first encountered error, if any.
    error: &'a OnceCell<E>,
    /// The fallible iterator, or `None` if an error has occurred and been
    /// stored.
    iter:  Option<I>,
}

impl<I, T, E> Iterator for ProcessResults<'_, I, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        // Return None if an error already occurred or the fallible iterator
        // returns None
        match self.iter.as_mut()?.next()? {
            Ok(val) => Some(val),
            Err(e) => {
                self.iter = None;
                let _ = self.error.set(e);
                None
            }
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        match &self.iter {
            Some(iter) => (0, iter.size_hint().1),
            None => (0, Some(0)),
        }
    }

    fn fold<B, F>(mut self, init: B, mut f: F) -> B
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> B, {
        let Some(iter) = &mut self.iter else {
            return init;
        };

        let res = iter.try_fold(init, |acc, res| match res {
            Ok(val) => ControlFlow::Continue(f(acc, val)),
            Err(e) => {
                let _ = self.error.set(e);
                ControlFlow::Break(acc)
            }
        });

        match res {
            ControlFlow::Continue(acc) | ControlFlow::Break(acc) => acc,
        }
    }

    fn try_fold<B, F, R>(&mut self, init: B, mut f: F) -> R
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> R,
        R: std::ops::Try<Output = B>, {
        let Some(iter) = &mut self.iter else {
            return R::from_output(init);
        };

        // Call the inner try_fold impl
        let res = iter.try_fold(init, |acc, res| match res {
            // The incoming result was `Ok`, so call the closure, breaking if an
            // error occurs. Any such error is wrapped in Err to indicate that
            // the outer try_fold should return an error
            Ok(val) => f(acc, val).branch().map_break(Err),
            // The incoming result was `Err`, so we must store this and abort
            // the iterator (as would happen with a next call). We do this by
            // issuing a Break, but we wrap the value in Ok to indicate that the
            // outer `try_fold` should not return an error
            Err(e) => {
                let _ = self.error.set(e);
                ControlFlow::Break(Ok(acc))
            }
        });

        match res {
            // No errors in iterator, no errors due to closure
            ControlFlow::Continue(val) => R::from_output(val),
            // Error encountered in iterator, not due to closure, so return Ok
            ControlFlow::Break(Ok(acc)) => {
                self.iter = None;
                R::from_output(acc)
            }
            // Error encountered in closure, not due to iterator, so return Err
            ControlFlow::Break(Err(err)) => R::from_residual(err),
        }
    }
}

impl<I, T, E> DoubleEndedIterator for ProcessResults<'_, I, E>
where
    I: DoubleEndedIterator<Item = Result<T, E>>,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        // Return None if an error already occurred or the fallible iterator
        // returns None
        match self.iter.as_mut()?.next_back()? {
            Ok(val) => Some(val),
            Err(e) => {
                self.iter = None;
                let _ = self.error.set(e);
                None
            }
        }
    }

    fn rfold<B, F>(mut self, init: B, mut f: F) -> B
    where
        F: FnMut(B, Self::Item) -> B, {
        let Some(iter) = &mut self.iter else {
            return init;
        };

        let res = iter.try_rfold(init, |acc, res| match res {
            Ok(val) => ControlFlow::Continue(f(acc, val)),
            Err(e) => {
                let _ = self.error.set(e);
                ControlFlow::Break(acc)
            }
        });

        match res {
            ControlFlow::Continue(acc) | ControlFlow::Break(acc) => acc,
        }
    }

    fn try_rfold<B, F, R>(&mut self, init: B, mut f: F) -> R
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> R,
        R: std::ops::Try<Output = B>, {
        let Some(iter) = &mut self.iter else {
            return R::from_output(init);
        };

        // Call the inner try_fold impl
        let res = iter.try_rfold(init, |acc, res| match res {
            // The incoming result was `Ok`, so call the closure, breaking if an
            // error occurs. Any such error is wrapped in Err to indicate that
            // the outer try_fold should return an error
            Ok(val) => f(acc, val).branch().map_break(Err),
            // The incoming result was `Err`, so we must store this and abort
            // the iterator (as would happen with a next call). We do this by
            // issuing a Break, but we wrap the value in Ok to indicate that the
            // outer `try_fold` should not return an error
            Err(e) => {
                let _ = self.error.set(e);
                ControlFlow::Break(Ok(acc))
            }
        });

        match res {
            // No errors in iterator, no errors due to closure
            ControlFlow::Continue(val) => R::from_output(val),
            // Error encountered in iterator, not due to closure, so return Ok
            ControlFlow::Break(Ok(acc)) => {
                self.iter = None;
                R::from_output(acc)
            }
            // Error encountered in closure, not due to iterator, so return Err
            ControlFlow::Break(Err(err)) => R::from_residual(err),
        }
    }
}

impl<T, I, E> FusedIterator for ProcessResults<'_, I, E> where I: FusedIterator<Item = Result<T, E>> {}

/// An extension trait providing [`process_results`], a method for robustly
/// handling iterators of results in a concise and ergonomic manner.
///
/// [`process_results`]: ProcessResultsExt::process_results
pub trait ProcessResultsExt<T, E>: Iterator<Item = Result<T, E>> + Sized {
    /// Processes the `Ok` values in an iterator of results using a closure,
    /// aborting upon and propagating the first encountered error.
    ///
    /// This allows for the iterator of results to be treated as an iterator of
    /// values within the passed closure, with the errors handled automatically.
    ///
    /// The closure `f` accepts a [`ProcessResults`] iterator, which contains
    /// the `Ok` values up until the first error.
    ///
    /// ## Errors
    ///
    /// Any errors encountered in `self` are propagated. Note that this function
    /// does not guarantee that the entire iterator is checked. If this is
    /// necessary, ensure that `f` fully consumes the iterator, such as using a
    /// call to [`last`] in the case that `self` is fused.
    ///
    /// ## Acknowledgements
    ///
    /// This is inspired by similar iterators in
    /// [Itertools](https://docs.rs/itertools/latest/itertools/) and
    /// [`iterr`](https://docs.rs/iterr/latest/iterr/).
    ///
    /// [`last`]: Iterator::last
    fn process_results<F, R>(self, f: F) -> Result<R, E>
    where
        F: FnOnce(ProcessResults<Self, E>) -> R, {
        let error = OnceCell::new();

        let value = f(ProcessResults {
            error: &error,
            iter:  Some(self),
        });

        match error.into_inner() {
            Some(err) => Err(err),
            None => Ok(value),
        }
    }

    /// A method called within [`process_results_many`] to turn a fallible
    /// iterator into an infallible iterator.
    fn or_stop(self, cx: &FallibleContext<E>) -> ProcessResults<'_, Self, E> {
        ProcessResults {
            error: &cx.error,
            iter:  Some(self),
        }
    }
}

impl<T, E, I: Iterator<Item = Result<T, E>>> ProcessResultsExt<T, E> for I {}

/// A method for handling multiple fallible iterators (or multiple fallible
/// steps in a single iterator) at once.
///
/// The [`process_results`] method is a good tool for handling fallible
/// iterators, using a closure to contain the logic to perform on the non-error
/// values of the iterator, and then automatically handling error propagation.
/// However, when there are many fallible iterators or many steps of
/// fallibility, using [`process_results`] can cause multiple nested closures,
/// reducing readability.
///
/// This standalone function provides a context allowing any fallible iterator
/// to be transformed into an infallible iterator by calling the [`or_stop`]
/// method on it. This methods yields the `Ok` items from the iterator until an
/// error is reached, after which `None` is returned. The first encountered
/// error is stored to the context and propagated when the closure exits.
///
/// ## Examples
///
/// Below is an example where two fallible iterators are zipped together:
///
/// ```
/// # use std::array::IntoIter;
/// # use zoe::iter_utils::{ProcessResultsExt, process_results_many};
///
/// let iter1 = [Ok(1), Ok(2), Ok(3), Ok(4)].into_iter();
/// let iter2 = [Ok('A'), Ok('B'), Ok('C')].into_iter();
/// # let iter1: IntoIter<Result<i32, ()>, 4> = iter1;
///
/// let zipped = process_results_many(|cx| {
///     let iter1 = iter1.or_stop(cx);
///     let iter2 = iter2.or_stop(cx);
///     iter1.zip(iter2).collect::<Vec<_>>()
/// });
///
/// assert_eq!(zipped, Ok(vec![(1, 'A'), (2, 'B'), (3, 'C')]));
///
///
/// let iter1 = [Ok(1), Ok(2), Ok(3), Ok(4)].into_iter();
/// let iter2 = [Ok('A'), Err("Failure"), Ok('C')].into_iter();
/// let zipped = process_results_many(|cx| {
///     let iter1 = iter1.or_stop(cx);
///     let iter2 = iter2.or_stop(cx);
///     iter1.zip(iter2).collect::<Vec<_>>()
/// });
///
/// assert_eq!(zipped, Err("Failure"));
/// ```
///
/// Another example involves two fallible `map` operations on a single iterator:
///
/// ```
/// # use zoe::iter_utils::{ProcessResultsExt, process_results_many};
///
/// let data = [1, 2, 3];
/// let idx_iter = ["1", "2", "3", "D"].into_iter();
///
/// let vals = process_results_many(|cx| {
///     idx_iter
///         .map(|s| s.parse::<usize>().map_err(|_| "Failed to parse index"))
///         .or_stop(cx)
///         .map(|idx| data.get(idx).ok_or("Index out of bounds"))
///         .or_stop(cx)
///         .collect::<Vec<_>>()
/// });
///
/// assert_eq!(vals, Err("Index out of bounds"));
/// ```
///
/// ## Limitations
///
/// Users must ensure that one of the fallible iterators returning `None` is
/// sufficient to cause the closure to exit, without unanticipated work or
/// side-effects. For example, using `zip_eq` from Itertools could cause
/// erroneous panics, or `zip_longest` could cause it to appear that the two
/// iterators are different lengths when in reality one encountered an error.
///
/// Furthermore, only the first error is propagated, despite any number of other
/// errors being potentially present in the iterator. If some of the errors in
/// the iterator get ignored or handled by the application, then it may be
/// necessary to perform this logic within the closure rather than waiting until
/// after `process_results_many`. Otherwise, some errors may be silently dropped
/// or never reached.
///
/// ## Errors
///
/// The first error returned by any of the iterators gets propagated, once the
/// closure completes.
///
/// [`process_results`]: ProcessResultsExt::process_results
/// [`or_stop`]: ProcessResultsExt::or_stop
pub fn process_results_many<F, T, E>(f: F) -> Result<T, E>
where
    F: for<'a> FnOnce(&'a FallibleContext<E>) -> T, {
    let cx = FallibleContext { error: OnceCell::new() };

    let value = f(&cx);

    match cx.error.into_inner() {
        Some(err) => Err(err),
        None => Ok(value),
    }
}

/// A context for use within [`process_results_many`], able to convert fallible
/// iterators into infallible iterators, recording and propagating the first
/// encountered error.
pub struct FallibleContext<E> {
    error: OnceCell<E>,
}
