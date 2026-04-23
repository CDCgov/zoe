//! ## Iterator Utilities
//!
//! This module provides miscellaneous iterators and tools for working with
//! iterators.

use std::{iter::FusedIterator, ops::ControlFlow};

#[cfg(feature = "rand")]
pub mod sampling;

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
    error: &'a mut Result<(), E>,
    iter:  Option<I>,
}

impl<I, T, E> Iterator for ProcessResults<'_, I, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        match self.iter.as_mut()?.next() {
            Some(Ok(val)) => Some(val),
            Some(Err(e)) => {
                self.iter = None;
                *self.error = Err(e);
                None
            }
            None => {
                self.iter = None;
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
                *self.error = Err(e);
                ControlFlow::Break(acc)
            }
        });

        match res {
            ControlFlow::Continue(acc) => acc,
            ControlFlow::Break(acc) => {
                self.iter = None;
                acc
            }
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
                *self.error = Err(e);
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
        match self.iter.as_mut()?.next_back() {
            Some(Ok(val)) => Some(val),
            Some(Err(e)) => {
                self.iter = None;
                *self.error = Err(e);
                None
            }
            None => {
                self.iter = None;
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
                *self.error = Err(e);
                ControlFlow::Break(acc)
            }
        });

        match res {
            ControlFlow::Continue(acc) => acc,
            ControlFlow::Break(acc) => {
                self.iter = None;
                acc
            }
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
                *self.error = Err(e);
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

impl<T, I, E> FusedIterator for ProcessResults<'_, I, E> where I: Iterator<Item = Result<T, E>> {}

/// An extension trait providing [`process_results`], a method for robustly
/// handling iterators of results in a concise and ergonomic manner.
///
/// [`process_results`]: ProcessResultsExt::process_results
pub trait ProcessResultsExt<T, E>: Iterator<Item = Result<T, E>> + Sized {
    /// Processes the `Ok` values in an iterator of results using a closure,
    /// aborting and propagating the first encountered error.
    ///
    /// This allows for the iterator of results to be treated as an iterator of
    /// values within the passed closure, with the errors handled automatically.
    ///
    /// The closure `f` accepts a [`ProcessResults`] iterator, which contains
    /// the `Ok` values up until the first error.
    ///
    /// ## Errors
    ///
    /// Any errors encountered in `iterable` are propagated. Note that this
    /// function does not guarantee that the entire iterator is checked. If this
    /// is necessary, ensure that `f` fully consumes the iterator, such as using
    /// a call to [`last`].
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
        let mut error = Ok(());

        let result = f(ProcessResults {
            error: &mut error,
            iter:  Some(self),
        });

        error.map(|()| result)
    }
}

impl<T, E, I: Iterator<Item = Result<T, E>>> ProcessResultsExt<T, E> for I {}

/// An iterator over subslices of length `window_size`, advancing forward by
/// `step` positions each iteration.
///
/// If `step < window_size`, then the subslices will be overlapping. In
/// particular, `step = 1` behaves like [`Windows`]. `step = window_size`
/// behaves like [`ChunksExact`].
///
/// Certain `step` and `window_size` combinations may not divide the input bytes
/// evenly. In this case, the final bytes are not yielded by the iterator. To
/// handle these, [`final_partial_window`] is provided, which returns a shorter
/// "partial" subslice that is formed by advancing forward again by `step` after
/// the final subslice. The return value will overlap with subslices previously
/// yielded by the iterator if `step < window_size`.
///
/// ## Parameters
///
/// - `'a`: The lifetime of the slice being iterated over
/// - `N`: The size of the windows to yield
/// - `T`: The type contained in the slice
///
/// [`Windows`]: std::slice::Windows
/// [`ChunksExact`]: std::slice::ChunksExact
/// [`final_partial_window`]: SteppedWindows::final_partial_window
/// [`ChunksExact::remainder`]: std::slice::ChunksExact::remainder
pub(crate) struct SteppedWindows<'a, const N: usize, T> {
    /// The remaining bytes to step over and yield from.
    slice:                &'a [T],
    /// The final partial window formed by advancing forward by `step` after the
    /// final complete window.
    final_partial_window: &'a [T],
    /// The amount to step forward by each iteration.
    step:                 usize,
}

impl<'a, T, const N: usize> SteppedWindows<'a, N, T> {
    /// Constructs a new [`SteppedWindows`] from `bytes`, the `step` size, and
    /// the `window_size`.
    ///
    /// ## Panics
    ///
    /// `step` and `N` must be non-zero.
    pub fn new(slice: &'a [T], step: usize) -> Self {
        assert!(step > 0, "`step` must be non-zero");
        const { assert!(N > 0, "`N` must be non-zero") };

        let final_partial_window_start = if let Some(max_window_start) = slice.len().checked_sub(N) {
            // Round down, since we cannot advance past max_window_start
            let num_steps_to_take = max_window_start / step;
            // Take 1 extra step to get the start of the remainder
            std::cmp::min((num_steps_to_take + 1) * step, slice.len())
        } else {
            0
        };

        Self {
            slice,
            final_partial_window: &slice[final_partial_window_start..],
            step,
        }
    }

    /// Returns the remainder of the original slice that is not going to be
    /// returned by the iterator. The returned slice has at most `window_size-1`
    /// elements.
    pub fn final_partial_window(&self) -> &'a [T] {
        self.final_partial_window
    }
}

impl<'a, T, const N: usize> Iterator for SteppedWindows<'a, N, T> {
    type Item = &'a [T; N];

    fn next(&mut self) -> Option<Self::Item> {
        let window = self.slice.first_chunk()?;
        self.slice = self.slice.get(self.step..).unwrap_or(&[]);
        Some(window)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = if let Some(max_window_start) = self.slice.len().checked_sub(N) {
            // max_window_start=0 implies 1 valid starting position for a window
            let num_valid_starts = max_window_start + 1;
            // Round up since we can always yield a window starting at 0. If
            // num_valid_starts=k*step, then k steps can be taken, so this also
            // holds.
            num_valid_starts.div_ceil(self.step)
        } else {
            0
        };

        (size, Some(size))
    }
}

impl<T, const N: usize> ExactSizeIterator for SteppedWindows<'_, N, T> {}
impl<T, const N: usize> FusedIterator for SteppedWindows<'_, N, T> {}
