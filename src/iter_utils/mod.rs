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

/// An iterator over subslices of length `window_size`, advancing forward by `step`
/// positions each iteration.
///
/// If `step` < `window_size`, then the subslices will be overlapping. In particular,
/// `step` = 1 behaves like [`Windows`]. `step` = `window_size` behaves like
/// [`ChunksExact`].
///
/// If the number of remaining elements is less than `window_size`, these elements will
/// will be omitted. However, they can be retrieved with the [`remainder`] function from
/// the iterator.
///
/// The algorithm is based on code found at (1).
///
/// ## Acknowledgements
///
/// This solution was chiefly derived from an answer found online for
/// "[Is there a smarter way to do this? Iterators](https://www.reddit.com/r/rust/comments/8x6ep7/is_there_a_smarter_way_to_do_this_iterators/)"
/// but also took inspiration from various iterators in the standard library.
///
/// [`Windows`]: std::slice::Windows
/// [`ChunksExact`]: std::slice::ChunksExact
/// [`remainder`]: SteppedWindows::remainder
pub(crate) struct SteppedWindows<'a, T: 'a> {
    v:           &'a [T],
    rem:         &'a [T],
    step:        usize,
    window_size: usize,
}

impl<'a, T: 'a> SteppedWindows<'a, T> {
    pub fn new(slice: &'a [T], stride: usize, size: usize) -> Self {
        let rem_start = if slice.len() >= size {
            std::cmp::min(((slice.len() - size) / stride + 1) * stride, slice.len())
        } else {
            0
        };

        Self {
            v:           slice,
            rem:         &slice[rem_start..],
            step:        stride,
            window_size: size,
        }
    }

    /// Returns the remainder of the original slice that is not going to be
    /// returned by the iterator. The returned slice has at most `window_size-1`
    /// elements.
    pub fn remainder(&self) -> &'a [T] {
        self.rem
    }
}

impl<'a, T: 'a> Iterator for SteppedWindows<'a, T> {
    type Item = &'a [T];

    fn next(&mut self) -> Option<Self::Item> {
        if self.window_size <= self.v.len() {
            let subslice = &self.v[..self.window_size];
            let advance_amount = std::cmp::min(self.step, self.v.len());
            self.v = &self.v[advance_amount..];

            Some(subslice)
        } else {
            None
        }
    }
}
