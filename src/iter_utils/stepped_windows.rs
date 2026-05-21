use std::iter::FusedIterator;

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
pub struct SteppedWindows<'a, const N: usize, T> {
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
    #[inline]
    #[must_use]
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
