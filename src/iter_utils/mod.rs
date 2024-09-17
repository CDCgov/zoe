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
/// ### Acknowledgements
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
