use crate::alignment::{Alignment, AlignmentStates};
use std::{
    mem::MaybeUninit,
    simd::{LaneCount, SupportedLaneCount, prelude::*},
};

/// Experimental backtracking matrix for alignment with affine gap scores.
#[derive(Clone, Debug)]
pub(crate) struct BacktrackMatrix {
    pub data: Vec<u8>,
    cols:     usize,
    cursor:   usize,
}

impl BacktrackMatrix {
    const UP: u8 = 1;
    const UP_EXTENDING: u8 = 2;
    const LEFT: u8 = 4;
    const LEFT_EXTENDING: u8 = 8;
    const STOP: u8 = 16;

    pub(crate) fn new(rows: usize, cols: usize) -> Self {
        BacktrackMatrix {
            data: vec![0u8; rows * cols],
            cols,
            cursor: 0,
        }
    }

    #[inline]
    pub(crate) fn up(&mut self) {
        self.data[self.cursor] |= Self::UP;
    }

    #[inline]
    pub(crate) fn left(&mut self) {
        self.data[self.cursor] |= Self::LEFT;
    }

    #[inline]
    pub(crate) fn up_extending(&mut self) {
        self.data[self.cursor] |= Self::UP_EXTENDING;
    }

    #[inline]
    pub(crate) fn left_extending(&mut self) {
        self.data[self.cursor] |= Self::LEFT_EXTENDING;
    }

    #[inline]
    pub(crate) fn stop(&mut self) {
        self.data[self.cursor] = Self::STOP;
    }
}

#[derive(Clone, Debug)]
pub(crate) struct BacktrackMatrixStriped<const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    pub data:  Vec<Simd<u8, N>>,
    num_vecs:  usize,
    v_cursor:  usize,
    curr_lane: usize,
}

impl<const N: usize> BacktrackMatrixStriped<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// Allocates lazily-initialized backtrack data which can later be converted
    /// to a [`BacktrackMatrixStriped`].
    ///
    /// The total number of SIMD vectors is given by `size`.
    #[inline]
    #[must_use]
    pub(crate) fn make_uninit_data(size: usize) -> Vec<MaybeUninit<Simd<u8, N>>> {
        let mut data = Vec::with_capacity(size);
        data.resize_with(size, MaybeUninit::uninit);
        data
    }

    /// Wraps a vector of SIMD vectors in a [`BacktrackMatrixStriped`] to
    /// facilitate backtracking.
    #[inline]
    #[must_use]
    pub(crate) fn new(data: Vec<Simd<u8, N>>, num_vecs: usize) -> Self {
        assert!(data.len().is_multiple_of(num_vecs));
        BacktrackMatrixStriped {
            data,
            num_vecs,
            v_cursor: 0,
            curr_lane: 0,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn debug_cell(&self, r: usize, c: usize) -> String {
        let v = c % self.num_vecs;
        let lane = (c - v) / self.num_vecs;
        let v = self.data[self.num_vecs * r + v];
        v.debug_cell(lane)
    }

    #[allow(dead_code)]
    pub(crate) fn print_row(&self, r: usize) {
        let mut bt = self.clone();
        let nc = self.num_vecs * N;

        print!("{r:02}: x");
        for c in 0..nc {
            bt.move_to(r, c);
            if bt.is_stop() {
                print!("o");
            } else if bt.is_up() {
                print!("^");
            } else if bt.is_left() {
                print!("<");
            } else if bt.is_up_extending() {
                print!(":");
            } else if bt.is_left_extending() {
                print!("-");
            } else {
                print!("\\");
            }
        }
        println!();
    }
}

/// Trait for SIMD backtracking operations on packed u8 values
pub(crate) trait SimdBacktrackFlags<const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    const UP: Simd<u8, N>;
    const UP_EXTENDING: Simd<u8, N>;
    const LEFT: Simd<u8, N>;
    const LEFT_EXTENDING: Simd<u8, N>;
    const STOP: Simd<u8, N>;

    /// Creates a new match variable
    fn simd_match() -> Simd<u8, N> {
        Simd::splat(0)
    }

    /// Add Up direction to lanes selected by mask
    fn simd_up(&mut self, mask: Mask<i8, N>);

    /// Add Left direction to lanes selected by mask
    fn simd_left(&mut self, mask: Mask<i8, N>);

    fn simd_correct_and_set_left(&mut self, mask: Mask<i8, N>);

    /// Apply Up Extending direction to lanes selected by mask
    fn simd_up_extending(&mut self, mask: Mask<i8, N>);

    /// Add Left Extending direction to lanes selected by mask
    fn simd_left_extending(&mut self, mask: Mask<i8, N>);

    /// Set Stop state to lanes selected by mask
    fn simd_stop(&mut self, mask: Mask<i8, N>);

    /// Debug representation of a cell at the given lane
    fn debug_cell(&self, lane: usize) -> String;
}

impl<const N: usize> SimdBacktrackFlags<N> for Simd<u8, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    const UP: Simd<u8, N> = Simd::splat(BacktrackMatrix::UP);
    const UP_EXTENDING: Simd<u8, N> = Simd::splat(BacktrackMatrix::UP_EXTENDING);
    const LEFT: Simd<u8, N> = Simd::splat(BacktrackMatrix::LEFT);
    const LEFT_EXTENDING: Simd<u8, N> = Simd::splat(BacktrackMatrix::LEFT_EXTENDING);
    const STOP: Simd<u8, N> = Simd::splat(BacktrackMatrix::STOP);

    #[inline]
    fn simd_up(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(*self | Self::UP, *self);
    }

    #[inline]
    fn simd_left(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(*self | Self::LEFT, *self);
    }

    #[inline]
    fn simd_correct_and_set_left(&mut self, mask: Mask<i8, N>) {
        *self = mask.select((*self & Self::UP_EXTENDING) | Self::LEFT, *self);
    }

    #[inline]
    fn simd_up_extending(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(*self | Self::UP_EXTENDING, *self);
    }

    #[inline]
    fn simd_left_extending(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(*self | Self::LEFT_EXTENDING, *self);
    }

    #[inline]
    fn simd_stop(&mut self, mask: Mask<i8, N>) {
        *self = mask.select(Self::STOP, *self);
    }

    fn debug_cell(&self, lane: usize) -> String {
        let cell = self[lane];

        let mut states = String::new();
        if cell & BacktrackMatrix::STOP > 0 {
            states.push('o');
        }

        if cell & BacktrackMatrix::UP > 0 {
            states.push('^');
        }

        if cell & BacktrackMatrix::LEFT > 0 {
            states.push('<');
        }

        if cell & BacktrackMatrix::UP_EXTENDING > 0 {
            states.push(':');
        }

        if cell & BacktrackMatrix::LEFT_EXTENDING > 0 {
            states.push('-');
        }

        if cell == 0 {
            states.push('\\');
        }

        states
    }
}

/// Trait for backtracking operations on a backtracking matrix
pub(crate) trait BackTrackable {
    /// Moves the cursor to the specified row and column
    fn move_to(&mut self, r: usize, c: usize);

    /// Checks if the current cell is in the **up** state
    fn is_up(&self) -> bool;

    /// Checks if the current cell is in the **left** state
    fn is_left(&self) -> bool;

    /// Checks if the current cell is in the **left extending** state
    fn is_left_extending(&self) -> bool;

    /// Checks if the current cell is in the **up extending** state
    fn is_up_extending(&self) -> bool;

    /// Checks if the current cell is in the **stop** state
    fn is_stop(&self) -> bool;

    /// Creates an [`Alignment`] from a backtrack matrix and other information.
    fn to_alignment<T>(
        &mut self, score: T, mut r_end: usize, mut c_end: usize, ref_len: usize, query_len: usize,
    ) -> Alignment<T> {
        let mut states = AlignmentStates::new();
        let mut op = 0;

        self.move_to(r_end, c_end); // 0-based move to max

        // 1-based as though we had a padded matrix
        r_end += 1;
        c_end += 1;

        let (mut r, mut c) = (r_end, c_end);

        // soft clip 3'
        states.soft_clip(query_len - c);

        while !self.is_stop() && r > 0 && c > 0 {
            if op == b'D' && self.is_up_extending() {
                op = b'D';
                r -= 1;
            } else if op == b'I' && self.is_left_extending() {
                op = b'I';
                c -= 1;
            } else if self.is_up() {
                op = b'D';
                r -= 1;
            } else if self.is_left() {
                op = b'I';
                c -= 1;
            } else {
                op = b'M';
                r -= 1;
                c -= 1;
            }
            states.add_state(op);
            self.move_to(r.saturating_sub(1), c.saturating_sub(1)); // 0-based next position
        }

        // soft clip 5'
        states.soft_clip(c);
        states.make_reverse();

        // r and c are decremented and becomes 0-based
        Alignment {
            score,
            ref_range: r..r_end,
            query_range: c..c_end,
            states,
            ref_len,
            query_len,
        }
    }

    /// Prints the backtracking matrix in a human-readable format
    #[allow(dead_code)]
    fn print(&self);
}

impl BackTrackable for BacktrackMatrix {
    #[inline]
    fn move_to(&mut self, r: usize, c: usize) {
        self.cursor = self.cols * r + c;
    }

    #[inline]
    fn is_up(&self) -> bool {
        self.data[self.cursor] & Self::UP > 0
    }

    #[inline]
    fn is_left(&self) -> bool {
        self.data[self.cursor] & Self::LEFT > 0
    }

    #[inline]
    fn is_left_extending(&self) -> bool {
        self.data[self.cursor] & Self::LEFT_EXTENDING > 0
    }

    #[inline]
    fn is_up_extending(&self) -> bool {
        self.data[self.cursor] & Self::UP_EXTENDING > 0
    }

    #[inline]
    fn is_stop(&self) -> bool {
        self.data[self.cursor] & Self::STOP > 0
    }

    fn print(&self) {
        let mut bt = self.clone();
        let nr = self.data.len() / self.cols;
        let nc = self.cols;

        print!("    ");
        for _ in 0..=nc {
            print!("x");
        }
        println!();
        for r in 0..nr {
            print!("{r:02}: x");
            for c in 0..nc {
                bt.move_to(r, c);
                if bt.is_stop() {
                    print!("o");
                } else if bt.is_up() {
                    print!("^");
                } else if bt.is_left() {
                    print!("<");
                } else if bt.is_up_extending() {
                    print!(":");
                } else if bt.is_left_extending() {
                    print!("-");
                } else {
                    print!("\\");
                }
            }
            println!();
        }
    }
}

impl<const N: usize> BackTrackable for BacktrackMatrixStriped<N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn move_to(&mut self, r: usize, c: usize) {
        let v = c % self.num_vecs;
        self.curr_lane = (c - v) / self.num_vecs;
        self.v_cursor = self.num_vecs * r + v;
    }

    #[inline]
    fn is_up(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::UP) > 0
    }

    #[inline]
    fn is_left(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::LEFT) > 0
    }

    #[inline]
    fn is_left_extending(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::LEFT_EXTENDING) > 0
    }

    #[inline]
    fn is_up_extending(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::UP_EXTENDING) > 0
    }

    #[inline]
    fn is_stop(&self) -> bool {
        (self.data[self.v_cursor][self.curr_lane] & BacktrackMatrix::STOP) > 0
    }

    fn print(&self) {
        let mut bt = self.clone();
        let nr = self.data.len() / self.num_vecs;
        let nc = self.num_vecs * N;

        print!("    ");
        for _ in 0..=nc {
            print!("x");
        }
        println!();
        for r in 0..nr {
            print!("{r:02}: x");
            for c in 0..nc {
                bt.move_to(r, c);
                if bt.is_stop() {
                    print!("o");
                } else if bt.is_up() {
                    print!("^");
                } else if bt.is_left() {
                    print!("<");
                } else if bt.is_up_extending() {
                    print!(":");
                } else if bt.is_left_extending() {
                    print!("-");
                } else {
                    print!("\\");
                }
            }
            println!();
        }
    }
}
