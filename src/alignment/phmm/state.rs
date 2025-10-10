use crate::alignment::phmm::{PhmmError, PhmmNumber};

/// An enum representing the three states within each layer of a pHMM.
///
/// This is used for readability when indexing.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// access this enum.
///
/// </div>
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhmmState {
    Match  = 0,
    Delete = 1,
    Insert = 2,
}

/// An enum representing the three states within each layer of a pHMM, in
/// addition to `Enter`.
///
/// This is useful for local pHMMs.
///
/// <div class="warning note">
///
/// **Note**
///
/// You must enable the *alignment-diagnostics* feature in your `Cargo.toml` to
/// access this enum.
///
/// </div>
#[repr(u8)]
#[derive(Clone, Copy, Debug)]
pub enum PhmmStateOrEnter {
    Match  = 0,
    Delete = 1,
    Insert = 2,
    Enter  = 3,
}

impl From<PhmmState> for u8 {
    #[inline]
    fn from(value: PhmmState) -> Self {
        value as u8
    }
}

impl From<PhmmStateOrEnter> for u8 {
    #[inline]
    fn from(value: PhmmStateOrEnter) -> Self {
        value as u8
    }
}

impl From<PhmmState> for usize {
    #[inline]
    fn from(value: PhmmState) -> Self {
        value as usize
    }
}

impl From<PhmmStateOrEnter> for usize {
    #[inline]
    fn from(value: PhmmStateOrEnter) -> Self {
        value as usize
    }
}

impl From<PhmmState> for PhmmStateOrEnter {
    #[inline]
    fn from(value: PhmmState) -> Self {
        match value {
            PhmmState::Match => PhmmStateOrEnter::Match,
            PhmmState::Delete => PhmmStateOrEnter::Delete,
            PhmmState::Insert => PhmmStateOrEnter::Insert,
        }
    }
}

impl PhmmState {
    /// Gets a [`PhmmState`] from a [`PhmmStateOrEnter`], returning `None` for
    /// the [`Enter`] state.
    ///
    /// [`Enter`]: PhmmStateOrEnter::Enter
    #[inline]
    #[must_use]
    pub(crate) fn get_from(value: PhmmStateOrEnter) -> Option<PhmmState> {
        match value {
            PhmmStateOrEnter::Match => Some(PhmmState::Match),
            PhmmStateOrEnter::Delete => Some(PhmmState::Delete),
            PhmmStateOrEnter::Insert => Some(PhmmState::Insert),
            PhmmStateOrEnter::Enter => None,
        }
    }

    /// Converts a CIGAR-style operation to a [`PhmmState`].
    ///
    /// ## Errors
    ///
    /// The operation must be in `MDI=X`.
    #[inline]
    #[allow(dead_code)]
    pub(crate) fn from_op(op: u8) -> Result<Self, PhmmError> {
        match op {
            b'M' | b'=' | b'X' => Ok(PhmmState::Match),
            b'D' => Ok(PhmmState::Delete),
            b'I' => Ok(PhmmState::Insert),
            _ => Err(PhmmError::InvalidCigarOp),
        }
    }
}

/// A [`PhmmState`] or [`PhmmStateOrEnter`] represented as a u8.
///
/// `Match` corresponds to 0, `Delete` corresponds to 1, `Insert` corresponds to
/// 2, and `Enter` corresponds to 3.
#[repr(transparent)]
#[derive(Copy, Clone)]
pub(crate) struct PhmmTracebackState(u8);

impl PhmmTracebackState {
    /// Returns whether the stored state is `Match`.
    #[inline]
    pub fn is_match(self) -> bool {
        self.0 == 0
    }

    /// Returns whether the stored state is `Delete`.
    #[inline]
    pub fn is_delete(self) -> bool {
        self.0 == 1
    }

    /// Returns whether the stored state is `Insert`.
    #[inline]
    #[allow(dead_code)]
    pub fn is_insert(self) -> bool {
        self.0 == 2
    }

    /// Returns whether the stored state is `Enter`.
    #[inline]
    pub fn is_enter(self) -> bool {
        self.0 == 3
    }
}

impl From<PhmmState> for PhmmTracebackState {
    #[inline]
    fn from(value: PhmmState) -> Self {
        PhmmTracebackState(value as u8)
    }
}

impl From<PhmmStateOrEnter> for PhmmTracebackState {
    #[inline]
    fn from(value: PhmmStateOrEnter) -> Self {
        PhmmTracebackState(value as u8)
    }
}

/// A set of flags for storing the traceback information for a particular layer
/// in a pHMM.
///
/// Three values are stored: the previous state used to reach the match state,
/// the previous state used to reach the insert state, and the previous state
/// used to reach the delete state. Specifically, this is three
/// [`PhmmTracebackState`] values packed into a single `u8`.
#[repr(transparent)]
#[derive(Copy, Clone)]
pub(crate) struct PhmmBacktrackFlags(u8);

impl PhmmBacktrackFlags {
    /// Creates a new [`PhmmBacktrackFlags`] object with all the previous states
    /// set to `Match`.
    pub fn new() -> Self {
        Self(0)
    }

    /// Sets the previous state for `Match` to `prev_state`.
    ///
    /// This function can only be called once on a given [`PhmmBacktrackFlags`],
    /// otherwise erroneous behavior could occur.
    #[inline]
    #[allow(clippy::verbose_bit_mask)]
    pub fn set_match(&mut self, prev_state: impl Into<PhmmTracebackState>) {
        let prev_state = prev_state.into();

        debug_assert!(prev_state.0 < 4);
        debug_assert!(self.0 & 0b00_00_00_11 == 0);

        self.0 |= prev_state.0;
    }

    /// Sets the previous state for `Delete` to `prev_state`.
    ///
    /// This function can only be called once on a given [`PhmmBacktrackFlags`],
    /// otherwise erroneous behavior could occur.
    #[inline]
    pub fn set_delete(&mut self, prev_state: impl Into<PhmmTracebackState>) {
        let prev_state = prev_state.into();

        debug_assert!(prev_state.0 < 4);
        debug_assert!(self.0 & 0b00_00_11_00 == 0);

        self.0 |= prev_state.0.wrapping_shl(2);
    }

    /// Sets the previous state for `Insert` to `prev_state`.
    ///
    /// This function can only be called once on a given [`PhmmBacktrackFlags`],
    /// otherwise erroneous behavior could occur.
    #[inline]
    pub fn set_insert(&mut self, prev_state: impl Into<PhmmTracebackState>) {
        let prev_state = prev_state.into();

        debug_assert!(prev_state.0 < 4);
        debug_assert!(self.0 & 0b00_11_00_00 == 0);

        self.0 |= prev_state.0.wrapping_shl(4);
    }

    /// Gets the stored previous state for the given `current_state`.
    #[inline]
    #[must_use]
    pub fn get_prev_state(self, current_state: impl Into<PhmmTracebackState>) -> PhmmTracebackState {
        let current_state = current_state.into();

        // current_state is 0, 1, 2, or 3
        debug_assert!(current_state.0 < 4);
        // shift is 0, 2, 4, or 6 respectively
        let shift = current_state.0.wrapping_mul(2);
        // mask is 0b00000011, 0b00001100, 0b00110000, 0b11000000 respectively
        let mask = 1u8.wrapping_shl(u32::from(shift)).wrapping_mul(3);
        // Given self.flags as 0bABCDEFGH, the result is 0b000000GH, 0b000000EF,
        // 0b000000CD, 0b000000AB respectively.
        PhmmTracebackState((self.0 & mask) >> shift)
    }
}

/// Identifies the minimum score among three values (corresponding to `Match`,
/// `Delete`, and `Insert`) and returns the best state and score.
pub(crate) fn best_state<T: PhmmNumber>(match_val: T, delete_val: T, insert_val: T) -> (PhmmState, T) {
    use PhmmState::*;

    let mut argmin = Match;
    let mut min = match_val;

    for (state, val) in [(Delete, delete_val), (Insert, insert_val)] {
        if val < min {
            argmin = state;
            min = val;
        }
    }

    (argmin, min)
}

/// Identifies the minimum score among four values (corresponding to `Match`,
/// `Delete`, and `Insert`, and `Enter`) and returns the best state and score.
pub(crate) fn best_state_or_enter<T: PhmmNumber>(
    match_val: T, delete_val: T, insert_val: T, enter_val: T,
) -> (PhmmStateOrEnter, T) {
    use PhmmStateOrEnter::*;

    let mut argmin = Match;
    let mut min = match_val;

    for (state, val) in [(Delete, delete_val), (Insert, insert_val), (Enter, enter_val)] {
        if val < min {
            argmin = state;
            min = val;
        }
    }

    (argmin, min)
}
