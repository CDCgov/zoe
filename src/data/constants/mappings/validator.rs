use crate::data::array_types;
use std::ops::{Index, RangeInclusive};

/// A mapping from bytes to booleans, indicating whether the byte is valid or
/// not.
///
/// One example use-case of this struct is checking whether nucleotide or amino
/// acid sequences match a specified alphabet.
///
/// To build a custom validator, start with a [`ByteValidator::all`] (all bytes
/// are valid) or [`ByteValidator::none`] (no bytes are valid), then apply
/// cumulative overrides with methods such as [`add`], [`remove`],
/// [`add_range`], and [`remove_range`].
///
/// [`add`]: ByteValidator::add
/// [`remove`]: ByteValidator::remove
/// [`add_range`]: ByteValidator::add_range
/// [`remove_range`]: ByteValidator::remove_range
#[repr(transparent)]
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct ByteValidator([bool; 256]);

impl Default for ByteValidator {
    #[inline]
    fn default() -> Self {
        Self([false; 256])
    }
}

impl ByteValidator {
    /// Creates a new custom [`ByteValidator`] from an existing boolean array.
    #[inline]
    #[must_use]
    pub const fn new(arr: [bool; 256]) -> Self {
        Self(arr)
    }

    /// Creates a [`ByteValidator`] where all bytes are valid.
    #[inline]
    #[must_use]
    pub const fn all() -> Self {
        Self([true; 256])
    }

    /// Creates a [`ByteValidator`] where no bytes are valid.
    #[inline]
    #[must_use]
    pub const fn none() -> Self {
        Self([false; 256])
    }

    /// Marks the specified `bytes` as being valid.
    ///
    /// ## Panics
    ///
    /// The `bytes` must not contain duplicates.
    #[must_use]
    pub const fn add<const S: usize>(self, bytes: &[u8; S]) -> Self {
        assert!(
            array_types::is_unique(bytes),
            "Attempted to add a byte multiple times in the same call!"
        );

        self.set_bytes(bytes, true)
    }

    /// Marks the specified `bytes` as being invalid.
    ///
    /// ## Panics
    ///
    /// The `bytes` must not contain duplicates.
    #[must_use]
    pub const fn remove<const S: usize>(self, bytes: &[u8; S]) -> Self {
        assert!(
            array_types::is_unique(bytes),
            "Attempted to remove a byte multiple times in the same call!"
        );

        self.set_bytes(bytes, false)
    }

    /// Marks the specified `range` of bytes as being valid.
    #[must_use]
    pub const fn add_range(self, range: RangeInclusive<u8>) -> Self {
        self.set_range(range, true)
    }

    /// Marks the specified `range` of bytes as being invalid.
    #[must_use]
    pub const fn remove_range(self, range: RangeInclusive<u8>) -> Self {
        self.set_range(range, false)
    }

    /// Returns the inner array of a [`ByteValidator`].
    #[inline]
    #[must_use]
    pub const fn into_inner(self) -> [bool; 256] {
        self.0
    }

    /// Sets whether the `bytes` are valid.
    const fn set_bytes<const S: usize>(mut self, bytes: &[u8; S], valid: bool) -> Self {
        let mut i = 0;
        while i < S {
            self.set_byte(bytes[i], valid);
            i += 1;
        }
        self
    }

    /// Sets whether the bytes in `range` are valid.
    #[allow(clippy::cast_possible_truncation)]
    const fn set_range(mut self, range: RangeInclusive<u8>, valid: bool) -> Self {
        // Convert from inclusive range to exclusive range, bumping to usize in
        // case of overflow at end
        let from = *range.start() as usize..*range.end() as usize + 1;

        let mut byte = from.start;
        while byte < from.end {
            // Validity: this will not truncate since src is strictly less than
            // the original from.end()+1, which means it is at most from.end(),
            // which is a `u8`.
            self.set_byte(byte as u8, valid);

            byte += 1;
        }

        self
    }

    /// Sets whether a byte is valid.
    #[inline]
    pub(crate) const fn set_byte(&mut self, byte: u8, valid: bool) {
        self.0[byte as usize] = valid;
    }
}

impl Index<u8> for ByteValidator {
    type Output = bool;

    #[inline]
    fn index(&self, index: u8) -> &bool {
        &self.0[index as usize]
    }
}
