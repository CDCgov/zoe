use crate::data::array_types::{self, arr_max, position};
use std::ops::Index;

/// Represents a mapping between bytes and indices. For example, this could be a
/// map from DNA bases to profile indices, such as [`DNA_PROFILE_MAP`].
///
/// ## Type Parameters
/// * `KEYS` - The number of bytes being mapped (such as 5 for DNA including
///   *N*)
///
/// [`DNA_PROFILE_MAP`]: super::DNA_PROFILE_MAP
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct ByteIndexMap<const KEYS: usize> {
    pub(crate) index_map: [u8; 256],
    pub(crate) byte_keys: [u8; KEYS],
}

impl<const S: usize> ByteIndexMap<S> {
    /// Create a new [`ByteIndexMap`] struct to represent a mapping between
    /// bytes and indices. For example, this could be a map from DNA bases to
    /// profile indices. Any byte that is not specified in `byte_keys` is mapped
    /// to the same thing as `catch_all`.
    ///
    /// ## Panics
    /// No duplicates can be present in `byte_keys`. `catch_all` must be present
    /// in `byte_keys`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new(byte_keys: [u8; S], catch_all: u8) -> Self {
        assert!(array_types::is_unique(&byte_keys));

        let catch_all_index =
            position(&byte_keys, catch_all).expect("The catch_all must be present in the byte_keys.") as u8;
        let mut out = ByteIndexMap {
            index_map: [catch_all_index; 256],
            byte_keys,
        };

        let mut i = 0;
        while i < byte_keys.len() {
            // Truncation will not occur because i cannot exceed
            // byte_keys.len(), and index must contain unique u8 values
            out.set_byte(byte_keys[i], i as u8);

            i += 1;
        }
        out
    }

    /// Create a new [`ByteIndexMap`] struct to represent a mapping between
    /// bytes and indices. For example, this could be a map from DNA bases to
    /// profile indices. Both `byte_keys` and `catch_all` ignore case.
    ///
    /// ## Panics
    /// No duplicates can be present in `byte_keys`. `catch_all` must be present
    /// in `byte_keys`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new_ignoring_case(mut byte_keys: [u8; S], catch_all: u8) -> Self {
        byte_keys = array_types::make_uppercase(&byte_keys);
        assert!(array_types::is_unique(&byte_keys));

        let catch_all_index = position(&byte_keys, catch_all.to_ascii_uppercase())
            .expect("The catch_all must be present in the byte_keys.") as u8;
        let mut out = ByteIndexMap {
            index_map: [catch_all_index; 256],
            byte_keys,
        };

        let mut i = 0;
        while i < byte_keys.len() {
            // Truncation will not occur because i cannot exceed
            // byte_keys.len(), and byte_keys must contain unique u8 values
            out.set_byte_ignoring_case(byte_keys[i], i as u8);

            i += 1;
        }
        out
    }

    /// Increment all values in the map to support a different starting index.
    pub(crate) const fn update_starting_index(mut self, starting_index: u8) -> Self {
        assert!(arr_max(&self.index_map).unwrap() < u8::MAX - starting_index);
        let mut i = 0;
        while i < self.index_map.len() {
            self.index_map[i] += starting_index;
            i += 1;
        }
        self
    }

    /// Set the index for a byte.
    #[inline]
    const fn set_byte(&mut self, byte: u8, index: u8) {
        self.index_map[byte as usize] = index;
    }

    /// Set the index for a byte, ignoring case.
    #[inline]
    const fn set_byte_ignoring_case(&mut self, byte: u8, index: u8) {
        self.index_map[byte.to_ascii_lowercase() as usize] = index;
        self.index_map[byte.to_ascii_uppercase() as usize] = index;
    }

    /// Change the [`ByteIndexMap`] so that `new_byte` maps to the same thing as
    /// `byte_key`.
    #[inline]
    #[must_use]
    pub const fn add_synonym(mut self, new_key: u8, previous_key: u8) -> Self {
        self.set_byte(new_key, self.copy_index(previous_key));
        self
    }

    /// Change the [`ByteIndexMap`] so that `new_byte` maps to the same thing as
    /// `byte_key`, ignoring case.
    #[inline]
    #[must_use]
    pub const fn add_synonym_ignoring_case(mut self, new_key: u8, previous_key: u8) -> Self {
        self.set_byte_ignoring_case(new_key, self.copy_index(previous_key));
        self
    }

    /// Get the length of `byte_keys`.
    #[inline]
    #[must_use]
    #[allow(clippy::len_without_is_empty)]
    pub const fn len(&self) -> usize {
        self.byte_keys.len()
    }

    /// Convert a base `b` into an index.
    #[inline]
    #[must_use]
    pub const fn to_index(&self, b: u8) -> usize {
        self.index_map[b as usize] as usize
    }

    #[inline]
    #[must_use]
    const fn copy_index(&self, b: u8) -> u8 {
        self.index_map[b as usize]
    }
}

impl<const S: usize> Index<u8> for ByteIndexMap<S> {
    type Output = u8;

    #[inline]
    fn index(&self, index: u8) -> &u8 {
        &self.index_map[index as usize]
    }
}
