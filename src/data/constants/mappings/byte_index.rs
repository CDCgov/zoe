use crate::data::array_types::{self, position};
use std::ops::Index;

/// A mapping from bytes to other bytes. This can be used to represent mappings
/// from residues to u8 numbers, as well as mappings from residues to other
/// residues.
///
/// Also see [`ByteIndexMap`].
#[repr(transparent)]
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub(crate) struct ByteMap([u8; 256]);

impl Default for ByteMap {
    fn default() -> Self {
        Self([0; 256])
    }
}

/// Create a new [`ByteMap`] through syntax such as the following:
/// ```ignore
/// pub(crate) const DNA_MAP: ByteMap = byte_map! {
///     @ignoring_case
///     b"ACG" => (0..=2),
///     b"TU" => 3,
///     b"N" => 4,
///     b"-." => 5
///     b"RYSWKMBDHV" => 6,
///     _ => 7
/// };
/// ```
///
/// The first line is either `@ignoring_case` or `@matching_case` to specify
/// case sensitivity. All other lines are match lines. One of the lines (with a
/// left side of `_`) specifies the default value. All other lines have a byte
/// string on the left, and on the right side, there is either a single byte
/// (all of the bytes on the left map to this byte) or a range of bytes (all of
/// the bytes on the left map to the consecutive bytes in the range).
macro_rules! byte_map {
    // Entry point for macro: initialization and parse each match line
    (@$kind:tt $($src:tt => $dest:tt),* $(,)?) => {{
        let mut set = [false; 256];
        #[allow(unused_assignments)]
        let mut default_found = false;
        let mut array = ByteMap::new([0; 256]);
        $(
            byte_map!(@$kind array, set, default_found, $src, $dest);
        )*
        assert!(default_found, "No default value was set!");
        array
    }};

    // Parse the default match line
    (@$kind:tt $array:ident, $set:ident, $default_found:ident, _, $dest:expr) => {{
        assert!(!$default_found, "Attempted to set a default value more than once!");
        $default_found = true;
        let mut i = 0;
        while (i as usize) < 255 {
            if !$set[i as usize] {
                $array.set_byte(i, $dest);
            }
            i += 1;
        }
        if !$set[i as usize] {
            $array.set_byte(i, $dest);
        }
    }};

    // Parse an end-exclusive range
    (@$kind:tt $array:ident, $set:ident, $default_found:ident, $src:expr, ($start:literal .. $end:literal)) => {{
        let mut i = 0;
        while i < $src.len() {
            byte_map!(#set $kind, $array, $set, $src[i], $start + i as u8);
            i += 1;
        }
        assert!(i == $end, "Incorrect length of range passed!")
    }};

    // Parse an end-inclusive range
    (@$kind:tt $array:ident, $set:ident, $default_found:ident, $src:expr, ($start:literal ..= $end:literal)) => {{
        let mut i = 0;
        while (i as usize) < $src.len() {
            byte_map!(#set $kind, $array, $set, $src[i as usize], $start + i);
            i += 1;
        }
        assert!($start + i == $end + 1, "Incorrect length of range passed!")
    }};

    // Parse a regular match line
    (@$kind:tt $array:ident, $set:ident, $default_found:ident, $src:expr, $dest:expr) => {
        let mut i = 0;
        while i < $src.len() {
            byte_map!(#set $kind, $array, $set, $src[i], $dest);
            i += 1;
        }
    };

    // Set a byte, case-sensitive
    (#set matching_case, $array:ident, $set:ident, $src:expr, $dest:expr) => {
        assert!(!$set[$src as usize], "Attempted to set byte multiple times!");
        $set[$src as usize] = true;
        $array.set_byte($src, $dest)
    };

    // Set a byte, case-insensitive
    (#set ignoring_case, $array:ident, $set:ident, $src:expr, $dest:expr) => {
        assert!(!$set[$src as usize], "Attempted to set byte multiple times!");
        $set[$src.to_ascii_lowercase() as usize] = true;
        $set[$src.to_ascii_uppercase() as usize] = true;
        $array.set_byte_ignoring_case($src, $dest)
    };
}

impl ByteMap {
    /// Create a new custom [`ByteMap`] from an existing mapping.
    pub const fn new(mapping: [u8; 256]) -> Self {
        Self(mapping)
    }

    /// Set the index for a byte.
    #[inline]
    pub(crate) const fn set_byte(&mut self, byte: u8, index: u8) {
        self.0[byte as usize] = index;
    }

    /// Set the index for a byte, ignoring case.
    #[inline]
    pub(crate) const fn set_byte_ignoring_case(&mut self, byte: u8, index: u8) {
        self.set_byte(byte.to_ascii_lowercase(), index);
        self.set_byte(byte.to_ascii_uppercase(), index);
    }

    /// Change the [`ByteMap`] so that `new_key` maps to the same thing as
    /// `previous_key`.
    #[inline]
    #[must_use]
    pub const fn add_synonym(mut self, new_key: u8, previous_key: u8) -> Self {
        self.set_byte(new_key, self.copy_mapped_val(previous_key));
        self
    }

    /// Change the [`ByteMap`] so that `new_key` maps to the same thing as
    /// `previous_key`, ignoring case.
    #[inline]
    #[must_use]
    pub const fn add_synonym_ignoring_case(mut self, new_key: u8, previous_key: u8) -> Self {
        self.set_byte_ignoring_case(new_key, self.copy_mapped_val(previous_key));
        self
    }

    /// Converts a base `b` into an index.
    #[inline]
    #[must_use]
    pub const fn to_index(&self, b: u8) -> usize {
        self.0[b as usize] as usize
    }

    /// Copies the value to which `b` maps.
    #[inline]
    #[must_use]
    const fn copy_mapped_val(&self, b: u8) -> u8 {
        self.0[b as usize]
    }
}

/// Represents a mapping between bytes and indices. For example, this could be a
/// map from DNA bases to profile indices, such as [`DNA_PROFILE_MAP`].
///
/// ## Type Parameters
///
/// * `S` - The number of unique indices in the output of the mapping (such as 5
///   for DNA including *N*)
///
/// [`DNA_PROFILE_MAP`]: super::DNA_PROFILE_MAP
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct ByteIndexMap<const S: usize> {
    pub(crate) index_map: ByteMap,
    pub(crate) byte_keys: [u8; S],
}

impl<const S: usize> ByteIndexMap<S> {
    /// Create a new [`ByteIndexMap`] struct to represent a mapping between
    /// bytes and indices. For example, this could be a map from DNA bases to
    /// profile indices. Any byte that is not specified in `byte_keys` is mapped
    /// to the same thing as `catch_all`.
    ///
    /// ## Panics
    ///
    /// No duplicates can be present in `byte_keys`. `catch_all` must be present
    /// in `byte_keys`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new(byte_keys: [u8; S], catch_all: u8) -> Self {
        assert!(array_types::is_unique(&byte_keys));

        let catch_all_index =
            position(&byte_keys, catch_all).expect("The catch_all must be present in the byte_keys.") as u8;

        let mut out = Self {
            index_map: ByteMap([catch_all_index; 256]),
            byte_keys,
        };

        let mut i = 0;
        while i < byte_keys.len() {
            // Truncation will not occur because i cannot exceed
            // byte_keys.len(), and index must contain unique u8 values
            out.index_map.set_byte(byte_keys[i], i as u8);

            i += 1;
        }
        out
    }

    /// Create a new [`ByteIndexMap`] struct to represent a mapping between
    /// bytes and indices. For example, this could be a map from DNA bases to
    /// profile indices. Both `byte_keys` and `catch_all` ignore case.
    ///
    /// ## Panics
    ///
    /// No duplicates can be present in `byte_keys`. `catch_all` must be present
    /// in `byte_keys`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new_ignoring_case(mut byte_keys: [u8; S], catch_all: u8) -> Self {
        byte_keys = array_types::make_uppercase(&byte_keys);
        assert!(array_types::is_unique(&byte_keys));

        let catch_all_index = position(&byte_keys, catch_all.to_ascii_uppercase())
            .expect("The catch_all must be present in the byte_keys.") as u8;

        let mut out = Self {
            index_map: ByteMap([catch_all_index; 256]),
            byte_keys,
        };

        let mut i = 0;
        while i < byte_keys.len() {
            // Truncation will not occur because i cannot exceed
            // byte_keys.len(), and byte_keys must contain unique u8 values
            out.index_map.set_byte_ignoring_case(byte_keys[i], i as u8);

            i += 1;
        }
        out
    }

    /// Change the [`ByteIndexMap`] so that `new_key` maps to the same thing as
    /// `previous_key`.
    #[inline]
    #[must_use]
    pub const fn add_synonym(mut self, new_key: u8, previous_key: u8) -> Self {
        self.index_map = self.index_map.add_synonym(new_key, previous_key);
        self
    }

    /// Change the [`ByteIndexMap`] so that `new_key` maps to the same thing as
    /// `previous_key`, ignoring case.
    #[inline]
    #[must_use]
    pub const fn add_synonym_ignoring_case(mut self, new_key: u8, previous_key: u8) -> Self {
        self.index_map = self.index_map.add_synonym_ignoring_case(new_key, previous_key);
        self
    }

    /// Gets the number of unique indices in the output of this mapping.
    #[inline]
    #[must_use]
    #[allow(clippy::len_without_is_empty)]
    pub const fn len(&self) -> usize {
        S
    }

    // TODO: Consider using const fn in trait when this becomes available
    /// Converts a byte `b` into an index.
    #[inline]
    #[must_use]
    pub const fn to_index(&self, b: u8) -> usize {
        self.index_map.to_index(b)
    }

    /// Converts an index back into the corresponding byte.
    #[inline]
    #[must_use]
    pub const fn to_byte(&self, index: usize) -> u8 {
        self.byte_keys[index]
    }

    /// Returns whether a byte `b` is in `byte_keys`.
    #[inline]
    #[must_use]
    pub const fn in_byte_keys(&self, b: u8) -> bool {
        b == self.to_byte(self.to_index(b))
    }
}

impl Index<u8> for ByteMap {
    type Output = u8;

    #[inline]
    fn index(&self, index: u8) -> &u8 {
        &self.0[index as usize]
    }
}

impl<const S: usize> Index<u8> for ByteIndexMap<S> {
    type Output = u8;

    #[inline]
    fn index(&self, index: u8) -> &u8 {
        &self.index_map[index]
    }
}

pub(crate) use byte_map;
