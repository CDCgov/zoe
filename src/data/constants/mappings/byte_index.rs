use crate::data::array_types::{self, position};
use std::ops::Index;

/// A mapping from bytes to other bytes.
///
/// This struct has many potential use-cases:
///
/// 1. Representing a mapping from residues to other residues, such as for
///    switching alphabets.
/// 2. Mapping residues while also filtering them, by mapping invalid residues
///    to 0.
/// 3. Mapping residues to integers, as may be necessary for indexing or
///    encoding. If each integer has a single canonical byte that it represents,
///    consider using [`ByteIndexMap`] which also allows for converting back to
///    regular sequences again.
///
/// To build a custom map, start with a complete mapping such as
/// [`ByteMap::identity`] or [`ByteMap::all`], then apply cumulative overrides
/// with methods such as [`ByteMap::preserve`], [`ByteMap::many_to_one`],
/// [`ByteMap::map`], or [`ByteMap::indexing`].
#[repr(transparent)]
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct ByteMap([u8; 256]);

impl Default for ByteMap {
    #[inline]
    fn default() -> Self {
        Self([0; 256])
    }
}

impl ByteMap {
    /// Creates a new custom [`ByteMap`] from an existing mapping.
    #[inline]
    #[must_use]
    pub const fn new(mapping: [u8; 256]) -> Self {
        Self(mapping)
    }

    /// Creates a [`ByteMap`] that maps every byte to the same destination byte.
    #[inline]
    #[must_use]
    pub const fn all(dest: u8) -> Self {
        Self([dest; 256])
    }

    /// Creates a [`ByteMap`] that maps each byte to itself.
    #[must_use]
    pub const fn identity() -> Self {
        let mut mapping = [0; 256];
        let mut i = 0u8;
        loop {
            mapping[i as usize] = i;
            if i == u8::MAX {
                return Self(mapping);
            }
            i += 1;
        }
    }

    /// Returns the inner array of a [`ByteMap`].
    #[inline]
    #[must_use]
    pub const fn into_inner(self) -> [u8; 256] {
        self.0
    }

    /// Preserves the identity of the given bytes in the map.
    ///
    /// Any later calls can still override these mappings.
    #[must_use]
    pub const fn preserve<const S: usize>(self, bytes: &[u8; S]) -> Self {
        self.map(bytes, bytes)
    }

    /// Preserves the given bytes in both ASCII cases.
    ///
    /// ASCII alphabetic bytes preserve the case of the matched input byte.
    /// Non-alphabetic bytes are preserved unchanged.
    #[must_use]
    pub const fn preserve_both_cases<const S: usize>(mut self, bytes: &[u8; S]) -> Self {
        Self::assert_unique_ignore_case(bytes);

        let mut i = 0;
        while i < S {
            self.set_byte(
                bytes[i].to_ascii_lowercase(),
                Self::match_case(bytes[i], bytes[i].to_ascii_lowercase()),
            );
            self.set_byte(
                bytes[i].to_ascii_uppercase(),
                Self::match_case(bytes[i], bytes[i].to_ascii_uppercase()),
            );
            i += 1;
        }
        self
    }

    /// Maps many source bytes to one destination byte.
    #[must_use]
    pub const fn many_to_one<const S: usize>(mut self, from: &[u8; S], to: u8) -> Self {
        Self::assert_unique(from);

        let mut i = 0;
        while i < S {
            self.set_byte(from[i], to);
            i += 1;
        }
        self
    }

    /// Maps many source bytes to one destination byte while matching their
    /// uppercase and lowercase forms.
    ///
    /// The destination byte is used literally. Use [`Self::preserve_both_cases`]
    /// when the output should follow the matched input case.
    #[must_use]
    pub const fn many_to_one_ignore_case<const S: usize>(mut self, from: &[u8; S], to: u8) -> Self {
        Self::assert_unique_ignore_case(from);

        let mut i = 0;
        while i < S {
            self.set_byte(from[i].to_ascii_lowercase(), to);
            self.set_byte(from[i].to_ascii_uppercase(), to);
            i += 1;
        }
        self
    }

    /// Maps source bytes to destination bytes pairwise.
    #[must_use]
    pub const fn map<const S: usize>(mut self, from: &[u8; S], to: &[u8; S]) -> Self {
        Self::assert_unique(from);

        let mut i = 0;
        while i < S {
            self.set_byte(from[i], to[i]);
            i += 1;
        }
        self
    }

    /// Maps source bytes to destination bytes pairwise while matching their
    /// uppercase and lowercase forms.
    ///
    /// Destination bytes are used literally. Use [`Self::preserve_both_cases`] when
    /// the output should follow the matched input case.
    #[must_use]
    pub const fn map_ignore_case<const S: usize>(mut self, from: &[u8; S], to: &[u8; S]) -> Self {
        Self::assert_unique_ignore_case(from);

        let mut i = 0;
        while i < S {
            self.set_byte(from[i].to_ascii_lowercase(), to[i]);
            self.set_byte(from[i].to_ascii_uppercase(), to[i]);
            i += 1;
        }
        self
    }

    /// Maps the provided bytes to ascending indexes.
    ///
    /// This maps `bytes[0]` to `0`, `bytes[1]` to `1`, and so on.
    ///
    /// ## Panics
    ///
    /// `bytes` must contain at most 256 unique values.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn indexing<const S: usize>(mut self, bytes: &[u8; S]) -> Self {
        assert!(S <= 256, "Cannot index more than 256 distinct bytes!");
        Self::assert_unique(bytes);

        let mut i = 0;
        while i < S {
            self.set_byte(bytes[i], i as u8);
            i += 1;
        }
        self
    }

    /// Maps the provided bytes to ascending indexes while matching their
    /// uppercase and lowercase forms.
    ///
    /// ## Panics
    ///
    /// `bytes` must contain at most 256 values that are unique after
    /// lowercasing them.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn indexing_ignore_case<const S: usize>(mut self, bytes: &[u8; S]) -> Self {
        assert!(S <= 256, "Cannot index more than 256 distinct bytes!");
        Self::assert_unique_ignore_case(bytes);

        let mut i = 0;
        while i < S {
            self.set_byte(bytes[i].to_ascii_lowercase(), i as u8);
            self.set_byte(bytes[i].to_ascii_uppercase(), i as u8);
            i += 1;
        }
        self
    }

    /// Sets the destination byte for a source byte.
    #[inline]
    pub(crate) const fn set_byte(&mut self, src: u8, dest: u8) {
        self.0[src as usize] = dest;
    }

    /// Converts a base `b` into an index.
    #[inline]
    #[must_use]
    pub const fn to_index(&self, b: u8) -> usize {
        self.0[b as usize] as usize
    }

    /// Converts the case of `byte` to match `match_case_of`.
    #[inline]
    #[must_use]
    const fn match_case(byte: u8, match_case_of: u8) -> u8 {
        if match_case_of.is_ascii_lowercase() {
            byte.to_ascii_lowercase()
        } else {
            byte.to_ascii_uppercase()
        }
    }

    /// Ensures a byte slice does not contain duplicates.
    #[inline]
    const fn assert_unique<const S: usize>(bytes: &[u8; S]) {
        assert!(
            array_types::is_unique(bytes),
            "Attempted to map a byte multiple times in the same call!"
        );
    }

    /// Ensures a byte slice does not contain duplicates after lowercasing them.
    #[inline]
    const fn assert_unique_ignore_case<const S: usize>(bytes: &[u8; S]) {
        let mut i = 0;
        while i < S {
            let mut j = i + 1;
            while j < S {
                assert!(
                    !bytes[i].eq_ignore_ascii_case(&bytes[j]),
                    "Attempted to map a byte multiple times in the same call!"
                );
                j += 1;
            }
            i += 1;
        }
    }

    /// Copies the value to which `b` maps.
    #[inline]
    #[must_use]
    const fn copy_mapped_val(&self, b: u8) -> u8 {
        self.0[b as usize]
    }
}

/// Represents a mapping between bytes and indices.
///
/// For example, this could be a map from DNA bases to profile indices, such as
/// [`DNA_PROFILE_MAP`].
///
/// ## Type Parameters
///
/// `S`: The number of unique indices in the output of the mapping (such as 5
/// for DNA including *N*)
///
/// [`DNA_PROFILE_MAP`]: super::DNA_PROFILE_MAP
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct ByteIndexMap<const S: usize> {
    index_map: ByteMap,
    byte_keys: [u8; S],
}

impl<const S: usize> ByteIndexMap<S> {
    /// Creates a new [`ByteIndexMap`] struct to represent a mapping between
    /// bytes and indices.
    ///
    /// For example, this could be a map from DNA bases to profile indices. Any
    /// byte that is not specified in `byte_keys` is mapped to the same thing as
    /// `catch_all`.
    ///
    /// ## Panics
    ///
    /// No duplicates can be present in `byte_keys`. `catch_all` must be present
    /// in `byte_keys`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new(byte_keys: [u8; S], catch_all: u8) -> Self {
        let catch_all_index =
            position(&byte_keys, catch_all).expect("The catch_all must be present in the byte_keys.") as u8;

        Self {
            // `indexing` validates that `byte_keys` are unique.
            index_map: ByteMap::all(catch_all_index).indexing(&byte_keys),
            byte_keys,
        }
    }

    /// Creates a new [`ByteIndexMap`] struct to represent a mapping between
    /// bytes and indices, ignoring case.
    ///
    /// For example, this could be a map from DNA bases to profile indices. Both
    /// `byte_keys` and `catch_all` ignore case.
    ///
    /// ## Panics
    ///
    /// No duplicates can be present in `byte_keys`. `catch_all` must be present
    /// in `byte_keys`.
    #[allow(clippy::cast_possible_truncation)]
    #[must_use]
    pub const fn new_ignoring_case(mut byte_keys: [u8; S], catch_all: u8) -> Self {
        byte_keys = array_types::make_uppercase(&byte_keys);

        let catch_all_index = position(&byte_keys, catch_all.to_ascii_uppercase())
            .expect("The catch_all must be present in the byte_keys.") as u8;

        Self {
            // `indexing_ignore_case` validates that `byte_keys` are unique after
            // lowercasing them.
            index_map: ByteMap::all(catch_all_index).indexing_ignore_case(&byte_keys),
            byte_keys,
        }
    }

    /// Changes the [`ByteIndexMap`] so that `new_key` maps to the same thing as
    /// `previous_key`.
    #[inline]
    #[must_use]
    pub const fn add_synonym(mut self, new_key: u8, previous_key: u8) -> Self {
        self.index_map.set_byte(new_key, self.index_map.copy_mapped_val(previous_key));
        self
    }

    /// Sets the index for a byte, ignoring case.
    #[inline]
    pub(crate) const fn set_byte_ignoring_case(&mut self, byte: u8, index: u8) {
        self.index_map.set_byte(byte.to_ascii_lowercase(), index);
        self.index_map.set_byte(byte.to_ascii_uppercase(), index);
    }

    /// Adds source bytes as synonyms for the canonical byte, ignoring case.
    #[inline]
    #[must_use]
    pub const fn add_synonym_ignore_case(mut self, new_key: u8, previous_key: u8) -> Self {
        self.set_byte_ignoring_case(new_key, self.index_map.copy_mapped_val(previous_key));
        self
    }

    /// Gets the number of unique indices in the output of this mapping.
    #[inline]
    #[must_use]
    #[allow(clippy::len_without_is_empty)]
    pub const fn len(&self) -> usize {
        S
    }

    /// Gets the byte keys for this mapping.
    #[inline]
    #[must_use]
    pub const fn byte_keys(&self) -> &[u8; S] {
        &self.byte_keys
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
