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
/// TODO: Add more doc links
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

    /// Returns the inner array of a [`ByteMap`].
    #[inline]
    #[must_use]
    pub const fn into_inner(self) -> [u8; 256] {
        self.0
    }

    /// Sets the destination byte for a source byte.
    #[inline]
    pub(crate) const fn set_byte(&mut self, src: u8, dest: u8) {
        self.0[src as usize] = dest;
    }

    /// Changes the [`ByteMap`] so that `new_key` maps to the same thing as
    /// `previous_key`.
    #[inline]
    #[must_use]
    pub const fn add_synonym(mut self, new_key: u8, previous_key: u8) -> Self {
        self.set_byte(new_key, self.copy_mapped_val(previous_key));
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
            out.index_map.set_byte(byte_keys[i].to_ascii_lowercase(), i as u8);
            out.index_map.set_byte(byte_keys[i].to_ascii_uppercase(), i as u8);

            i += 1;
        }
        out
    }

    /// Changes the [`ByteIndexMap`] so that `new_key` maps to the same thing as
    /// `previous_key`.
    #[inline]
    #[must_use]
    pub const fn add_synonym(mut self, new_key: u8, previous_key: u8) -> Self {
        self.index_map = self.index_map.add_synonym(new_key, previous_key);
        self
    }

    /// Sets the index for a byte, ignoring case.
    #[inline]
    pub(crate) const fn set_byte_ignoring_case(&mut self, byte: u8, index: u8) {
        self.index_map.set_byte(byte.to_ascii_lowercase(), index);
        self.index_map.set_byte(byte.to_ascii_uppercase(), index);
    }

    /// Changes the mapping so that `new_key` maps to the same thing as
    /// `previous_key`, ignoring case.
    #[inline]
    #[must_use]
    pub const fn add_synonym_ignoring_case(mut self, new_key: u8, previous_key: u8) -> Self {
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

/// An enum specifying the case of bytes to match in [`ByteMapBuilder`].
#[derive(Clone, Copy, Default)]
pub enum FromCase {
    #[default]
    Exact,
    Any,
}

/// An enum specifying the case of bytes to map to in [`ByteMapBuilder`].
#[derive(Clone, Copy, Default)]
pub enum ToCase {
    /// Uses the case directly as is specified in the destination bytes
    #[default]
    Exact,
    /// Maps to lowercase bytes
    Lower,
    /// Maps to uppercase bytes
    Upper,
    /// Maps to bytes of the same case as the input/source bytes
    Preserve,
}

/// A builder for a [`ByteMap`].
///
/// This struct provides convenient and robust methods for constructing a
/// [`ByteMap`] for a particular set of specifications, all at compile time.
/// Using this has the following steps (all of which are optional, except the
/// first and last):
///
/// 1. Call the [`new`] method to initialize the builder
/// 2. Call [`handle_case`] to enable case-insensitive matching and/or case
///    conversion
/// 3. Define the mapping with [`map`], [`map_many`], and/or [`map_to_self`]
/// 4. Optionally override the default behavior (mapping any other byte to
///    itself) using a function such as [`default_to_byte`]
/// 5. Call the [`build`] method to get a [`ByteMap`]
///
/// [`new`]: ByteMapBuilder::new
/// [`handle_case`]: ByteMapBuilder::handle_case
/// [`map`]: ByteMapBuilder::map
/// [`map_many`]: ByteMapBuilder::map_many
/// [`map_to_self`]: ByteMapBuilder::map_to_self
/// [`default_to_byte`]: ByteMapBuilder::default_to_byte
/// [`build`]: ByteMapBuilder::build
pub struct ByteMapBuilder {
    /// The values that have been set so far
    vals:               [Option<u8>; 256],
    /// The behavior for which case of bytes to match
    from_case:          FromCase,
    /// The behavior for which case of bytes to map to
    to_case:            ToCase,
    /// Whether [`handle_case`] has already been called
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    handle_case_called: bool,
}

impl ByteMapBuilder {
    /// Initializes a new (empty) [`ByteMapBuilder`].
    ///
    /// Unless [`handle_case`] is used to change the behavior, all matching will
    /// be done case-sensitively.
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    #[must_use]
    pub const fn new() -> Self {
        ByteMapBuilder {
            vals:               [None; 256],
            from_case:          FromCase::Exact,
            to_case:            ToCase::Exact,
            handle_case_called: false,
        }
    }

    /// Sets the case-sensitivity and case conversion of the builder.
    ///
    /// By default, [`ByteMapBuilder`] will match bases case-sensitively. To
    /// change this, pass [`FromCase::Any`] for `from_case`. Similarly, the
    /// bytes that are mapped to will by default be the same case as they are
    /// specified. To convert to a different case, use [`ToCase::Lower`] or
    /// [`ToCase::Upper`]. Alternatively, to match the case of the input byte,
    /// use [`ToCase::Preserve`].
    ///
    /// The default behavior is equivalent to calling this with
    /// [`FromCase::Exact`] and [`ToCase::Exact`].
    ///
    /// ## Panics
    ///
    /// This must be called directly after [`new`].
    ///
    /// [`new`]: ByteMapBuilder::new
    #[must_use]
    pub const fn handle_case(mut self, from_case: FromCase, to_case: ToCase) -> ByteMapBuilder {
        assert!(
            !self.handle_case_called,
            "handle_case should only be called once!"
        );

        let mut i = 0;
        while i < 256 {
            assert!(
                self.vals[i].is_none(),
                "handle_case must be called directly after new!"
            );
            i += 1;
        }

        self.handle_case_called = true;
        self.from_case = from_case;
        self.to_case = to_case;
        self
    }

    /// Converts the case of `byte` to the same case as `match_case_of`.
    #[must_use]
    const fn match_case(byte: u8, match_case_of: u8) -> u8 {
        if match_case_of.is_ascii_lowercase() {
            byte.to_ascii_lowercase()
        } else {
            byte.to_ascii_uppercase()
        }
    }

    /// Sets a byte in the [`ByteMapBuilder`].
    #[must_use]
    const fn set_byte_exact_case(mut self, from: u8, dest: u8) -> Self {
        self.vals[from as usize] = Some(dest);
        self
    }

    /// Maps one byte `from` to another byte `dest`.
    ///
    /// This respects the case behavior specified with [`handle_case`].
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    #[must_use]
    const fn map_one(self, from: u8, dest: u8) -> Self {
        match (self.from_case, self.to_case) {
            (FromCase::Exact, ToCase::Exact) => self.set_byte_exact_case(from, dest),
            (FromCase::Exact, ToCase::Lower) => self.set_byte_exact_case(from, dest.to_ascii_lowercase()),
            (FromCase::Exact, ToCase::Upper) => self.set_byte_exact_case(from, dest.to_ascii_uppercase()),
            (FromCase::Exact, ToCase::Preserve) => self.set_byte_exact_case(from, Self::match_case(dest, from)),
            (FromCase::Any, ToCase::Exact) => self
                .set_byte_exact_case(from.to_ascii_lowercase(), dest)
                .set_byte_exact_case(from.to_ascii_uppercase(), dest),
            (FromCase::Any, ToCase::Lower) => self
                .set_byte_exact_case(from.to_ascii_lowercase(), dest.to_ascii_lowercase())
                .set_byte_exact_case(from.to_ascii_uppercase(), dest.to_ascii_lowercase()),
            (FromCase::Any, ToCase::Upper) => self
                .set_byte_exact_case(from.to_ascii_lowercase(), dest.to_ascii_uppercase())
                .set_byte_exact_case(from.to_ascii_uppercase(), dest.to_ascii_uppercase()),
            (FromCase::Any, ToCase::Preserve) => self
                .set_byte_exact_case(from.to_ascii_lowercase(), dest.to_ascii_lowercase())
                .set_byte_exact_case(from.to_ascii_uppercase(), dest.to_ascii_uppercase()),
        }
    }

    /// Maps many bytes `from` to a single byte `dest`.
    ///
    /// This respects the case behavior specified with [`handle_case`].
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    #[must_use]
    pub const fn map_many<const S: usize>(mut self, from: &[u8; S], dest: u8) -> Self {
        let mut i = 0;
        while i < S {
            self = self.map_one(from[i], dest);
            i += 1;
        }
        self
    }

    /// Maps many bytes `from` to corresponding bytes `dest`.
    ///
    /// This respects the case behavior specified with [`handle_case`].
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    #[must_use]
    pub const fn map<const S: usize>(mut self, from: &[u8; S], dest: &[u8; S]) -> Self {
        let mut i = 0;
        while i < S {
            self = self.map_one(from[i], dest[i]);
            i += 1;
        }
        self
    }

    /// Maps many bytes to themselves.
    ///
    /// This respects the case behavior specified with [`handle_case`].
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    #[must_use]
    pub const fn map_to_self<const S: usize>(mut self, bytes: &[u8; S]) -> Self {
        let mut i = 0;
        while i < S {
            self = self.map_one(bytes[i], bytes[i]);
            i += 1;
        }
        self
    }

    /// Maps any unspecified bytes to themselves, but in uppercase.
    #[must_use]
    pub const fn default_to_self_upper(mut self) -> ByteMapBuilder {
        let mut i = 0u8;
        loop {
            if self.vals[i as usize].is_none() {
                self.vals[i as usize] = Some(i.to_ascii_uppercase());
            }
            if i == u8::MAX {
                return self;
            }
            i += 1;
        }
    }

    /// Maps any unspecified bytes to themselves, but in lowercase.
    #[must_use]
    pub const fn default_to_self_lower(mut self) -> ByteMapBuilder {
        let mut i = 0u8;
        loop {
            if self.vals[i as usize].is_none() {
                self.vals[i as usize] = Some(i.to_ascii_lowercase());
            }
            if i == u8::MAX {
                return self;
            }
            i += 1;
        }
    }

    /// Maps any unspecified bytes to a default byte.
    ///
    /// This does NOT use the behavior specified with [`handle_case`]. To
    /// preserve the case of the unspecified byte, use
    /// [`default_to_byte_preserve_case`].
    ///
    /// [`handle_case`]: ByteMapBuilder::handle_case
    /// [`default_to_byte_preserve_case`]: ByteMapBuilder::default_to_byte_preserve_case
    #[must_use]
    pub const fn default_to_byte(mut self, byte: u8) -> ByteMapBuilder {
        let mut i = 0u8;
        loop {
            if self.vals[i as usize].is_none() {
                self.vals[i as usize] = Some(byte);
            }
            if i == u8::MAX {
                return self;
            }
            i += 1;
        }
    }

    /// Maps any unspecified bytes to a default byte, matching its case.
    #[must_use]
    pub const fn default_to_byte_preserve_case(mut self, byte: u8) -> ByteMapBuilder {
        let mut i = 0u8;
        loop {
            if self.vals[i as usize].is_none() {
                self.vals[i as usize] = Some(Self::match_case(byte, i));
            }
            if i == u8::MAX {
                return self;
            }
            i += 1;
        }
    }

    /// Builds a [`ByteMap`] from the builder.
    ///
    /// If no default behavior was specified, then any unset bytes are mapped to
    /// themselves.
    #[must_use]
    pub const fn build(self) -> ByteMap {
        let mut mapping = [0; 256];
        let mut i = 0;
        loop {
            mapping[i as usize] = match self.vals[i as usize] {
                Some(val) => val,
                None => i,
            };
            if i == u8::MAX {
                return ByteMap::new(mapping);
            }
            i += 1;
        }
    }
}

impl Default for ByteMapBuilder {
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}
