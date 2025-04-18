use crate::{data::ByteIndexMap, math::Uint};
use std::ops::{Add, AddAssign};

/// Count statistics for use with [`ByteIndexMap`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct ByteIndexCounts<T: Uint, const S: usize> {
    mapping: &'static ByteIndexMap<S>,
    inner:   [T; S],
}

impl<T: Uint, const S: usize> ByteIndexCounts<T, S> {
    /// Creates a new [`ByteIndexCounts`] object with counts initialized to 0.
    #[inline]
    #[must_use]
    pub fn new(mapping: &'static ByteIndexMap<S>) -> Self {
        Self {
            mapping,
            inner: [T::ZERO; S],
        }
    }

    /// Increments the counts using the bytes in `seq`.
    #[inline]
    pub fn tally_from_seq<Q: AsRef<[u8]>>(&mut self, seq: Q) {
        for byte in seq.as_ref() {
            *self += *byte;
        }
    }

    /// Retrives the counts as an array of size S.
    #[inline]
    #[must_use]
    pub fn into_inner(self) -> [T; S] {
        self.inner
    }

    /// Retrieves the count of a byte.
    ///
    /// ## Panics
    ///
    /// The byte must be present in the `byte_keys` of the associated mapping.
    pub fn get_count(&self, byte: u8) -> T {
        let index = self.mapping.to_index(byte);
        if byte == self.mapping.to_byte(index) {
            return self.inner[index];
        }
        panic!("The specified byte is valid for the mapping")
    }
}

impl<T, const S: usize> Add<u8> for ByteIndexCounts<T, S>
where
    T: Uint,
{
    type Output = Self;

    #[inline]
    fn add(mut self, other: u8) -> Self {
        self.inner[self.mapping.to_index(other)] += T::ONE;
        self
    }
}

impl<T, const S: usize> AddAssign<u8> for ByteIndexCounts<T, S>
where
    T: Uint,
{
    #[inline]
    fn add_assign(&mut self, other: u8) {
        self.inner[self.mapping.to_index(other)] += T::ONE;
    }
}
