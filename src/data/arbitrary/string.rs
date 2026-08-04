//! A specification struct for generating arbitrary [`String`] values.

use crate::data::arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecsView, Case, VecSpecs};
use arbitrary::{Result, Unstructured};

/// Specifications for generating an arbitrary [`String`].
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct StringSpecs {
    /// The character set to which the `u8` bytes must belong.
    pub set: ByteSet,

    /// The case of the ASCII characters.
    pub case: Case,

    // This holds even with invalid UTF-8 based on
    // https://doc.rust-lang.org/nightly/std/str/struct.Utf8Chunk.html#method.invalid,
    // since invalid UTF-8 chunks are at most 3 bytes and the replacement
    // character is 3 bytes
    /// The minimum length of the [`String`], in bytes.
    pub min_len: usize,

    /// The exact length of the [`String`] to generate, in bytes.
    ///
    /// If set, this ignores the `min_len` and `max_len` fields. If a non-ASCII
    /// `set` is chosen, then this requirement might fail to be satisfied during
    /// lossy conversion.
    pub len: Option<usize>,

    /// The maximum length of the [`String`], in bytes.
    ///
    /// If a non-ASCII `set` is chosen, then this requirement might fail to be
    /// satisfied during lossy conversion.
    pub max_len: usize,
}

impl Default for StringSpecs {
    fn default() -> Self {
        Self {
            set:     ByteSet::Any,
            case:    Case::Any,
            min_len: 0,
            len:     None,
            max_len: usize::MAX,
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for StringSpecs {
    type Output = String;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let byte_specs = ByteSpecsView {
            set:  &self.set,
            case: self.case,
        };
        let vec_specs = VecSpecs {
            element_specs: byte_specs,
            min_len:       self.min_len,
            len:           self.len,
            max_len:       self.max_len,
        };
        Ok(String::from_utf8_lossy(&vec_specs.make_arbitrary(u)?).to_string())
    }
}
