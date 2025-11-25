//! A specification struct for generating arbitrary [`String`] values.

use crate::data::arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecs, Case, VecSpecs};
use arbitrary::{Result, Unstructured};

/// Specifications for generating an arbitrary [`String`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct StringSpecs {
    /// The character set to which the `u8` bytes must belong.
    pub set: ByteSet,

    /// The case of the ASCII characters.
    pub case: Case,
}

impl From<ByteSpecs> for StringSpecs {
    #[inline]
    fn from(value: ByteSpecs) -> Self {
        Self {
            set:  value.set,
            case: value.case,
        }
    }
}

impl From<StringSpecs> for ByteSpecs {
    #[inline]
    fn from(value: StringSpecs) -> Self {
        Self {
            set:  value.set,
            case: value.case,
        }
    }
}

impl<'a> ArbitrarySpecs<'a> for StringSpecs {
    type Output = String;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let byte_specs = ByteSpecs::from(*self);
        let vec_specs = VecSpecs {
            element_specs: byte_specs,
            min_len:       0,
            len:           None,
            max_len:       usize::MAX,
        };
        Ok(String::from_utf8_lossy(&vec_specs.make_arbitrary(u)?).to_string())
    }
}
