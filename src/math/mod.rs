mod cast;
mod float;
mod integer;

pub use crate::simd::SimdAnyInt;
pub use float::*;
pub use integer::*;

#[cfg(any(feature = "fuzzing", test))]
mod float_compare;
#[cfg(any(feature = "fuzzing", test))]
pub use float_compare::*;

#[cfg(feature = "rand")]
mod rand;
#[cfg(feature = "rand")]
pub use rand::*;
