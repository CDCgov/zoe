/// Search and/or replace bytes.
mod bytes;
/// Search for fuzzy substrings.
mod inexact;
/// Search byte substrings where needle is a repeating character.
mod k_repeating;
/// Types for representing restricted range searches
mod range_search;
/// Search byte substrings.
mod substring;

pub use bytes::*;
pub use inexact::*;
pub use k_repeating::*;
pub use range_search::{RangeSearch, ToRangeSearch};
pub use substring::*;
