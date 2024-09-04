/// Search and/or replace bytes.
mod bytes;
/// Search for fuzzy substrings
mod inexact;
/// Search byte substrings
mod substring;

pub use bytes::*;
pub use inexact::*;
pub use substring::*;
