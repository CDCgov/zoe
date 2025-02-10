//! ## Sequence search and/or replacement.
//!
//! This module provides functions for performing searches and/or replacements
//! within byte slices or vectors. In addition to many customized search
//! functions, two general-purpose search algorithms are [`substring_match`] and
//! [`substring_match_simd`]. The latter is more efficient for larger needles.
//!
//! Similarly, for fuzzy matching (permitting some number of mismatches between
//! the desired and located needle), [`fuzzy_substring_match`] and
//! [`fuzzy_substring_match_simd`] are provided.
//!
//! The [`ByteSubstring`] trait provides methods which wrap these underlying
//! algorithms, such as [`contains_substring`], [`find_substring`], and
//! [`find_fuzzy_substring`]. These methods return a range instead of just the
//! starting index.
//!
//! ## Restricting the search range
//!
//! Sometimes, you may want to search a particular region of a string. There are
//! two main scenarios:
//! 1. The string search should return an index/range with respect to the
//!    subsequence. For example, searching `ACGT` for `G` in the range `1..3`
//!    would return a starting index of 1.
//! 2. The string search should return an index/range with respect to the full,
//!    original sequence. Using the same example, the starting index would be 2.
//!
//! The first scenario is easily resolved by creating a byte slice or a
//! [view](crate::data#views), and then searching within that. For example:
//! ```
//! # use zoe::prelude::*;
//! let haystack = b"CACATAACGTACAGTTTTACACAGTTAGGT";
//! let needle = b"CACA";
//! let position = haystack[10..25].find_substring(needle);
//! assert_eq!(position, Some(9..13))
//! ```
//!
//! The second scenario could also be solved the same way, but it would require
//! shifting the output by the starting index. An alternative is the
//! [`RangeSearch`] struct, which implements [`ByteSubstring`] and automatically
//! adjusts the range afterwards. An example usage is below. Notice that the
//! potential match at the beginning of the haystack is not returned, since it
//! is not in the provided range.
//! ```
//! # use zoe::{search::ToRangeSearch, prelude::*};
//! let haystack = b"CACATAACGTACAGTTTTACACAGTTAGGT";
//! let needle = b"CACA";
//! let position = haystack.search_in(10..25).find_substring(needle);
//! assert_eq!(position, Some(19..23))
//! ```
//!
//! [`RangeSearch`] also provides additional convenience methods. One use-case
//! is when a needle must be found in at most the first `n` elements of the
//! haystack, or at most the last `n` elements of the haystack. Ordinarily, one
//! would have to ensure the range is in bounds before slicing or using
//! `search_in`. To simplify this, *Zoe* provides `search_in_first` and
//! `search_in_last`, as shown below. Once again, the second occurrence of the
//! needle is returned.
//! ```
//! # use zoe::{search::ToRangeSearch, prelude::*};
//! let haystack = b"CACATAACGTACAGTTTTACACAGTTAGGT";
//! let needle = b"CACA";
//! let position = haystack.search_in_last(15).find_substring(needle);
//! assert_eq!(position, Some(19..23))
//! ```
//!
//! [`substring_match`]: crate::search::substring_match
//! [`substring_match_simd`]: crate::search::substring_match_simd
//! [`fuzzy_substring_match`]: crate::search::fuzzy_substring_match
//! [`fuzzy_substring_match_simd`]: crate::search::fuzzy_substring_match_simd
//! [`ByteSubstring`]: crate::search::ByteSubstring
//! [`contains_substring`]: crate::search::ByteSubstring::contains_substring
//! [`find_substring`]: crate::search::ByteSubstring::find_substring
//! [`find_fuzzy_substring`]: crate::search::ByteSubstring::find_fuzzy_substring
//! [`RangeSearch`]: crate::search::RangeSearch

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
