//! Structs and traits for indexing inside of pHMMs and query sequences,
//! required for internal or advanced operations on pHMMs.
//!
//! ## Indexing Types
//!
//! Indexing into pHMMs is a difficult problem. The $n$th layer of the pHMM can
//! be ambiguous, depending on whether the BEGIN state is counted. The first
//! match state with emissions is 0 in some contexts and 1 other times.
//! Similarly, with query sequences, index 0 in the dynamic programming table
//! refers to aligning nothing, whereas indexing the sequence at 0 gives the
//! first character. This module intends to provide an abstraction to prevent
//! *Zoe* and advanced users from having bugs related to this, and to promote
//! self-documenting code.
//!
//! The trait [`PhmmIndex`] represents an index in a pHMM (or something related
//! to a pHMM, like a [`SemiLocalModule`]). Anything implementing
//! [`PhmmIndexable`] can be indexed using a [`PhmmIndex`]. The following types
//! can be used as indices:
//!
//! - [`DpIndex`]: An index into a pHMM with respect to a dynamic programming
//!   table, where 0 represents the BEGIN state, 1 represents the first match
//!   state with emissions, etc.
//! - [`SeqIndex`]: An index into a pHMM with respect to the reference
//!   coordinates, where 0 represents the first reference position (or first
//!   match state with emissions), 1 represents the second, etc.
//! - [`Begin`]: An index representing the BEGIN state of the pHMM, equivalent
//!   to `DpIndex(0)`.
//! - [`FirstMatch`]: An index representing the first match state with emissions
//!   or the first reference position, equivalent to `DpIndex(1)` or
//!   `SeqIndex(0)`.
//! - [`LastMatch`]: An index representing the last match state with emissions
//!   or the last reference position, whose value depends on the length of the
//!   pHMM.
//! - [`End`]: An index representing the END state of the pHMM, whose value
//!   depends on the length of the pHMM.
//!
//! The same principle is also used for the query sequence, where [`QueryIndex`]
//! can index into a [`QueryIndexable`]:
//!
//! - [`DpIndex`]: An index into a query with respect to a dynamic programming
//!   table, where 0 represents aligning nothing, 1 represents aligning the
//!   first residue, etc.
//! - [`SeqIndex`]: An index into a query with respect to the sequence
//!   coordinates, where 0 represents the first residue, 1 represents the
//!   second, etc.
//! - [`NoBases`]: An index representing aligning nothing so far, equivalent to
//!   `DpIndex(0)`.
//! - [`FirstBase`]: An index representing the first residue in the sequence,
//!   equivalent to `DpIndex(1)` or `SeqIndex(0)`.
//! - [`LastBase`]: An index representing the last residue in the sequence,
//!   whose value depends on the length of the sequence.
//!
//! [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule

mod accessors;
mod indices;
mod ranges;

pub use accessors::*;
pub use indices::*;
pub use ranges::*;
