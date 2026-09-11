//! Structs and traits for indexing into pHMMs, useful for diagnostics or
//! advanced pHMM usage.
//!
//! ## Indexing Types
//!
//! There are nuances when indexing for dynamic programming alignment problems
//! because the sequence coordinates are not the same as the indices into the DP
//! table. A dynamic programming table for alignment uses index 0 to hold "no
//! bases aligned", and index 1 to hold sequence coordinate 0. Indexing from the
//! end can also be similarly challenging and context-dependent.
//!
//! The same problem exists for pHMMs. Indexing the layers of a pHMM can be
//! ambiguous, depending on whether the BEGIN state is counted. The first match
//! state with emissions is index 0 in some contexts and index 1 other times.
//!
//! This module intends to provide an abstraction to prevent *Zoe* and advanced
//! users from having bugs related to this, and to promote self-documenting
//! code.
//!
//! The trait [`AlnIndex`] represents an index in a pHMM, something related to a
//! pHMM like a [`SemiLocalModule`], or a sequence. Anything implementing
//! [`AlnIndexable`] can be indexed using a [`AlnIndex`]. The following types
//! can be used as indices:
//!
//! - [`DpIndex`]: An index with respect to a dynamic programming table. For
//!   sequences, 0 represents no residues aligned and 1 represents the first
//!   residue. For pHMMs, 0 represents the BEGIN state and 1 represents the
//!   first match state with emissions.
//! - [`SeqIndex`]: An index with respect to the sequence coordinates. For
//!   sequences, 0 represents the first residue. For pHMMs, 0 represents the
//!   first match state with emissions.
//! - [`Begin`]: An index representing no residues aligned or the BEGIN state of
//!   a pHMM, equivalent to `DpIndex(0)`.
//! - [`FirstResidue`]: An index representing the first residue or the first
//!   match state with emissions, equivalent to `DpIndex(1)` or `SeqIndex(0)`.
//! - [`LastResidue`]: An index representing the last residue or the last match
//!   state with emissions, whose value depends on the length of the
//!   sequence/pHMM.
//! - [`End`]: An index after [`LastResidue`]. This represents the END state of
//!   a pHMM, and also occurs as an exclusive end bound on ranges.
//!
//! [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule

mod accessors;
mod aln_index;
mod aln_indexable;
mod ranges;

pub use accessors::*;
pub use aln_index::*;
pub use aln_indexable::*;
pub use ranges::*;
