//! Collections of k-mers, such as sets and counters.
//!
//! [`KmerSet`] and [`KmerCounter`] offer wrappers around a hash set and hash
//! map respectively, with methods specialized to work on k-mers.
//!
//! For small k-mers, and especially with the two-bit k-mer encoder, the same
//! functionality can be achieved without hashing using [`IndexedKmerSet`] and
//! [`IndexedKmerCounter`]. These use a `Vec<bool>` and `Vec<usize>`
//! respectively to hold the k-mers, where the encoded k-mer is used as an
//! index.

mod counter;
mod indexed_counter;
mod indexed_set;
mod set;

pub use counter::*;
pub use indexed_counter::*;
pub use indexed_set::*;
pub use set::*;
