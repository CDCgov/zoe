//! ## K-mer processing, searching, and counting.
//!
//! This module provides structs and methods for handling common k-mer
//! operations efficiently using integer encodings. A k-mer is a short
//! subsequence of nucleotides, which are represented in *Zoe* by [`Kmer`].
//!
//! *Zoe* has two structs to store k-mers:
//! * A [`KmerSet`], which is a [`HashSet`] for k-mers
//! * A [`KmerCounter`] to store k-mers and their counts (using a [`HashMap`])
//!
//! These structs can be populated from individual k-mers, sequences, or k-mers
//! with up to `N` mismatches. They can be then be queried, iterated over, or
//! used to search within a sequence. See [`KmerSet`] and [`KmerCounter`] for
//! more details.
//!
//! ## Generic Parameters
//!
//! Many structs related to k-mers, including [`KmerSet`] and [`KmerCounter`],
//! are generic over a maximum possible k-mer length `MAX_LEN`. This is used to
//! determine (at compile-time) the appropriate integer type to store the
//! encoded k-mer. The actual k-mer length can be set to a different value at
//! runtime as long as it is less than `MAX_LEN`. For guidance on picking the
//! appropriate `MAX_LEN`, see [`SupportedKmerLen`].
//!
//! [`KmerSet`] and [`KmerCounter`] are also generic over the type of
//! [`KmerEncoder`]. *Zoe* currently provides a single encoder,
//! [`ThreeBitKmerEncoder`], which uses three bits to store each base. It allows
//! for `A`, `C`, `G`, `T`, and `N` to all be represented. It does not preserve
//! case or the distinction between `T` and `U`. `N` is used as a catch-all for
//! bases that are not `ACGTUNacgtun`.
//!
//! For convenience, the type aliases [`ThreeBitKmerSet`] and
//! [`ThreeBitKmerCounter`] are provided.
//!
//! <div class="warning important">
//!
//! **Important**
//!
//! When performing more specialized k-mer operations, you may need to directly
//! encode and decode k-mers, rather than relying on the methods in [`KmerSet`]
//! or [`KmerCounter`]. It is important to use the same [`KmerEncoder`] to both
//! encode and decode the k-mers.
//!
//! </div>
//!
//! ## Examples
//!
//! Count the 3-mers in a sequence:
//! ```
//! # use zoe::{kmer::encoders::three_bit::ThreeBitKmerCounter, prelude::*};
//! let sequence = b"GGCCACCAAGGCCA";
//! let mut kmer_counter = ThreeBitKmerCounter::<3>::new(3).unwrap();
//! kmer_counter.tally_from_sequence(sequence);
//! for (kmer, count) in kmer_counter {
//!     println!("{kmer}\t{count}");
//! }
//! ```
//!
//! Search for the 17-mers of a primer within a sequence, with up to one
//! mismatch:
//! ```
//! # use zoe::{kmer::encoders::three_bit::ThreeBitKmerSet, prelude::*};
//! let primer = b"TGATAGTTTTAGAGTTAGGTAG";
//! let sequence = b"TGCCCGTAACGTACAGTTTTACAGTTAGGTACCC";
//! let mut kmer_set = ThreeBitKmerSet::<17>::new(17).unwrap();
//! kmer_set.insert_from_sequence_with_variants::<1>(primer);
//! let kmer_pos = kmer_set.find_in_seq(sequence);
//! assert_eq!(kmer_pos, Some(14..31));
//! ```
//!
//! This can be equivalent performed using:
//! ```
//! # use zoe::{kmer::{encoders::three_bit::ThreeBitKmerSet, FindKmers}, prelude::*};
//! let primer: Nucleotides = b"TGATAGTTTTAGAGTTAGGTAG".into();
//! let sequence: Nucleotides = b"TGCCCGTAACGTACAGTTTTACAGTTAGGTACCC".into();
//! let mut kmer_set = ThreeBitKmerSet::<17>::new(17).unwrap();
//! kmer_set.insert_from_sequence_with_variants::<1>(primer);
//! let kmer_pos = sequence.find_kmers(&kmer_set);
//! assert_eq!(kmer_pos, Some(14..31));
//! ```
//!
//! [`HashSet`]: std::collections::HashSet
//! [`HashMap`]: std::collections::HashMap
//! [`Display`]: std::fmt::Display
//! [`Kmer`]: Kmer
//! [`KmerEncoder`]: KmerEncoder
//! [`KmerSet`]: KmerSet
//! [`KmerCounter`]: KmerCounter
//! [`ThreeBitKmerEncoder`]: encoders::three_bit::ThreeBitKmerEncoder
//! [`ThreeBitEncodedKmer`]: encoders::three_bit::ThreeBitEncodedKmer
//! [`ThreeBitKmerSet`]: encoders::three_bit::ThreeBitKmerSet
//! [`ThreeBitKmerCounter`]: encoders::three_bit::ThreeBitKmerCounter
//! [`SupportedKmerLen`]: SupportedKmerLen
//! [`insert_from_sequence`]: KmerSet::insert_from_sequence
//! [`find_in_seq`]: FindKmersInSeq::find_in_seq
//! [`find_in_seq_rev`]: FindKmersInSeq::find_in_seq_rev
//! [`FindKmers`]: FindKmers
//! [`find_kmers`]: FindKmers::find_kmers
//! [`find_kmers_rev`]: FindKmers::find_kmers_rev

pub mod encoders;
mod errors;
pub mod kmer_counter;
pub mod kmer_set;
mod traits;
mod type_mappings;

pub use errors::*;
pub use traits::*;
pub use type_mappings::*;

pub use encoders::KmerEncoder;
pub use kmer_counter::KmerCounter;
pub use kmer_set::KmerSet;

use crate::{data::CheckSequence, prelude::Nucleotides};

/// A k-mer, stored as ASCII bytes in an array of size `MAX_LEN`.
///
/// Since k-mers are short sequences, *Zoe* stores them on the stack using an
/// array, rather than allocating a vector. This struct stores the array as well
/// as the length of the k-mer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Kmer<const MAX_LEN: usize> {
    /// The length of the k-mer
    length: u8,

    /// The bytes in the k-mer, guaranteed to be valid ASCII.
    ///
    /// The k-mer is stored in the first `length` indices. The other indices are
    /// 0.
    buffer: [u8; MAX_LEN],
}

impl<const MAX_LEN: usize> Ord for Kmer<MAX_LEN> {
    #[inline]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // All unoccupied indices are 0, so this implementation is valid
        self.as_ref().cmp(other.as_ref())
    }
}

impl<const MAX_LEN: usize> PartialOrd for Kmer<MAX_LEN> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const MAX_LEN: usize> std::fmt::Display for Kmer<MAX_LEN> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Safety: A Kmer instance is guaranteed to be valid ASCII
        f.write_str(unsafe { std::str::from_utf8_unchecked(self.as_ref()) })
    }
}

impl<const MAX_LEN: usize> Kmer<MAX_LEN> {
    /// Creates a new [`Kmer`] from the provided bases.
    ///
    /// ## Panics
    ///
    /// Panics if `bases` does not contain valid ASCII, or if the length of
    /// `bases` is less than 2 or longer than `MAX_LEN`.
    #[inline]
    #[must_use]
    pub fn new<S: AsRef<[u8]>>(bases: S) -> Self {
        let bases = bases.as_ref();
        let length = bases.len();

        assert!(bases.is_ascii_simd::<16>());
        assert!(bases.len() >= 2);
        let mut buffer = [0; MAX_LEN];
        buffer[..length].copy_from_slice(bases);

        // Safety: buffer is initialized with 0 (valid ASCII) and is filled with
        // values from bases, which was verified to be valid ASCII
        unsafe { Kmer::new_unchecked(length, buffer) }
    }

    /// Creates a new [`Kmer`] from a length and buffer.
    ///
    /// ## Safety
    ///
    /// `buffer` must contain valid ASCII, and length must be between 2 and
    /// `MAX_LEN`.
    ///
    /// ## Validity
    ///
    /// All values in buffer after the first `length` values must be 0.
    #[inline]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub unsafe fn new_unchecked(length: usize, buffer: [u8; MAX_LEN]) -> Self {
        Self {
            // This will not wrap since we require all k-mers to be of length
            // less than 256
            length: length as u8,
            buffer,
        }
    }

    /// Returns the length of the k-mer.
    #[inline]
    #[must_use]
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.length as usize
    }

    /// Returns a vector of the residues in the k-mer.
    #[inline]
    #[must_use]
    pub fn to_vec(&self) -> Vec<u8> {
        self.as_ref().to_vec()
    }

    /// Returns the k-mer as [`Nucleotides`].
    #[inline]
    #[must_use]
    pub fn to_nucleotides(&self) -> Nucleotides {
        Nucleotides(self.to_vec())
    }
}

impl<const MAX_LEN: usize> AsRef<[u8]> for Kmer<MAX_LEN> {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.buffer[..self.len()]
    }
}
