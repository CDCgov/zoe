//! ## K-mer processing, searching, and counting.
//!
//! This module provides structs and methods for handling common k-mer
//! operations efficiently using integer encodings. A k-mer is a short
//! subsequence of nucleotides, which are represented in *Zoe* by [`Kmer`].
//!
//! *Zoe* has two structs to store k-mers:
//! * A [`KmerSet`], which is a [`HashSet`] for k-mers
//! * A [`KmerCounter`] to store k-mers and their counts
//!
//! Several convenience methods are provided, such as:
//! * [`insert_from_sequence`], which quickly inserts/counts all overlapping
//!   k-mers from a sequence
//! * [`find_in_seq`] and [`find_in_seq_rev`], which search a sequence for any
//!   of the stored k-mers. The [`FindKmers`] trait provides similar methods
//!   [`find_kmers`] and [`find_kmers_rev`] directly on the sequence, and
//!   accepting the kmer collection as an argument.
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
//! encode and decode the k-mer.
//!
//! </div>
//!
//! ## Examples
//!
//! Count the 3-mers in a sequence:
//! ```
//! # use zoe::{kmer::ThreeBitKmerCounter, prelude::*};
//! let sequence = b"GGCCACCAAGGCCA";
//! let mut kmer_counter = ThreeBitKmerCounter::<3>::new(3).unwrap();
//! kmer_counter.insert_from_sequence(sequence);
//! for (kmer, count) in kmer_counter {
//!     println!("{kmer}\t{count}");
//! }
//! ```
//!
//! Search for the 17-mers of a primer within a sequence, with up to one
//! mismatch:
//! ```
//! # use zoe::{kmer::ThreeBitKmerSet, prelude::*};
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
//! # use zoe::{kmer::{ThreeBitKmerSet, FindKmers}, prelude::*};
//! let primer: Nucleotides = b"TGATAGTTTTAGAGTTAGGTAG".into();
//! let sequence: Nucleotides = b"TGCCCGTAACGTACAGTTTTACAGTTAGGTACCC".into();
//! let mut kmer_set = ThreeBitKmerSet::<17>::new(17).unwrap();
//! kmer_set.insert_from_sequence_with_variants::<1>(primer);
//! let kmer_pos = sequence.find_kmers(&kmer_set);
//! assert_eq!(kmer_pos, Some(14..31));
//! ```
//!
//! [`HashSet`]: std::collections::HashSet
//! [`Display`]: std::fmt::Display
//! [`Kmer`]: Kmer
//! [`KmerEncoder`]: KmerEncoder
//! [`KmerSet`]: KmerSet
//! [`KmerCounter`]: KmerCounter
//! [`ThreeBitKmerEncoder`]: ThreeBitKmerEncoder
//! [`ThreeBitEncodedKmer`]: ThreeBitEncodedKmer
//! [`ThreeBitKmerSet`]: ThreeBitKmerSet
//! [`ThreeBitKmerCounter`]: ThreeBitKmerCounter
//! [`SupportedKmerLen`]: SupportedKmerLen
//! [`insert_from_sequence`]: KmerSet::insert_from_sequence
//! [`find_in_seq`]: KmerCollectionContains::find_in_seq
//! [`find_in_seq_rev`]: KmerCollectionContains::find_in_seq_rev
//! [`FindKmers`]: FindKmers
//! [`find_kmers`]: FindKmers::find_kmers
//! [`find_kmers_rev`]: FindKmers::find_kmers_rev

mod encoder;
mod errors;
mod kmer_counter;
mod kmer_set;
mod three_bit;
mod traits;
mod type_mappings;

pub use encoder::*;
pub use errors::*;
pub use kmer_counter::*;
pub use kmer_set::*;
pub use three_bit::*;
pub use traits::*;
pub use type_mappings::*;
