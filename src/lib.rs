#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]
#![warn(clippy::all, clippy::pedantic)]
#![allow(
    clippy::module_name_repetitions,
    clippy::similar_names,
    clippy::wildcard_imports,
    clippy::enum_glob_use,
    clippy::comparison_chain,
    stable_features
)]
#![cfg_attr(test, feature(test))]
#![feature(portable_simd, const_fn_floating_point_arithmetic, let_chains, lazy_cell, slice_as_chunks)]
#![forbid(incomplete_features)]

pub mod alignment;
/// Composition and consensus functions.
pub mod composition;
pub mod data;
pub mod distance;
pub mod kmer;
/// Mathematical utilities.
pub mod math;
pub mod search;

#[cfg(feature = "rand")]
/// Generate sequences and other data.
pub(crate) mod generate;
/// Iterator utilities.
pub(crate) mod iter_utils;

/// Common structures and traits re-exported
pub mod prelude {
    pub use crate::alignment::PairwiseSequence;
    pub use crate::composition::{AlignmentComposition, CreateConsensus, GcContent, NucleotideCounts, ToBaseCounts};
    pub use crate::data::{
        StdForSequences,
        err::OrFail,
        records::{
            fasta::FastaReader,
            fastq::{FastQ, FastQReader, FastQView, FastQViewMut},
        },
        types::{
            amino_acids::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
            nucleotides::{
                IsValidDNA, Nucleotides, NucleotidesView, NucleotidesViewMut, RecodeDNAStrat, RecodeNucleotides,
                RefineDNAStrat, RetainNucleotides, Translate,
            },
            phred::{QualityScores, QualityScoresView, QualityScoresViewMut, QualityStats},
        },
        view_traits::{DataOwned, DataView, DataViewMut, Len, Slice, SliceMut},
    };
    #[cfg(feature = "rand")]
    pub use crate::generate::rand_sequence;
    pub use crate::kmer::{EncodedKmerCollection, FindKmers, KmerCollectionContains, KmerCounter, KmerEncoder, KmerSet};
    pub use crate::search::{ByteSubstring, ByteSubstringMut};
}

pub(crate) use crate::data::extension::simd;

/// The default SIMD lanes for Zoe.
const DEFAULT_SIMD_LANES: usize = 32;
