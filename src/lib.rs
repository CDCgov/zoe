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
#![feature(portable_simd)]
// Stabilized (in at least one channel) but left here in case it helps backwards
// compatibility. Can be removed in the future.
#![feature(const_fn_floating_point_arithmetic, let_chains, lazy_cell, slice_as_chunks)]
#![forbid(incomplete_features)]
#![deny(rustdoc::broken_intra_doc_links)]

/// Mathematical utilities.
pub mod math;

pub mod alignment;
/// Composition and consensus functions.
pub mod composition;
pub mod data;
pub mod distance;
pub mod kmer;

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
                CheckNucleotides, IsValidDNA, Nucleotides, NucleotidesView, NucleotidesViewMut, RecodeDNAStrat,
                RecodeNucleotides, RefineDNAStrat, RetainNucleotides, Translate,
            },
            phred::{QualityScores, QualityScoresView, QualityScoresViewMut, QualityStats},
        },
        views::{DataOwned, DataView, DataViewMut, Len, Restrict, Slice, SliceMut},
    };
    #[cfg(feature = "rand")]
    pub use crate::generate::rand_sequence;
    pub use crate::kmer::{EncodedKmerCollection, FindKmers, KmerCollectionContains, KmerCounter, KmerEncoder, KmerSet};
    pub use crate::search::{ByteSubstring, ByteSubstringMut};
}

pub(crate) use crate::data::extension::simd;

mod private {
    use crate::{
        data::{
            cigar::{Cigar, CigarView, CigarViewMut},
            sam::{SamData, SamDataView, SamDataViewMut},
        },
        prelude::*,
    };
    use std::{
        hash::BuildHasher,
        simd::{LaneCount, Simd, SimdElement, SupportedLaneCount},
    };

    macro_rules! sealed {
        ($($t:ty),* $(,)?) => { $( impl Sealed for $t {} )* };
    }

    pub trait Sealed {}
    sealed!(String, &String, &str, str);
    impl<T> Sealed for Vec<T> {}
    impl<T> Sealed for &[T] {}
    impl<T> Sealed for [T] {}
    impl<T, const N: usize> Sealed for &[T; N] {}
    impl<T, const N: usize> Sealed for [T; N] {}
    impl Sealed for crate::search::RangeSearch<'_> {}
    impl<T: SimdElement, const N: usize> Sealed for Simd<T, N> where LaneCount<N>: SupportedLaneCount {}
    impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> Sealed for KmerCounter<MAX_LEN, E, S> {}
    impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> Sealed for KmerSet<MAX_LEN, E, S> {}
    sealed!(
        AminoAcids,
        AminoAcidsView<'_>,
        AminoAcidsViewMut<'_>,
        Nucleotides,
        NucleotidesView<'_>,
        NucleotidesViewMut<'_>,
        FastQ,
        FastQView<'_>,
        FastQViewMut<'_>,
        QualityScores,
        QualityScoresView<'_>,
        QualityScoresViewMut<'_>,
        Cigar,
        CigarView<'_>,
        CigarViewMut<'_>,
        SamData,
        SamDataView<'_>,
        SamDataViewMut<'_>,
    );
    sealed!(f32, f64, u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize);
}

/// The default SIMD lanes for Zoe.
const DEFAULT_SIMD_LANES: usize = 32;
