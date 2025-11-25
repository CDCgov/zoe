//! Arbitrary implementations and specification structs for k-mer related types
//! in the [`kmer`](crate::kmer) module.

use crate::{
    data::arbitrary::{ArbitrarySpecs, ByteSet, ByteSpecs, Case, VecSpecs},
    kmer::{
        Kmer, KmerEncoder, SupportedKmerLen,
        encoders::three_bit::{ThreeBitKmerEncoder, ThreeBitKmerLen, ThreeBitKmerSet},
    },
    prelude::Len,
};
use arbitrary::{Arbitrary, Result, Unstructured};
use std::hash::BuildHasher;

impl<'a, const MAX_LEN: usize> Arbitrary<'a> for Kmer<MAX_LEN> {
    /// Generates an arbitrary [`Kmer`] from the given unstructured data.
    ///
    /// This ensures that the k-mer contains a buffer and length which agree
    /// with each other, and that the bytes are ASCII.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let kmer_specs = VecSpecs {
            element_specs: ByteSpecs {
                set:  ByteSet::Ascii,
                case: Case::Any,
            },
            min_len:       2,
            max_len:       MAX_LEN,
            len:           None,
        };

        Ok(Kmer::new(kmer_specs.make_arbitrary(u)?))
    }
}

impl<'a, const MAX_LEN: usize> Arbitrary<'a> for ThreeBitKmerEncoder<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    /// Generates an arbitrary [`ThreeBitKmerEncoder`] from the given
    /// unstructured data.
    ///
    /// This ensures that the encoder contains valid state; the only arbitrary
    /// aspect is the k-mer length chosen.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(ThreeBitKmerEncoder::new(u.int_in_range(2..=MAX_LEN)?).unwrap())
    }
}

impl<'a, const MAX_LEN: usize, S: BuildHasher + Default> Arbitrary<'a> for ThreeBitKmerSet<MAX_LEN, S>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    /// Generates an arbitrary [`ThreeBitKmerSet`] from the given unstructured
    /// data.
    ///
    /// This ensures that the set contains encoded k-mers of the proper size.
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let mut set = ThreeBitKmerSet::with_hasher(u.int_in_range(2..=MAX_LEN)?, S::default()).unwrap();

        let kmer_specs = VecSpecs {
            element_specs: ByteSpecs {
                set:  ByteSet::Ascii,
                case: Case::Any,
            },
            min_len:       0,
            len:           Some(set.len()),
            max_len:       usize::MAX,
        };

        set.insert_from_iter(kmer_specs.make_arbitrary_iter(u).filter_map(std::result::Result::ok));

        Ok(set)
    }
}
