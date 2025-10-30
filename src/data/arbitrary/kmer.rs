use crate::{
    data::arbitrary::{AsciiByte, SameSizeVecs, VecBounded},
    kmer::{
        EncodedKmerCollection, Kmer, KmerEncoder, SupportedKmerLen,
        encoders::three_bit::{ThreeBitKmerEncoder, ThreeBitKmerLen, ThreeBitKmerSet},
    },
};
use arbitrary::{Arbitrary, Result, Unstructured};
use std::hash::BuildHasher;

impl<'a, const MAX_LEN: usize> Arbitrary<'a> for Kmer<MAX_LEN> {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let kmer = VecBounded::<2, MAX_LEN, AsciiByte>::arbitrary(u)?
            .0
            .into_iter()
            .map(|x| x.0)
            .collect::<Vec<_>>();

        Ok(Kmer::new(kmer))
    }
}

impl<'a, const MAX_LEN: usize> Arbitrary<'a> for ThreeBitKmerEncoder<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(ThreeBitKmerEncoder::new(u.int_in_range(2..=MAX_LEN)?).unwrap())
    }
}

impl<'a, const MAX_LEN: usize, S: BuildHasher + Default> Arbitrary<'a> for ThreeBitKmerSet<MAX_LEN, S>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let mut set = ThreeBitKmerSet::with_hasher(u.int_in_range(2..=MAX_LEN)?, S::default()).unwrap();
        let kmers = SameSizeVecs::<2, MAX_LEN, AsciiByte>::arbitrary_with_len(u, set.kmer_length())?.vecs;
        set.insert_from_iter(kmers.into_iter().map(|x| x.0.into_iter().map(|x| x.0).collect::<Vec<_>>()));
        Ok(set)
    }
}
