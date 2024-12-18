use crate::kmer::{
    encoder::KmerEncoder, errors::KmerError, kmer_counter::KmerCounter, kmer_set::KmerSet,
    three_bit::encoder::ThreeBitKmerEncoder,
};
use std::{
    collections::HashMap,
    hash::{BuildHasher, RandomState},
    ops::Index,
};

use super::{
    encoder::EncodedKmer,
    int_mappings::{KmerLen, SupportedThreeBitKmerLen},
};

/// A [`KmerCounter`] utilizing [`ThreeBitKmerEncoder`] as its encoder. This
/// allows for `A`, `C`, `G`, `T`, and `N` to all be represented. This counter
/// does not preserve case or the distinction between `T` and `U`. `N` is used
/// as a catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerCounter<const MAX_LEN: usize, S = RandomState>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    map:     HashMap<EncodedKmer<MAX_LEN>, usize, S>,
    encoder: ThreeBitKmerEncoder<MAX_LEN>,
}

impl<const MAX_LEN: usize> ThreeBitKmerCounter<MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    /// Creates a new [`ThreeBitKmerCounter`] with the specified k-mer length.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is 0 or greater than 21.
    #[inline]
    pub fn new(kmer_length: usize) -> Result<Self, KmerError> {
        Ok(Self {
            map:     HashMap::default(),
            encoder: ThreeBitKmerEncoder::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, S> ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    /// Creates a new [`ThreeBitKmerCounter`] with the specified k-mer length and hasher.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is 0 or greater than 21.
    #[inline]
    pub fn with_hasher(kmer_length: usize, hasher: S) -> Result<Self, KmerError> {
        Ok(Self {
            map:     HashMap::with_hasher(hasher),
            encoder: ThreeBitKmerEncoder::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, S: BuildHasher> KmerSet for ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type EncodedKmer = EncodedKmer<MAX_LEN>;
    type Encoder = ThreeBitKmerEncoder<MAX_LEN>;

    /// Get the encoder used for the [`ThreeBitKmerCounter`].
    #[inline]
    fn get_encoder(&self) -> &Self::Encoder {
        &self.encoder
    }

    /// Increment the count of an encoded k-mer in this counter. The encoded
    /// k-mer must have been generated using the encoder associated with this
    /// [`ThreeBitKmerCounter`].
    #[inline]
    fn insert_encoded_kmer(&mut self, encoded_kmer: Self::EncodedKmer) {
        *self.map.entry(encoded_kmer).or_default() += 1;
    }

    /// Return whether an encoded k-mer is present in this counter. The encoded
    /// k-mer must have been generated using the encoder associated with this
    /// [`ThreeBitKmerCounter`].
    #[inline]
    fn contains_encoded(&self, kmer: EncodedKmer<MAX_LEN>) -> bool {
        self.map.contains_key(&kmer)
    }
}

impl<const MAX_LEN: usize, S: BuildHasher> Index<EncodedKmer<MAX_LEN>> for ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Output = usize;

    #[inline]
    fn index(&self, index: EncodedKmer<MAX_LEN>) -> &Self::Output {
        &self.map[&index]
    }
}

impl<const MAX_LEN: usize, S: BuildHasher> KmerCounter for ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    /// Get the count of an encoded k-mer. If the k-mer is not present in the
    /// counter, then `0` is returned. The encoded k-mer must have been generated
    /// using the encoder associated with this [`ThreeBitKmerCounter`].
    #[inline]
    fn get_encoded(&self, kmer: Self::EncodedKmer) -> usize {
        self.map.get(&kmer).copied().unwrap_or_default()
    }
}

impl<const MAX_LEN: usize, S> IntoIterator for ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = (EncodedKmer<MAX_LEN>, usize);
    type IntoIter = ThreeBitKmerCounterIntoIter<MAX_LEN, S>;

    #[inline]
    fn into_iter(self) -> ThreeBitKmerCounterIntoIter<MAX_LEN, S> {
        ThreeBitKmerCounterIntoIter {
            map_into_iter: self.map.into_iter(),
        }
    }
}

/// An iterator over a [`ThreeBitKmerCounter`] yielding decoded k-mers and their
/// counts.
pub struct ThreeBitKmerCounterIntoIter<const MAX_LEN: usize, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    pub(crate) map_into_iter: <HashMap<EncodedKmer<MAX_LEN>, usize, S> as IntoIterator>::IntoIter,
}

impl<const MAX_LEN: usize, S> Iterator for ThreeBitKmerCounterIntoIter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = (EncodedKmer<MAX_LEN>, usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.map_into_iter.next()
    }
}
