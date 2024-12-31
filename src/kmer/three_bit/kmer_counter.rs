use super::{
    encoder::ThreeBitEncodedKmer,
    len_mappings::{KmerLen, SupportedThreeBitKmerLen},
};
use crate::kmer::{
    encoder::{Kmer, KmerEncoder},
    errors::KmerError,
    kmer_counter::KmerCounter,
    kmer_set::KmerSet,
    three_bit::encoder::ThreeBitKmerEncoder,
};
use std::{
    collections::{hash_map, HashMap},
    hash::{BuildHasher, RandomState},
    ops::Index,
};

/// A [`KmerCounter`] utilizing [`ThreeBitKmerEncoder`] as its encoder. This
/// allows for `A`, `C`, `G`, `T`, and `N` to all be represented. This counter
/// does not preserve case or the distinction between `T` and `U`. `N` is used
/// as a catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerCounter<const MAX_LEN: usize, S = RandomState>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    map:     HashMap<ThreeBitEncodedKmer<MAX_LEN>, u64, S>,
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
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
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
    /// Creates a new [`ThreeBitKmerCounter`] with the specified k-mer length
    /// and hasher.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn with_hasher(kmer_length: usize, hasher: S) -> Result<Self, KmerError> {
        Ok(Self {
            map:     HashMap::with_hasher(hasher),
            encoder: ThreeBitKmerEncoder::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, S: BuildHasher> KmerSet<MAX_LEN> for ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type EncodedKmer = ThreeBitEncodedKmer<MAX_LEN>;
    type Encoder = ThreeBitKmerEncoder<MAX_LEN>;
    type DecodedIter<'a>
        = ThreeBitKmerCounterDecodedIter<'a, MAX_LEN>
    where
        S: 'a;
    type EncodedIter<'a>
        = hash_map::Iter<'a, ThreeBitEncodedKmer<MAX_LEN>, u64>
    where
        S: 'a;

    /// Get the encoder used for the [`ThreeBitKmerCounter`].
    #[inline]
    fn encoder(&self) -> &Self::Encoder {
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
    fn contains_encoded(&self, kmer: ThreeBitEncodedKmer<MAX_LEN>) -> bool {
        self.map.contains_key(&kmer)
    }

    /// Iterate over the decoded k-mers and counts in the counter.
    #[inline]
    fn iter_decoded(&self) -> Self::DecodedIter<'_> {
        Self::DecodedIter {
            map_into_iter: self.map.iter(),
            encoder:       &self.encoder,
        }
    }

    /// Iterate over the encoded k-mers and counts in the counter.
    #[inline]
    fn iter_encoded(&self) -> Self::EncodedIter<'_> {
        self.map.iter()
    }
}

impl<const MAX_LEN: usize, S: BuildHasher> Index<ThreeBitEncodedKmer<MAX_LEN>> for ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Output = u64;

    #[inline]
    fn index(&self, index: ThreeBitEncodedKmer<MAX_LEN>) -> &Self::Output {
        &self.map[&index]
    }
}

impl<const MAX_LEN: usize, S: BuildHasher> KmerCounter<MAX_LEN> for ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    /// Get the count of an encoded k-mer. If the k-mer is not present in the
    /// counter, then `0` is returned. The encoded k-mer must have been generated
    /// using the encoder associated with this [`ThreeBitKmerCounter`].
    #[inline]
    fn get_encoded(&self, kmer: Self::EncodedKmer) -> u64 {
        self.map.get(&kmer).copied().unwrap_or_default()
    }
}

impl<const MAX_LEN: usize, S> IntoIterator for ThreeBitKmerCounter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = (Kmer<MAX_LEN>, u64);
    type IntoIter = ThreeBitKmerCounterDecodedIntoIter<MAX_LEN, S>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            map_into_iter: self.map.into_iter(),
            encoder:       self.encoder,
        }
    }
}

/// An iterator over a [`ThreeBitKmerCounter`] yielding decoded k-mers and their
/// counts. The iterator consumes the original counter.
pub struct ThreeBitKmerCounterDecodedIntoIter<const MAX_LEN: usize, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    pub(crate) map_into_iter: <HashMap<ThreeBitEncodedKmer<MAX_LEN>, u64, S> as IntoIterator>::IntoIter,
    pub(crate) encoder:       ThreeBitKmerEncoder<MAX_LEN>,
}

impl<const MAX_LEN: usize, S> Iterator for ThreeBitKmerCounterDecodedIntoIter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = (Kmer<MAX_LEN>, u64);

    #[inline]
    fn next(&mut self) -> Option<(Kmer<MAX_LEN>, u64)> {
        self.map_into_iter.next().map(|(x, c)| (self.encoder.decode_kmer(x), c))
    }
}

/// An iterator over a [`ThreeBitKmerCounter`] yielding decoded k-mers and their
/// counts. The iterator takes the original counter by reference.
pub struct ThreeBitKmerCounterDecodedIter<'a, const MAX_LEN: usize>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    pub(crate) map_into_iter: hash_map::Iter<'a, ThreeBitEncodedKmer<MAX_LEN>, u64>,
    pub(crate) encoder:       &'a ThreeBitKmerEncoder<MAX_LEN>,
}

impl<const MAX_LEN: usize> Iterator for ThreeBitKmerCounterDecodedIter<'_, MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = (Kmer<MAX_LEN>, u64);

    #[inline]
    fn next(&mut self) -> Option<(Kmer<MAX_LEN>, u64)> {
        self.map_into_iter.next().map(|(x, &c)| (self.encoder.decode_kmer(*x), c))
    }
}
