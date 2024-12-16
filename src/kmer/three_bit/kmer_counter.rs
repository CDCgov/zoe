use crate::kmer::{
    encoder::KmerEncoder, errors::KmerError, kmer_counter::KmerCounter, kmer_set::KmerSet,
    three_bit::encoder::ThreeBitKmerEncoder,
};
use std::{
    collections::HashMap,
    hash::{BuildHasher, RandomState},
    ops::Index,
};

use super::encoder::EncodedKmer;

/// A [`KmerCounter`] utilizing [`ThreeBitKmerEncoder`] as its encoder. This
/// allows for `A`, `C`, `G`, `T`, and `N` to all be represented. This counter
/// does not preserve case or the distinction between `T` and `U`. `N` is used
/// as a catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerCounter<S = RandomState> {
    map:     HashMap<EncodedKmer, usize, S>,
    encoder: ThreeBitKmerEncoder,
}

impl ThreeBitKmerCounter {
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

impl<S> ThreeBitKmerCounter<S> {
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

impl<S: BuildHasher> KmerSet for ThreeBitKmerCounter<S> {
    type Kmer = EncodedKmer;
    type Encoder = ThreeBitKmerEncoder;

    /// Get the encoder used for the [`ThreeBitKmerCounter`].
    #[inline]
    fn get_encoder(&self) -> &ThreeBitKmerEncoder {
        &self.encoder
    }

    /// Increment the count of an encoded k-mer in this counter. The encoded
    /// k-mer must have been generated using the encoder associated with this
    /// [`ThreeBitKmerCounter`].
    #[inline]
    fn insert_encoded_kmer(&mut self, encoded_kmer: EncodedKmer) {
        *self.map.entry(encoded_kmer).or_default() += 1;
    }

    /// Return whether an encoded k-mer is present in this counter. The encoded
    /// k-mer must have been generated using the encoder associated with this
    /// [`ThreeBitKmerCounter`].
    #[inline]
    fn contains_encoded(&self, kmer: EncodedKmer) -> bool {
        self.map.contains_key(&kmer)
    }
}

impl<S: BuildHasher> Index<EncodedKmer> for ThreeBitKmerCounter<S> {
    type Output = usize;

    #[inline]
    fn index(&self, index: EncodedKmer) -> &Self::Output {
        &self.map[&index]
    }
}

impl<S: BuildHasher> KmerCounter for ThreeBitKmerCounter<S> {
    /// Get the count of an encoded k-mer. If the k-mer is not present in the
    /// counter, then `0` is returned. The encoded k-mer must have been generated
    /// using the encoder associated with this [`ThreeBitKmerCounter`].
    #[inline]
    fn get_encoded(&self, kmer: Self::Kmer) -> usize {
        self.map.get(&kmer).copied().unwrap_or_default()
    }
}

impl<S> IntoIterator for ThreeBitKmerCounter<S> {
    type Item = (EncodedKmer, usize);
    type IntoIter = ThreeBitKmerCounterIntoIter<S>;

    #[inline]
    fn into_iter(self) -> ThreeBitKmerCounterIntoIter<S> {
        ThreeBitKmerCounterIntoIter {
            map_into_iter: self.map.into_iter(),
        }
    }
}

/// An iterator over a [`ThreeBitKmerCounter`] yielding decoded k-mers and their
/// counts.
pub struct ThreeBitKmerCounterIntoIter<S> {
    pub(crate) map_into_iter: <HashMap<EncodedKmer, usize, S> as IntoIterator>::IntoIter,
}

impl<S> Iterator for ThreeBitKmerCounterIntoIter<S> {
    type Item = (EncodedKmer, usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.map_into_iter.next()
    }
}
