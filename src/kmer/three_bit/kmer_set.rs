use crate::{
    data::types::Uint,
    kmer::{encoder::KmerEncoder, errors::KmerError, kmer_set::KmerSet, three_bit::encoder::ThreeBitKmerEncoder},
};
use std::{
    collections::HashSet,
    hash::{BuildHasher, RandomState},
};

use super::encoder::EncodedKmer;

/// A [`KmerSet`] utilizing [`ThreeBitKmerEncoder`] as its encoder. This allows
/// for `A`, `C`, `G`, `T`, and `N` to all be represented. This set does not
/// preserve case or the distinction between `T` and `U`. `N` is used as a
/// catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerSet<T: Uint, S = RandomState> {
    set:     HashSet<EncodedKmer<T>, S>,
    encoder: ThreeBitKmerEncoder<T>,
}

impl<T: Uint> ThreeBitKmerSet<T> {
    /// Creates a new [`ThreeBitKmerSet`] with the specified k-mer length.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is 0 or greater
    /// than 21.
    #[inline]
    pub fn new(kmer_length: usize) -> Result<Self, KmerError> {
        Ok(Self {
            set:     HashSet::default(),
            encoder: ThreeBitKmerEncoder::new(kmer_length)?,
        })
    }
}

impl<T: Uint, S> ThreeBitKmerSet<T, S> {
    /// Creates a new [`ThreeBitKmerSet`] with the specified k-mer length and
    /// hasher.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is 0 or greater
    /// than 21.
    #[inline]
    pub fn with_hasher(kmer_length: usize, hasher: S) -> Result<Self, KmerError> {
        Ok(Self {
            set:     HashSet::with_hasher(hasher),
            encoder: ThreeBitKmerEncoder::new(kmer_length)?,
        })
    }
}

impl<T: Uint, S: BuildHasher> KmerSet for ThreeBitKmerSet<T, S> {
    type EncodedKmer = EncodedKmer<T>;
    type Encoder = ThreeBitKmerEncoder<T>;

    /// Get the encoder used for the [`ThreeBitKmerSet`].
    #[inline]
    fn get_encoder(&self) -> &ThreeBitKmerEncoder<T> {
        &self.encoder
    }

    /// Insert an encoded k-mer into the set. The encoded k-mer must have been
    /// generated using the encoder associated with this [`ThreeBitKmerSet`].
    #[inline]
    fn insert_encoded_kmer(&mut self, encoded_kmer: EncodedKmer<T>) {
        self.set.insert(encoded_kmer);
    }

    /// Return whether an encoded k-mer is present in the set. The encoded k-mer
    /// must have been generated using the encoder associated with this
    /// [`ThreeBitKmerSet`].
    #[inline]
    fn contains_encoded(&self, kmer: EncodedKmer<T>) -> bool {
        self.set.contains(&kmer)
    }
}

impl<T: Uint, S> IntoIterator for ThreeBitKmerSet<T, S> {
    type Item = EncodedKmer<T>;
    type IntoIter = ThreeBitKmerSetIntoIter<T, S>;

    #[inline]
    fn into_iter(self) -> ThreeBitKmerSetIntoIter<T, S> {
        ThreeBitKmerSetIntoIter {
            set_into_iter: self.set.into_iter(),
        }
    }
}

/// An iterator over a [`ThreeBitKmerSet`] or [`ThreeBitOneMismatchKmerSet`]
/// yielding decoded k-mers.
///
/// [`ThreeBitOneMismatchKmerSet`]:
///     super::kmer_set_one_mismatch::ThreeBitOneMismatchKmerSet
pub struct ThreeBitKmerSetIntoIter<T, S> {
    pub(crate) set_into_iter: <HashSet<EncodedKmer<T>, S> as IntoIterator>::IntoIter,
}

impl<T: Uint, S> Iterator for ThreeBitKmerSetIntoIter<T, S> {
    type Item = EncodedKmer<T>;

    #[inline]
    fn next(&mut self) -> Option<EncodedKmer<T>> {
        self.set_into_iter.next()
    }
}
