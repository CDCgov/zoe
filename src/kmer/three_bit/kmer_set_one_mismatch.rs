use crate::{
    data::types::Uint,
    kmer::{
        encoder::KmerEncoder,
        errors::KmerError,
        kmer_set::KmerSet,
        three_bit::{encoder::ThreeBitKmerEncoder, kmer_set::ThreeBitKmerSetIntoIter},
    },
};
use std::{
    collections::HashSet,
    hash::{BuildHasher, RandomState},
};

use super::encoder::EncodedKmer;

/// A [`KmerSet`] utilizing [`ThreeBitKmerEncoder`] as its encoder, but where
/// all insertions add all k-mers of Hamming distance at most 1. This is useful
/// for fuzzy matching. The encoder allows for `A`, `C`, `G`, `T`, and `N` to
/// all be represented. This set does not preserve case or the distinction
/// between `T` and `U`. `N` is used as a catch-all for bases that are not
/// `ACGTUNacgtun`.
pub struct ThreeBitOneMismatchKmerSet<T: Uint, S = RandomState> {
    set:     HashSet<EncodedKmer<T>, S>,
    encoder: ThreeBitKmerEncoder<T>,
}

impl<T: Uint> ThreeBitOneMismatchKmerSet<T> {
    /// Creates a new [`ThreeBitOneMismatchKmerSet`] with the specified k-mer
    /// length.
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

impl<T: Uint, S> ThreeBitOneMismatchKmerSet<T, S> {
    /// Creates a new [`ThreeBitOneMismatchKmerSet`] with the specified k-mer
    /// length and hasher.
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

impl<T: Uint, S: BuildHasher> KmerSet for ThreeBitOneMismatchKmerSet<T, S> {
    type EncodedKmer = EncodedKmer<T>;
    type Encoder = ThreeBitKmerEncoder<T>;

    /// Get the encoder used for the [`ThreeBitOneMismatchKmerSet`].
    #[inline]
    fn get_encoder(&self) -> &Self::Encoder {
        &self.encoder
    }

    /// Insert an encoded k-mer into the set. This will insert all k-mers with
    /// Hamming distance at most 1 from the provided k-mer. The encoded k-mer
    /// must have been generated using the encoder associated with this
    /// [`ThreeBitOneMismatchKmerSet`].
    #[inline]
    fn insert_encoded_kmer(&mut self, encoded_kmer: EncodedKmer<T>) {
        self.set.insert(encoded_kmer);
        for variant in self.get_encoder().get_variants_one_mismatch(encoded_kmer) {
            self.set.insert(variant);
        }
    }

    /// Return whether an encoded k-mer is present in the set. The encoded k-mer
    /// must have been generated using the encoder associated with this
    /// [`ThreeBitOneMismatchKmerSet`].
    #[inline]
    fn contains_encoded(&self, kmer: EncodedKmer<T>) -> bool {
        self.set.contains(&kmer)
    }
}

impl<T: Uint, S> IntoIterator for ThreeBitOneMismatchKmerSet<T, S> {
    type Item = EncodedKmer<T>;
    type IntoIter = ThreeBitKmerSetIntoIter<T, S>;

    #[inline]
    fn into_iter(self) -> ThreeBitKmerSetIntoIter<T, S> {
        ThreeBitKmerSetIntoIter {
            set_into_iter: self.set.into_iter(),
        }
    }
}
