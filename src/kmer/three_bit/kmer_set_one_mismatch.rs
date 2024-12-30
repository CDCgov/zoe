use crate::kmer::{
    encoder::{Kmer, KmerEncoder},
    errors::KmerError,
    kmer_set::KmerSet,
    three_bit::{encoder::ThreeBitKmerEncoder, kmer_set::ThreeBitKmerSetDecodedIter},
};
use std::{
    collections::{hash_set, HashSet},
    hash::{BuildHasher, RandomState},
};

use super::{
    encoder::ThreeBitEncodedKmer,
    kmer_set::ThreeBitKmerSetDecodedIntoIter,
    len_mappings::{KmerLen, SupportedThreeBitKmerLen},
};

/// A [`KmerSet`] utilizing [`ThreeBitKmerEncoder`] as its encoder, but where
/// all insertions add all k-mers of Hamming distance at most 1. This is useful
/// for fuzzy matching. The encoder allows for `A`, `C`, `G`, `T`, and `N` to
/// all be represented. This set does not preserve case or the distinction
/// between `T` and `U`. `N` is used as a catch-all for bases that are not
/// `ACGTUNacgtun`.
pub struct ThreeBitOneMismatchKmerSet<const MAX_LEN: usize, S = RandomState>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    set:     HashSet<ThreeBitEncodedKmer<MAX_LEN>, S>,
    encoder: ThreeBitKmerEncoder<MAX_LEN>,
}

impl<const MAX_LEN: usize> ThreeBitOneMismatchKmerSet<MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    /// Creates a new [`ThreeBitOneMismatchKmerSet`] with the specified k-mer
    /// length.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn new(kmer_length: usize) -> Result<Self, KmerError> {
        Ok(Self {
            set:     HashSet::default(),
            encoder: ThreeBitKmerEncoder::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, S> ThreeBitOneMismatchKmerSet<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    /// Creates a new [`ThreeBitOneMismatchKmerSet`] with the specified k-mer
    /// length and hasher.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn with_hasher(kmer_length: usize, hasher: S) -> Result<Self, KmerError> {
        Ok(Self {
            set:     HashSet::with_hasher(hasher),
            encoder: ThreeBitKmerEncoder::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, S: BuildHasher> KmerSet<MAX_LEN> for ThreeBitOneMismatchKmerSet<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type EncodedKmer = ThreeBitEncodedKmer<MAX_LEN>;
    type Encoder = ThreeBitKmerEncoder<MAX_LEN>;
    type DecodedIter<'a>
        = ThreeBitKmerSetDecodedIter<'a, MAX_LEN>
    where
        S: 'a;
    type EncodedIter<'a>
        = hash_set::Iter<'a, ThreeBitEncodedKmer<MAX_LEN>>
    where
        S: 'a;

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
    fn insert_encoded_kmer(&mut self, encoded_kmer: Self::EncodedKmer) {
        self.set.insert(encoded_kmer);
        for variant in self.get_encoder().get_variants_one_mismatch(encoded_kmer) {
            self.set.insert(variant);
        }
    }

    /// Return whether an encoded k-mer is present in the set. The encoded k-mer
    /// must have been generated using the encoder associated with this
    /// [`ThreeBitOneMismatchKmerSet`].
    #[inline]
    fn contains_encoded(&self, kmer: Self::EncodedKmer) -> bool {
        self.set.contains(&kmer)
    }

    /// Iterate over the decoded k-mers in the set.
    #[inline]
    fn iter_decoded(&self) -> Self::DecodedIter<'_> {
        Self::DecodedIter {
            set_into_iter: self.set.iter(),
            encoder:       &self.encoder,
        }
    }

    /// Iterate over the encoded k-mers in the set.
    #[inline]
    fn iter_encoded(&self) -> Self::EncodedIter<'_> {
        self.set.iter()
    }
}

impl<const MAX_LEN: usize, S> IntoIterator for ThreeBitOneMismatchKmerSet<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = Kmer<MAX_LEN>;
    type IntoIter = ThreeBitKmerSetDecodedIntoIter<MAX_LEN, S>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            set_into_iter: self.set.into_iter(),
            encoder:       self.encoder,
        }
    }
}
