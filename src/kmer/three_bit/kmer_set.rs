use super::{
    encoder::EncodedKmer,
    int_mappings::{KmerLen, SupportedThreeBitKmerLen},
};
use crate::kmer::{encoder::KmerEncoder, errors::KmerError, kmer_set::KmerSet, three_bit::encoder::ThreeBitKmerEncoder};
use std::{
    collections::HashSet,
    hash::{BuildHasher, RandomState},
};

/// A [`KmerSet`] utilizing [`ThreeBitKmerEncoder`] as its encoder. This allows
/// for `A`, `C`, `G`, `T`, and `N` to all be represented. This set does not
/// preserve case or the distinction between `T` and `U`. `N` is used as a
/// catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerSet<const MAX_LEN: usize, S = RandomState>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    set:     HashSet<EncodedKmer<MAX_LEN>, S>,
    encoder: ThreeBitKmerEncoder<MAX_LEN>,
}

impl<const MAX_LEN: usize> ThreeBitKmerSet<MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
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

impl<const MAX_LEN: usize, S> ThreeBitKmerSet<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
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

impl<const MAX_LEN: usize, S: BuildHasher> KmerSet for ThreeBitKmerSet<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type EncodedKmer = EncodedKmer<MAX_LEN>;
    type Encoder = ThreeBitKmerEncoder<MAX_LEN>;

    /// Get the encoder used for the [`ThreeBitKmerSet`].
    #[inline]
    fn get_encoder(&self) -> &Self::Encoder {
        &self.encoder
    }

    /// Insert an encoded k-mer into the set. The encoded k-mer must have been
    /// generated using the encoder associated with this [`ThreeBitKmerSet`].
    #[inline]
    fn insert_encoded_kmer(&mut self, encoded_kmer: Self::EncodedKmer) {
        self.set.insert(encoded_kmer);
    }

    /// Return whether an encoded k-mer is present in the set. The encoded k-mer
    /// must have been generated using the encoder associated with this
    /// [`ThreeBitKmerSet`].
    #[inline]
    fn contains_encoded(&self, kmer: Self::EncodedKmer) -> bool {
        self.set.contains(&kmer)
    }
}

impl<const MAX_LEN: usize, S> IntoIterator for ThreeBitKmerSet<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = EncodedKmer<MAX_LEN>;
    type IntoIter = ThreeBitKmerSetIntoIter<MAX_LEN, S>;

    #[inline]
    fn into_iter(self) -> ThreeBitKmerSetIntoIter<MAX_LEN, S> {
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
pub struct ThreeBitKmerSetIntoIter<const MAX_LEN: usize, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    pub(crate) set_into_iter: <HashSet<EncodedKmer<MAX_LEN>, S> as IntoIterator>::IntoIter,
}

impl<const MAX_LEN: usize, S> Iterator for ThreeBitKmerSetIntoIter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = EncodedKmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<EncodedKmer<MAX_LEN>> {
        self.set_into_iter.next()
    }
}
