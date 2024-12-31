use super::{
    encoder::ThreeBitEncodedKmer,
    len_mappings::{KmerLen, SupportedThreeBitKmerLen},
};
use crate::kmer::{
    encoder::{Kmer, KmerEncoder},
    errors::KmerError,
    kmer_set::KmerSet,
    three_bit::encoder::ThreeBitKmerEncoder,
};
use std::{
    collections::{hash_set, HashSet},
    hash::{BuildHasher, RandomState},
};

/// A [`KmerSet`] utilizing [`ThreeBitKmerEncoder`] as its encoder. This allows
/// for `A`, `C`, `G`, `T`, and `N` to all be represented. This set does not
/// preserve case or the distinction between `T` and `U`. `N` is used as a
/// catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerSet<const MAX_LEN: usize, S = RandomState>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    set:     HashSet<ThreeBitEncodedKmer<MAX_LEN>, S>,
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

impl<const MAX_LEN: usize, S> ThreeBitKmerSet<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    /// Creates a new [`ThreeBitKmerSet`] with the specified k-mer length and
    /// hasher.
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

impl<const MAX_LEN: usize, S: BuildHasher> KmerSet<MAX_LEN> for ThreeBitKmerSet<MAX_LEN, S>
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

    /// Get the encoder used for the [`ThreeBitKmerSet`].
    #[inline]
    fn encoder(&self) -> &Self::Encoder {
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

impl<const MAX_LEN: usize, S> IntoIterator for ThreeBitKmerSet<MAX_LEN, S>
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

/// An iterator over a [`ThreeBitKmerSet`] or [`ThreeBitOneMismatchKmerSet`]
/// yielding decoded k-mers. The iterator consumes the original set.
///
/// [`ThreeBitOneMismatchKmerSet`]:
///     super::kmer_set_one_mismatch::ThreeBitOneMismatchKmerSet
pub struct ThreeBitKmerSetDecodedIntoIter<const MAX_LEN: usize, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    pub(crate) set_into_iter: <HashSet<ThreeBitEncodedKmer<MAX_LEN>, S> as IntoIterator>::IntoIter,
    pub(crate) encoder:       ThreeBitKmerEncoder<MAX_LEN>,
}

impl<const MAX_LEN: usize, S> Iterator for ThreeBitKmerSetDecodedIntoIter<MAX_LEN, S>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = Kmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Kmer<MAX_LEN>> {
        self.set_into_iter.next().map(|x| self.encoder.decode_kmer(x))
    }
}

/// An iterator over a [`ThreeBitKmerSet`] or [`ThreeBitOneMismatchKmerSet`]
/// yielding decoded k-mers. The iterator takes the original set by reference.
///
/// [`ThreeBitOneMismatchKmerSet`]:
///     super::kmer_set_one_mismatch::ThreeBitOneMismatchKmerSet
pub struct ThreeBitKmerSetDecodedIter<'a, const MAX_LEN: usize>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen, {
    pub(crate) set_into_iter: hash_set::Iter<'a, ThreeBitEncodedKmer<MAX_LEN>>,
    pub(crate) encoder:       &'a ThreeBitKmerEncoder<MAX_LEN>,
}

impl<const MAX_LEN: usize> Iterator for ThreeBitKmerSetDecodedIter<'_, MAX_LEN>
where
    KmerLen<MAX_LEN>: SupportedThreeBitKmerLen,
{
    type Item = Kmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Kmer<MAX_LEN>> {
        self.set_into_iter.next().map(|x| self.encoder.decode_kmer(*x))
    }
}
