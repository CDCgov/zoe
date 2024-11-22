use crate::kmer::{encoder::KmerEncoder, errors::KmerError, kmer_set::KmerSet, three_bit::encoder::ThreeBitKmerEncoder};
use std::{
    collections::HashSet,
    hash::{BuildHasher, RandomState},
};

/// A [`KmerSet`] utilizing [`ThreeBitKmerEncoder`] as its encoder. This allows
/// for `A`, `C`, `G`, `T`, and `N` to all be represented. This set does not
/// preserve case or the distinction between `T` and `U`. `N` is used as a
/// catch-all for bases that are not `ACGTUNacgtun`.
pub struct ThreeBitKmerSet<S = RandomState> {
    set:     HashSet<u64, S>,
    encoder: ThreeBitKmerEncoder,
}

impl ThreeBitKmerSet {
    /// Creates a new [`ThreeBitKmerSet`] with the specified kmer length.
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

impl<S> ThreeBitKmerSet<S> {
    /// Creates a new [`ThreeBitKmerSet`] with the specified kmer length and
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

impl<S: BuildHasher> KmerSet for ThreeBitKmerSet<S> {
    type Kmer = u64;
    type Encoder = ThreeBitKmerEncoder;

    /// Get the encoder used for the [`ThreeBitKmerSet`].
    #[inline]
    fn get_encoder(&self) -> &ThreeBitKmerEncoder {
        &self.encoder
    }

    /// Insert an encoded kmer into the set. The encoded kmer must have been
    /// generated using the encoder associated with this [`ThreeBitKmerSet`].
    #[inline]
    fn insert_encoded_kmer(&mut self, encoded_kmer: u64) {
        self.set.insert(encoded_kmer);
    }

    /// Return whether an encoded kmer is present in the set. The encoded kmer
    /// must have been generated using the encoder associated with this
    /// [`ThreeBitKmerSet`].
    #[inline]
    fn contains_encoded(&self, kmer: u64) -> bool {
        self.set.contains(&kmer)
    }
}

impl<S> IntoIterator for ThreeBitKmerSet<S> {
    type Item = Vec<u8>;
    type IntoIter = ThreeBitKmerSetIntoIter<S>;

    #[inline]
    fn into_iter(self) -> ThreeBitKmerSetIntoIter<S> {
        ThreeBitKmerSetIntoIter {
            encoder:       self.encoder,
            set_into_iter: self.set.into_iter(),
        }
    }
}

/// An iterator over a [`ThreeBitKmerSet`] or [`ThreeBitOneMismatchKmerSet`]
/// yielding decoded kmers.
///
/// [`ThreeBitOneMismatchKmerSet`]:
///     super::kmer_set_one_mismatch::ThreeBitOneMismatchKmerSet
pub struct ThreeBitKmerSetIntoIter<S> {
    pub(crate) encoder:       ThreeBitKmerEncoder,
    pub(crate) set_into_iter: <HashSet<u64, S> as IntoIterator>::IntoIter,
}

impl<S> Iterator for ThreeBitKmerSetIntoIter<S> {
    type Item = Vec<u8>;

    #[inline]
    fn next(&mut self) -> Option<Vec<u8>> {
        Some(self.encoder.decode_kmer(self.set_into_iter.next()?))
    }
}
