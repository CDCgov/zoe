use crate::{
    kmer::{Kmer, KmerEncoder, KmerError, KmerLen, SupportedKmerLen},
    prelude::Len,
};
use std::{
    collections::{HashSet, hash_set},
    hash::{BuildHasher, RandomState},
    ops::Range,
};

use super::{MismatchNumber, SupportedMismatchNumber};

/// A [`KmerSet`] holds a set of encoded k-mers and provides methods for
/// efficiently using them. For instance, a [`KmerSet`] can be used to find the
/// leftmost or rightmost occurrence of k-mers in a sequence, using
/// [`find_kmers`] or [`find_kmers_rev`]. A [`KmerSet`] can also be constructed
/// using [`insert_from_sequence`] or [`insert_from_sequence_with_variants`].
///
/// [`find_kmers`]: KmerSet::find_kmers
/// [`find_kmers_rev`]: KmerSet::find_kmers_rev
/// [`insert_from_sequence`]: KmerSet::insert_from_sequence
/// [`insert_from_sequence_with_variants`]:
///     KmerSet::insert_from_sequence_with_variants
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct KmerSet<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S = RandomState>
where
    S: BuildHasher, {
    set:     HashSet<E::EncodedKmer, S>,
    encoder: E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> KmerSet<MAX_LEN, E> {
    /// Creates a new [`KmerSet`] with the specified k-mer length.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn new(kmer_length: usize) -> Result<Self, KmerError> {
        Ok(Self {
            set:     HashSet::default(),
            encoder: E::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerSet<MAX_LEN, E, S> {
    /// Creates a new [`KmerSet`] with the specified k-mer length and hasher.
    ///
    /// # Errors
    ///
    /// Returns [`KmerError::InvalidLength`] if `kmer_length` is less than 2 or
    /// greater than `MAX_LEN`.
    #[inline]
    pub fn with_hasher(kmer_length: usize, hasher: S) -> Result<Self, KmerError> {
        Ok(Self {
            set:     HashSet::with_hasher(hasher),
            encoder: E::new(kmer_length)?,
        })
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> KmerSet<MAX_LEN, E, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    /// Get the encoder associated with this [`KmerSet`].
    #[inline]
    pub fn encoder(&self) -> &E {
        &self.encoder
    }

    /// Get the length of the k-mers being stored in the set.
    #[inline]
    pub fn kmer_length(&self) -> usize {
        self.encoder.kmer_length()
    }

    /// Insert an already encoded k-mer into the set. The encoded k-mer must
    /// have been generated using the [`KmerEncoder`] associated with this
    /// [`KmerSet`].
    #[inline]
    pub fn insert_encoded_kmer(&mut self, encoded_kmer: E::EncodedKmer) {
        self.set.insert(encoded_kmer);
    }

    /// Insert a k-mer into the set. The bases and k-mer length are assumed to be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. Consider
    /// [`insert_kmer_checked`] when it is not known whether the bases and k-mer
    /// length will be valid.
    ///
    /// [`insert_kmer_checked`]: KmerSet::insert_kmer_checked
    #[inline]
    pub fn insert_kmer(&mut self, kmer: impl AsRef<[u8]>) {
        self.insert_encoded_kmer(self.encoder.encode_kmer(kmer));
    }

    /// Insert a k-mer into the set. If the bases and k-mer length are not valid
    /// for the [`KmerEncoder`] associated with this [`KmerSet`], then `false`
    /// is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_checked(&mut self, kmer: impl AsRef<[u8]>) -> bool {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer(encoded_kmer);
        true
    }

    /// Return whether an already encoded k-mer is present in the set. The
    /// encoded k-mer must have been generated using the [`KmerEncoder`]
    /// associated with this [`KmerSet`].
    #[inline]
    pub fn contains_encoded(&self, kmer: E::EncodedKmer) -> bool {
        self.set.contains(&kmer)
    }

    /// Return whether a k-mer is present in the set. The bases and k-mer length
    /// are assumed to be valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`]. Consider [`contains_checked`] when it is not known whether
    /// the bases and k-mer length will be valid.
    ///
    /// [`contains_checked`]: KmerSet::contains_checked
    #[inline]
    pub fn contains(&self, kmer: impl AsRef<[u8]>) -> bool {
        self.contains_encoded(self.encoder.encode_kmer(kmer))
    }

    /// Return whether a k-mer is present in the set. If the bases and k-mer
    /// length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`], then `None` is returned.
    #[inline]
    pub fn contains_checked(&self, kmer: impl AsRef<[u8]>) -> Option<bool> {
        Some(self.contains_encoded(self.encoder.encode_kmer_checked(kmer)?))
    }

    /// Iterate over the decoded k-mers in the set.
    #[inline]
    pub fn iter_decoded(&self) -> KmerSetDecodedIter<MAX_LEN, E> {
        KmerSetDecodedIter {
            set_into_iter: self.set.iter(),
            encoder:       &self.encoder,
        }
    }

    /// Iterate over the encoded k-mers in the set.
    #[inline]
    pub fn iter_encoded(&self) -> hash_set::Iter<'_, E::EncodedKmer> {
        self.set.iter()
    }

    /// Insert all k-mers into the set with at most `N` mismatches compared to the
    /// provided, already encoded k-mer. The original k-mer is also inserted.
    /// The encoded k-mer must have been generated using the encoder associated
    /// with this [`KmerSet`].
    ///
    /// `N` must be a supported number of mismatches. See
    /// [`SupportedMismatchNumber`] for more details.
    #[inline]
    pub fn insert_encoded_kmer_with_variants<const N: usize>(&mut self, encoded_kmer: E::EncodedKmer)
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        for variant in self.encoder.get_variants::<N>(encoded_kmer) {
            self.insert_encoded_kmer(variant);
        }
    }

    /// Insert all k-mers into the set with at most N mismatches compared to the
    /// provided k-mer. The original k-mer is also inserted. The bases and k-mer
    /// length are assumed to be valid for the [`KmerEncoder`] associated with
    /// this [`KmerSet`]. Consider [`insert_kmer_with_variants_checked`] when it
    /// is not known whether the bases and k-mer length will be valid.
    ///
    /// [`insert_kmer_with_variants_checked`]:
    ///     KmerSet::insert_kmer_with_variants_checked
    #[inline]
    pub fn insert_kmer_with_variants<const N: usize>(&mut self, kmer: impl AsRef<[u8]>)
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        self.insert_encoded_kmer_with_variants::<N>(self.encoder.encode_kmer(kmer));
    }

    /// Insert all k-mers into the set with at most N mismatches compared to the
    /// provided k-mer. The original k-mer is also inserted. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`], then `false` is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_with_variants_checked<const N: usize>(&mut self, kmer: impl AsRef<[u8]>) -> bool
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer_with_variants::<N>(encoded_kmer);
        true
    }

    /// Insert all k-mers from a sequence into the [`KmerSet`]. The bases in the
    /// sequence must be valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`].
    #[inline]
    pub fn insert_from_sequence(&mut self, seq: impl AsRef<[u8]>) {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer(encoded_kmer);
        }
    }

    /// Insert all k-mers from a sequence into the [`KmerSet`], in addition to
    /// all k-mers with up to N mismatches from those in the sequence. The bases
    /// in the sequence must be valid for the [`KmerEncoder`] associated with
    /// this [`KmerSet`].
    ///
    /// ```
    /// # use zoe::kmer::ThreeBitKmerSet;
    /// let mut set = ThreeBitKmerSet::<8>::new(8).unwrap();
    /// let seq = b"GATAGGGGATTGT";
    /// set.insert_from_sequence_with_variants::<2>(seq);
    /// ```
    #[inline]
    pub fn insert_from_sequence_with_variants<const N: usize>(&mut self, seq: impl AsRef<[u8]>)
    where
        MismatchNumber<N>: SupportedMismatchNumber<MAX_LEN, E>, {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer_with_variants::<N>(encoded_kmer);
        }
    }

    /// Return the indices of the leftmost occurrence of any of the k-mers in
    /// this set within a provided sequence. The bases in the sequence must be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. If no
    /// occurrence is found, then `None` is returned.
    pub fn find_kmers(&self, seq: impl AsRef<[u8]>) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder.iter_from_sequence(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                return Some(i..i + self.kmer_length());
            }
        }
        None
    }

    /// Return the indices of the rightmost occurrence of any of the k-mers in
    /// this set within a provided sequence. The bases in the sequence must be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. If no
    /// occurrence is found, then `None` is returned.
    pub fn find_kmers_rev(&self, seq: impl AsRef<[u8]>) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder.iter_from_sequence_rev(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                let end = seq.as_ref().len() - i;
                return Some(end - self.kmer_length()..end);
            }
        }
        None
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S: BuildHasher> IntoIterator for KmerSet<MAX_LEN, E, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    type Item = Kmer<MAX_LEN>;
    type IntoIter = KmerSetDecodedIntoIter<MAX_LEN, E, S>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            set_into_iter: self.set.into_iter(),
            encoder:       self.encoder,
        }
    }
}

/// An iterator over a [`KmerSet`] yielding decoded k-mers. The iterator
/// consumes the original set.
pub struct KmerSetDecodedIntoIter<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen, {
    pub(crate) set_into_iter: <HashSet<E::EncodedKmer, S> as IntoIterator>::IntoIter,
    pub(crate) encoder:       E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> Iterator for KmerSetDecodedIntoIter<MAX_LEN, E, S>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    type Item = Kmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Kmer<MAX_LEN>> {
        self.set_into_iter.next().map(|x| self.encoder.decode_kmer(x))
    }
}

/// An iterator over a [`KmerSet`] yielding decoded k-mers. The iterator takes
/// the original set by reference.
pub struct KmerSetDecodedIter<'a, const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen, {
    pub(crate) set_into_iter: hash_set::Iter<'a, E::EncodedKmer>,
    pub(crate) encoder:       &'a E,
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> Iterator for KmerSetDecodedIter<'_, MAX_LEN, E>
where
    KmerLen<MAX_LEN, E>: SupportedKmerLen,
{
    type Item = Kmer<MAX_LEN>;

    #[inline]
    fn next(&mut self) -> Option<Kmer<MAX_LEN>> {
        self.set_into_iter.next().map(|x| self.encoder.decode_kmer(*x))
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> Len for KmerSet<MAX_LEN, E> {
    #[inline]
    fn is_empty(&self) -> bool {
        self.set.is_empty()
    }

    #[inline]
    fn len(&self) -> usize {
        self.set.len()
    }
}
