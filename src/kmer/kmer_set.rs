use crate::kmer::{Kmer, KmerEncoder, KmerError, KmerLen, SupportedKmerLen};
use std::{
    collections::{HashSet, hash_set},
    hash::{BuildHasher, RandomState},
    ops::Range,
};

pub struct KmerSet<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S = RandomState> {
    set:         HashSet<E::EncodedKmer, S>,
    pub encoder: E,
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

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> KmerSet<MAX_LEN, E, S> {
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
    /// Get the length of the k-mers being stored in the set.
    #[inline]
    pub fn get_kmer_length(&self) -> usize {
        self.encoder.get_kmer_length()
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
    pub fn insert_kmer<Q: AsRef<[u8]>>(&mut self, kmer: Q) {
        self.insert_encoded_kmer(self.encoder.encode_kmer(kmer));
    }

    /// Insert a k-mer into the set. If the bases and k-mer length are not valid
    /// for the [`KmerEncoder`] associated with this [`KmerSet`], then `false`
    /// is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_checked<Q: AsRef<[u8]>>(&mut self, kmer: Q) -> bool {
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
    pub fn contains<Q: AsRef<[u8]>>(&self, kmer: Q) -> bool {
        self.contains_encoded(self.encoder.encode_kmer(kmer))
    }

    /// Return whether a k-mer is present in the set. If the bases and k-mer
    /// length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`], then `None` is returned.
    #[inline]
    pub fn contains_checked<Q: AsRef<[u8]>>(&self, kmer: Q) -> Option<bool> {
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

    /// Insert all k-mers into the set with at most one mismatch compared to the
    /// provided, already encoded k-mer. The original k-mer is also inserted.
    /// The encoded k-mer must have been generated using the encoder associated
    /// with this [`KmerSet`].
    #[inline]
    pub fn insert_encoded_kmer_one_mismatch(&mut self, encoded_kmer: E::EncodedKmer) {
        self.insert_encoded_kmer(encoded_kmer);
        for variant in self.encoder.get_variants_one_mismatch(encoded_kmer) {
            self.insert_encoded_kmer(variant);
        }
    }

    /// Insert all k-mers into the set with at most one mismatch compared to the
    /// provided k-mer. The original k-mer is also inserted. The bases and k-mer
    /// length are assumed to be valid for the [`KmerEncoder`] associated with
    /// this [`KmerSet`]. Consider [`insert_kmer_one_mismatch_checked`] when it
    /// is not known whether the bases and k-mer length will be valid.
    ///
    /// [`insert_kmer_one_mismatch_checked`]:
    ///     KmerSet::insert_kmer_one_mismatch_checked
    #[inline]
    pub fn insert_kmer_one_mismatch<Q: AsRef<[u8]>>(&mut self, kmer: Q) {
        self.insert_encoded_kmer_one_mismatch(self.encoder.encode_kmer(kmer));
    }

    /// Insert all k-mers into the set with at most one mismatch compared to the
    /// provided k-mer. The original k-mer is also inserted. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`], then `false` is returned and no insertion is performed.
    #[inline]
    pub fn insert_kmer_one_mismatch_checked<Q: AsRef<[u8]>>(&mut self, kmer: Q) -> bool {
        let Some(encoded_kmer) = self.encoder.encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer_one_mismatch(encoded_kmer);
        true
    }

    /// Insert all k-mers from a sequence into the [`KmerSet`]. The bases in the
    /// sequence must be valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`].
    #[inline]
    pub fn insert_from_sequence<Q: AsRef<[u8]>>(&mut self, seq: Q) {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer(encoded_kmer);
        }
    }

    /// Insert all k-mers from a sequence into the [`KmerSet`], in addition to
    /// all k-mers with up to one mismatch from those in the sequence. The bases
    /// in the sequence must be valid for the [`KmerEncoder`] associated with
    /// this [`KmerSet`].
    #[inline]
    pub fn insert_from_sequence_one_mismatch<Q: AsRef<[u8]>>(&mut self, seq: Q) {
        for encoded_kmer in self.encoder.iter_from_sequence(&seq) {
            self.insert_encoded_kmer_one_mismatch(encoded_kmer);
        }
    }

    /// Return the indices of the leftmost occurrence of any of the k-mers in
    /// this set within a provided sequence. The bases in the sequence must be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. If no
    /// occurrence is found, then `None` is returned.
    pub fn find_kmers<Q: AsRef<[u8]>>(&self, seq: Q) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder.iter_from_sequence(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                return Some(i..i + self.get_kmer_length());
            }
        }
        None
    }

    /// Return the indices of the rightmost occurrence of any of the k-mers in
    /// this set within a provided sequence. The bases in the sequence must be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. If no
    /// occurrence is found, then `None` is returned.
    pub fn find_kmers_rev<Q: AsRef<[u8]>>(&self, seq: Q) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder.iter_from_sequence_rev(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                let end = seq.as_ref().len() - i;
                return Some(end - self.get_kmer_length()..end);
            }
        }
        None
    }
}

impl<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>, S> IntoIterator for KmerSet<MAX_LEN, E, S>
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
