use crate::kmer::encoder::KmerEncoder;
use std::{hash::Hash, ops::Range};

/// [`KmerSet`] represents a set of encoded k-mers along with relevant methods
/// that can be applied using it. For instance, a [`KmerSet`] can be used to
/// find the leftmost or rightmost occurrence of k-mers in a sequence, using
/// [`find_kmers`] and [`find_kmers_rev`]. A [`KmerSet`] can also be constructed
/// using [`insert_from_sequence`].
///
/// [`find_kmers`]: KmerSet::find_kmers
/// [`find_kmers_rev`]: KmerSet::find_kmers_rev
/// [`insert_from_sequence`]: KmerSet::insert_from_sequence
pub trait KmerSet<const MAX_LEN: usize>: IntoIterator {
    /// The type of an encoded k-mer.
    type EncodedKmer: Hash;
    /// The type of the encoder.
    type Encoder: KmerEncoder<MAX_LEN, EncodedKmer = Self::EncodedKmer>;
    /// The type of the non-consuming iterator over decoded k-mers.
    type DecodedIter<'a>
    where
        Self: 'a;
    /// The type of the non-consuming iterator over encoded k-mers.
    type EncodedIter<'a>
    where
        Self: 'a;

    /// Get the encoder used for the [`KmerSet`].
    fn get_encoder(&self) -> &Self::Encoder;

    /// Insert an encoded k-mer into the set. The encoded k-mer must have been
    /// generated using the encoder associated with this [`KmerSet`].
    fn insert_encoded_kmer(&mut self, encoded_kmer: Self::EncodedKmer);

    /// Return whether an encoded k-mer is present in the set. The encoded k-mer
    /// must have been generated using the encoder associated with this
    /// [`KmerSet`].
    fn contains_encoded(&self, kmer: Self::EncodedKmer) -> bool;

    /// Iterate over the decoded k-mers in the set.
    fn iter_decoded(&self) -> Self::DecodedIter<'_>;

    /// Iterate over the encoded k-mers in the set.
    fn iter_encoded(&self) -> Self::EncodedIter<'_>;

    /// Get the length of the k-mers being stored in the set.
    #[inline]
    fn get_kmer_length(&self) -> usize {
        self.get_encoder().get_kmer_length()
    }

    /// Insert a k-mer into the set. The bases and k-mer length are assumed to be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. Consider
    /// [`insert_kmer_checked`] when it is not known whether the bases and k-mer
    /// length will be valid.
    ///
    /// [`insert_kmer_checked`]: KmerSet::insert_kmer_checked
    #[inline]
    fn insert_kmer<S: AsRef<[u8]>>(&mut self, kmer: S) {
        self.insert_encoded_kmer(self.get_encoder().encode_kmer(kmer));
    }

    /// Insert a k-mer into the set. If the bases and k-mer length are not valid
    /// for the [`KmerEncoder`] associated with this [`KmerSet`], then `false`
    /// is returned and no insertion is performed.
    #[inline]
    fn insert_kmer_checked<S: AsRef<[u8]>>(&mut self, kmer: S) -> bool {
        let Some(encoded_kmer) = self.get_encoder().encode_kmer_checked(kmer) else {
            return false;
        };
        self.insert_encoded_kmer(encoded_kmer);
        true
    }

    /// Insert all k-mers from a sequence into the [`KmerSet`]. The bases in the
    /// sequence must be valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`].
    #[inline]
    fn insert_from_sequence<S: AsRef<[u8]>>(&mut self, seq: S) {
        for encoded_kmer in self.get_encoder().iter_from_sequence(&seq) {
            self.insert_encoded_kmer(encoded_kmer);
        }
    }

    /// Return whether a k-mer is present in the set. The bases and k-mer length
    /// are assumed to be valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`]. Consider [`contains_checked`] when it is not known whether
    /// the bases and k-mer length will be valid.
    ///
    /// [`contains_checked`]: KmerSet::contains_checked
    #[inline]
    fn contains<S: AsRef<[u8]>>(&self, kmer: S) -> bool {
        self.contains_encoded(self.get_encoder().encode_kmer(kmer))
    }

    /// Return whether a k-mer is present in the set. If the bases and k-mer
    /// length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`], then `None` is returned.
    #[inline]
    fn contains_checked<S: AsRef<[u8]>>(&self, kmer: S) -> Option<bool> {
        Some(self.contains_encoded(self.get_encoder().encode_kmer_checked(kmer)?))
    }

    /// Return the indices of the leftmost occurrence of any of the k-mers in
    /// this set within a provided sequence. The bases in the sequence must be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. If no
    /// occurrence is found, then `None` is returned.
    fn find_kmers<S: AsRef<[u8]>>(&self, seq: S) -> Option<Range<usize>> {
        for (i, kmer) in self.get_encoder().iter_from_sequence(&seq).enumerate() {
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
    fn find_kmers_rev<S: AsRef<[u8]>>(&self, seq: S) -> Option<Range<usize>> {
        for (i, kmer) in self.get_encoder().iter_from_sequence_rev(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                let end = seq.as_ref().len() - i;
                return Some(end - self.get_kmer_length() - i..end);
            }
        }
        None
    }
}
