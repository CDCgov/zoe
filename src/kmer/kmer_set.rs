use crate::kmer::encoder::KmerEncoder;
use std::hash::Hash;

/// [`KmerSet`] represents a set of encoded k-mers along with relevant methods
/// that can be applied using it. For instance, a [`KmerSet`] can be used to
/// find the leftmost or rightmost occurrence of k-mers in a sequence, using
/// [`find_kmers`] and [`find_kmers_rev`]. A [`KmerSet`] can also be constructed
/// in efficient ways, such as using [`insert_from_sequence`].
///
/// [`find_kmers`]: KmerSet::find_kmers
/// [`find_kmers_rev`]: KmerSet::find_kmers_rev
/// [`insert_from_sequence`]: KmerSet::insert_from_sequence
pub trait KmerSet: IntoIterator {
    /// The type of an encoded kmer.
    type EncodedKmer: Hash;
    /// The type of the encoder.
    type Encoder: KmerEncoder<EncodedKmer = Self::EncodedKmer>;

    /// Get the encoder used for the [`KmerSet`].
    fn get_encoder(&self) -> &Self::Encoder;

    /// Insert an encoded k-mer into the set. The encoded k-mer must have been
    /// generated using the encoder associated with this [`KmerSet`].
    fn insert_encoded_kmer(&mut self, encoded_kmer: Self::EncodedKmer);

    /// Return whether an encoded k-mer is present in the set. The encoded k-mer
    /// must have been generated using the encoder associated with this
    /// [`KmerSet`].
    fn contains_encoded(&self, kmer: Self::EncodedKmer) -> bool;

    /// Insert a k-mer into the set. The bases and k-mer length are assumed to be
    /// valid for the [`KmerEncoder`] associated with this [`KmerSet`]. Consider
    /// [`insert_kmer_checked`] when it is not known whether the bases and k-mer
    /// length will be valid.
    ///
    /// [`insert_kmer_checked`]: KmerSet::insert_kmer_checked
    #[inline]
    fn insert_kmer(&mut self, kmer: &[u8]) {
        self.insert_encoded_kmer(self.get_encoder().encode_kmer(kmer));
    }

    /// Insert a k-mer into the set. If the bases and k-mer length are not valid
    /// for the [`KmerEncoder`] associated with this [`KmerSet`], then `false`
    /// is returned and no insertion is performed.
    #[inline]
    fn insert_kmer_checked(&mut self, kmer: &[u8]) -> bool {
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
        for encoded_kmer in self.get_encoder().iter_from_sequence(seq.as_ref()) {
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
    fn contains(&self, kmer: &[u8]) -> bool {
        self.contains_encoded(self.get_encoder().encode_kmer(kmer))
    }

    /// Return whether a k-mer is present in the set. If the bases and k-mer
    /// length are not valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`], then `None` is returned.
    #[inline]
    fn contains_checked(&self, kmer: &[u8]) -> Option<bool> {
        Some(self.contains_encoded(self.get_encoder().encode_kmer_checked(kmer)?))
    }

    /// Return the starting index of the leftmost occurrence of any of the k-mers
    /// in this set within a provided sequence. The bases in the
    /// sequence must be valid for the [`KmerEncoder`] associated with this
    /// [`KmerSet`]. If no occurrence is found, then `None` is returned.
    fn find_kmers<S: AsRef<[u8]>>(&self, seq: S) -> Option<usize> {
        for (i, kmer) in self.get_encoder().iter_from_sequence(seq.as_ref()).enumerate() {
            if self.contains_encoded(kmer) {
                return Some(i);
            }
        }
        None
    }

    /// Return the starting index of the rightmost occurrence of any of the
    /// k-mers in this set within a provided sequence. The bases in the sequence
    /// must be valid for the [`KmerEncoder`] associated with this [`KmerSet`].
    /// If no occurrence is found, then `None` is returned.
    fn find_kmers_rev<S: AsRef<[u8]>>(&self, seq: S) -> Option<usize> {
        for (i, kmer) in self.get_encoder().iter_from_sequence_rev(seq.as_ref()).enumerate() {
            if self.contains_encoded(kmer) {
                return Some(seq.as_ref().len() - self.get_encoder().get_kmer_length() - i);
            }
        }
        None
    }
}
