use crate::data::nucleotides::NucleotidesReadable;

use super::KmerEncoder;
use std::ops::Range;

pub trait EncodedKmerCollection<const MAX_LEN: usize> {
    type Encoder: KmerEncoder<MAX_LEN, EncodedKmer = Self::EncodedKmer>;
    type EncodedKmer;

    /// Gets the encoder associated with this collection.
    #[must_use]
    fn encoder(&self) -> &Self::Encoder;

    /// Gets the length of the k-mers being stored in this collection.
    #[inline]
    #[must_use]
    fn kmer_length(&self) -> usize {
        self.encoder().kmer_length()
    }
}

/// Trait for k-mer collections with a `contains` method (i.e., that can act as
/// a set).
pub trait KmerCollectionContains<const MAX_LEN: usize>: EncodedKmerCollection<MAX_LEN> {
    /// Return whether an already encoded k-mer is present in the collection.
    /// The encoded k-mer must have been generated using the [`KmerEncoder`]
    /// associated with this collection.
    #[must_use]
    fn contains_encoded(&self, kmer: Self::EncodedKmer) -> bool;

    /// Return whether a k-mer is present in the collection. The bases and k-mer
    /// length are assumed to be valid for the [`KmerEncoder`] associated with
    /// this collection. Consider [`contains_checked`] when it is not known
    /// whether the bases and k-mer length will be valid.
    ///
    /// [`contains_checked`]: KmerCollectionContains::contains_checked
    #[inline]
    #[must_use]
    fn contains(&self, kmer: impl AsRef<[u8]>) -> bool {
        self.contains_encoded(self.encoder().encode_kmer(kmer))
    }

    /// Return whether a k-mer is present in the collection. If the bases and
    /// k-mer length are not valid for the [`KmerEncoder`] associated with this
    /// collection, then `None` is returned.
    #[inline]
    #[must_use]
    fn contains_checked(&self, kmer: impl AsRef<[u8]>) -> Option<bool> {
        Some(self.contains_encoded(self.encoder().encode_kmer_checked(kmer)?))
    }

    /// Return the indices of the leftmost occurrence of any of the k-mers in
    /// this collection within a provided sequence. The bases in the sequence
    /// must be valid for the [`KmerEncoder`] associated with this collection.
    /// If no occurrence is found, then `None` is returned.
    #[must_use]
    fn find_in_seq(&self, seq: impl AsRef<[u8]>) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder().iter_from_sequence(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                return Some(i..i + self.kmer_length());
            }
        }
        None
    }

    /// Return the indices of the rightmost occurrence of any of the k-mers in
    /// this collection within a provided sequence. The bases in the sequence
    /// must be valid for the [`KmerEncoder`] associated with this collection.
    /// If no occurrence is found, then `None` is returned.
    #[must_use]
    fn find_in_seq_rev(&self, seq: impl AsRef<[u8]>) -> Option<Range<usize>> {
        for (i, kmer) in self.encoder().iter_from_sequence_rev(&seq).enumerate() {
            if self.contains_encoded(kmer) {
                let end = seq.as_ref().len() - i;
                return Some(end - self.kmer_length()..end);
            }
        }
        None
    }
}

/// Provides methods for searching for a kmer in a collection.
pub trait FindKmers<const MAX_LEN: usize>: NucleotidesReadable {
    /// Return the indices of the leftmost occurrence of any of the k-mers in a
    /// collection. The bases in the sequence must be valid for the
    /// [`KmerEncoder`] associated with the collection. If no occurrence is
    /// found, then `None` is returned.
    #[inline]
    #[must_use]
    fn find_kmers<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        kmers.find_in_seq(self.nucleotide_bytes())
    }

    /// Return the indices of the rightmost occurrence of any of the k-mers in a
    /// collection. The bases in the sequence must be valid for the
    /// [`KmerEncoder`] associated with the collection. If no occurrence is
    /// found, then `None` is returned.
    #[inline]
    #[must_use]
    fn find_kmers_rev<T: KmerCollectionContains<MAX_LEN>>(&self, kmers: &T) -> Option<Range<usize>> {
        kmers.find_in_seq_rev(self.nucleotide_bytes())
    }
}

impl<const MAX_LEN: usize, T: NucleotidesReadable> FindKmers<MAX_LEN> for T {}
