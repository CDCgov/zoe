use crate::kmer::{encoder::KmerEncoder, kmer_set::KmerSet};
use std::ops::Index;

/// [`KmerCounter`] is a container for storing counts of encoded kmers. A
/// [`KmerCounter`] can also be viewed as a multiset, hence why every
/// [`KmerCounter`] also implements [`KmerSet`] for inserting and checking the
/// existence of kmers. For example, [`insert_kmer`] has the impact of
/// incrementing the count of a kmer (or setting it to 1 if it was not
/// previously present).
///
/// [`insert_kmer`]: KmerSet::insert_kmer
pub trait KmerCounter: KmerSet + Index<Self::Kmer> {
    /// Get the count of an encoded kmer. If the kmer is not present in the
    /// counter, then `0` is returned. The encoded kmer must have been generated
    /// using the encoder associated with this [`KmerCounter`].
    fn get_encoded(&self, kmer: Self::Kmer) -> usize;

    /// Get the count of a kmer. If the kmer is not present in the counter, then
    /// `0` is returned. The bases and kmer length are assumed to be valid for
    /// the encoder associated with this [`KmerCounter`]. Consider
    /// [`get_checked`] when it is not known whether the bases and kmer length
    /// will be valid.
    ///
    /// [`get_checked`]: KmerCounter::get_checked
    #[inline]
    fn get(&self, kmer: &[u8]) -> usize {
        self.get_encoded(self.get_encoder().encode_kmer(kmer))
    }

    /// Get the count of a kmer. If the kmer is not present in the counter, then
    /// `0` is returned. If the bases and kmer length are not valid for the
    /// encoder associated with this [`KmerCounter`], then `None` is returned.
    #[inline]
    fn get_checked(&self, kmer: &[u8]) -> Option<usize> {
        Some(self.get_encoded(self.get_encoder().encode_kmer_checked(kmer)?))
    }
}
