use crate::kmer::{encoder::KmerEncoder, kmer_set::KmerSet};
use std::ops::Index;

/// [`KmerCounter`] is a container for storing counts of encoded k-mers. A
/// [`KmerCounter`] can also be viewed as a multiset, hence why every
/// [`KmerCounter`] also implements [`KmerSet`] for inserting and checking the
/// existence of k-mers. For example, [`insert_kmer`] has the impact of
/// incrementing the count of a k-mer (or setting it to 1 if it was not
/// previously present).
///
/// [`insert_kmer`]: KmerSet::insert_kmer
pub trait KmerCounter<const MAX_LEN: usize>: KmerSet<MAX_LEN> + Index<Self::EncodedKmer> {
    /// Get the count of an encoded k-mer. If the k-mer is not present in the
    /// counter, then `0` is returned. The encoded k-mer must have been
    /// generated using the encoder associated with this [`KmerCounter`].
    fn get_encoded(&self, kmer: Self::EncodedKmer) -> u64;

    /// Get the count of a k-mer. If the k-mer is not present in the counter,
    /// then `0` is returned. The bases and k-mer length are assumed to be valid
    /// for the encoder associated with this [`KmerCounter`]. Consider
    /// [`get_checked`] when it is not known whether the bases and k-mer length
    /// will be valid.
    ///
    /// [`get_checked`]: KmerCounter::get_checked
    #[inline]
    fn get<S: AsRef<[u8]>>(&self, kmer: S) -> u64 {
        self.get_encoded(self.get_encoder().encode_kmer(kmer))
    }

    /// Get the count of a k-mer. If the k-mer is not present in the counter, then
    /// `0` is returned. If the bases and k-mer length are not valid for the
    /// encoder associated with this [`KmerCounter`], then `None` is returned.
    #[inline]
    fn get_checked<S: AsRef<[u8]>>(&self, kmer: S) -> Option<u64> {
        Some(self.get_encoded(self.get_encoder().encode_kmer_checked(kmer)?))
    }
}
