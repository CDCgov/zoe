use crate::kmer::{
    KmerEncoder, SupportedKmerLen,
    encoders::three_bit::{ThreeBitKmerEncoder, ThreeBitKmerLen, ThreeBitMismatchIter, ThreeBitOneMismatchIter},
};

/// A trait for specializing the implementation of [`KmerEncoder::get_variants`]
/// based on the number of mismatches.
pub trait GetVariants<const N: usize, const MAX_LEN: usize>: KmerEncoder<MAX_LEN> {
    /// The iterator for generating the specified variants.
    type Iter: Iterator<Item = Self::EncodedKmer>;

    /// Forms an iterator over the variants of `encoded_kmer`. Consider using
    /// [`KmerEncoder::get_variants`] instead.
    #[must_use]
    fn variants(&self, encoded_kmer: Self::EncodedKmer) -> Self::Iter;
}

impl<const MAX_LEN: usize> GetVariants<1, MAX_LEN> for ThreeBitKmerEncoder<MAX_LEN>
where
    ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
{
    type Iter = ThreeBitOneMismatchIter<MAX_LEN>;

    #[inline]
    fn variants(&self, encoded_kmer: Self::EncodedKmer) -> Self::Iter {
        ThreeBitOneMismatchIter::new(encoded_kmer, self)
    }
}

/// A macro for implementing [`GetVariants`] for [`ThreeBitKmerEncoder`] for
/// values of `N` greater than 1.
macro_rules! impl_three_bit_get_variants {
    ($($n:literal),* $(,)?) => {
        $(
            impl<const MAX_LEN: usize> GetVariants<$n, MAX_LEN>
                for ThreeBitKmerEncoder<MAX_LEN>
            where
                ThreeBitKmerLen<MAX_LEN>: SupportedKmerLen,
            {
                type Iter = ThreeBitMismatchIter<MAX_LEN, $n>;

                fn variants(&self, encoded_kmer: Self::EncodedKmer) -> Self::Iter {
                    ThreeBitMismatchIter::new(encoded_kmer, self)
                }
            }
        )*
    };
}

impl_three_bit_get_variants!(2, 3, 4, 5, 6, 7, 8, 9, 10);
