use crate::{
    kmer::{
        KmerEncoder,
        encoders::three_bit::{ThreeBitKmerEncoder, ThreeBitMismatchIter, ThreeBitMismatchNumber, ThreeBitOneMismatchIter},
    },
    math::Uint,
};
use std::marker::PhantomData;

/// A zero-size struct used to hold both a `MAX_LEN` argument and a
/// [`KmerEncoder`].
///
/// For valid values of `MAX_LEN`, [`SupportedKmerLen`] will be implemented, and
/// will provide the appropriate integer type to use for the encoded k-mer. See
/// [`SupportedKmerLen`] for more details.
pub struct KmerLen<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>>(PhantomData<E>);

/// Statically guarantees that a max k-mer length is marked as supported for the
/// given [`KmerEncoder`], as well as provides the appropriate integer type.
///
/// Supported max k-mer lengths for [`ThreeBitKmerEncoder`] include:
///
/// - `2`: uses a `u8`
/// - `3..=5`: uses a `u16`
/// - `6..=10`: uses a `u32`
/// - `11..=21`: uses a `u64`
/// - `22..=42`: uses a `u128`
///
/// This must never be implemented on any [`KmerLen`] where `MAX_LEN > 255`,
/// since [`Kmer`] uses a `u8` to store the kmer length.
///
/// [`ThreeBitKmerEncoder`]: crate::kmer::encoders::three_bit::ThreeBitKmerEncoder
/// [`Kmer`]: super::Kmer
pub trait SupportedKmerLen {
    /// The integer type used to store the encoded k-mer.
    type T: Uint;
}

/// A macro for implementing [`SupportedKmerLen`].
macro_rules! impl_supported_kmer_len {
    ($encoder:ident, $($max:expr => $type:ty),*) => {
        $(
            impl SupportedKmerLen for KmerLen<$max, $encoder<$max>> {
                type T = $type;
            }
        )*
    }
}

impl_supported_kmer_len! {
    ThreeBitKmerEncoder,
    2 => u8,
    3 => u16,
    4 => u16,
    5 => u16,
    6 => u32,
    7 => u32,
    8 => u32,
    9 => u32,
    10 => u32,
    11 => u64,
    12 => u64,
    13 => u64,
    14 => u64,
    15 => u64,
    16 => u64,
    17 => u64,
    18 => u64,
    19 => u64,
    20 => u64,
    21 => u64,
    22 => u128,
    23 => u128,
    24 => u128,
    25 => u128,
    26 => u128,
    27 => u128,
    28 => u128,
    29 => u128,
    30 => u128,
    31 => u128,
    32 => u128,
    33 => u128,
    34 => u128,
    35 => u128,
    36 => u128,
    37 => u128,
    38 => u128,
    39 => u128,
    40 => u128,
    41 => u128,
    42 => u128
}

/// A type alias for easily getting the encoded k-mer integer type from
/// `MAX_LEN` and the [`KmerEncoder`].
pub(crate) type MaxLenToType<const MAX_LEN: usize, E> = <KmerLen<MAX_LEN, E> as SupportedKmerLen>::T;

/// Statically guarantees that a number of mismatches is supported for
/// [`KmerEncoder::get_variants`], as well as provides the appropriate iterator
/// type.
///
/// Supported number of mismatches for [`ThreeBitKmerEncoder`] include 1-10. 1
/// uses a specialized implementation.
///
/// [`ThreeBitKmerEncoder`]:
///     crate::kmer::encoders::three_bit::ThreeBitKmerEncoder
pub trait SupportedMismatchNumber<const MAX_LEN: usize, E: KmerEncoder<MAX_LEN>> {
    /// The iterator type for the k-mer variants.
    type MismatchIter: Iterator<Item = E::EncodedKmer>;

    /// Returns the iterator over the k-mer variants, given the encoded k-mer
    /// and the encoder.
    fn get_iterator(encoded_kmer: E::EncodedKmer, encoder: &E) -> Self::MismatchIter;
}

macro_rules! impl_select_mismatch_iter {
    ($encoder:ident, $mismatch:ident, $($n:expr => $iter:ty),*) => {
        $(
            impl<const MAX_LEN: usize> SupportedMismatchNumber<MAX_LEN, $encoder<MAX_LEN>> for $mismatch<$n>
                where
                    KmerLen<MAX_LEN, $encoder<MAX_LEN>>: SupportedKmerLen,
                {
                    type MismatchIter = $iter;

                    #[inline]
                    fn get_iterator(
                        encoded_kmer: <$encoder<MAX_LEN> as KmerEncoder<MAX_LEN>>::EncodedKmer,
                        encoder: &$encoder<MAX_LEN>,
                    ) -> Self::MismatchIter {
                        Self::MismatchIter::new(encoded_kmer, &encoder)
                    }
                }
        )*
    }
}

impl_select_mismatch_iter! {
    ThreeBitKmerEncoder,
    ThreeBitMismatchNumber,
    1 => ThreeBitOneMismatchIter<MAX_LEN>,
    2 => ThreeBitMismatchIter<MAX_LEN, 2>,
    3 => ThreeBitMismatchIter<MAX_LEN, 3>,
    4 => ThreeBitMismatchIter<MAX_LEN, 4>,
    5 => ThreeBitMismatchIter<MAX_LEN, 5>,
    6 => ThreeBitMismatchIter<MAX_LEN, 6>,
    7 => ThreeBitMismatchIter<MAX_LEN, 7>,
    8 => ThreeBitMismatchIter<MAX_LEN, 8>,
    9 => ThreeBitMismatchIter<MAX_LEN, 9>,
    10 => ThreeBitMismatchIter<MAX_LEN, 10>
}
