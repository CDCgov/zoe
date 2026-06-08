//! Helper traits/structs for supporting compile time dispatch and type mappings

use crate::{
    kmer::{
        KmerEncoder,
        encoders::{three_bit::ThreeBitKmerEncoder, two_bit::TwoBitKmerEncoder},
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
/// Supported max k-mer lengths for [`ThreeBitKmerEncoder`] are:
///
/// - `2`: uses a `u8`
/// - `3..=5`: uses a `u16`
/// - `6..=10`: uses a `u32`
/// - `11..=21`: uses a `u64`
/// - `22..=42`: uses a `u128`
///
/// For a [`TwoBitKmerEncoder`], they are:
///
/// - `2..=4`: uses a `u8`
/// - `5..=8`: uses a `u16`
/// - `9..=16`: uses a `u32`
/// - `17..=32`: uses a `u64`
/// - `33..=64`: uses a `u128`
///
/// This must never be implemented on any [`KmerLen`] where `MAX_LEN > 255`,
/// since [`Kmer`] uses a `u8` to store the k-mer length.
///
/// [`ThreeBitKmerEncoder`]:
///     crate::kmer::encoders::three_bit::ThreeBitKmerEncoder
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

impl_supported_kmer_len! {
    TwoBitKmerEncoder,
    2 => u8,
    3 => u8,
    4 => u8,
    5 => u16,
    6 => u16,
    7 => u16,
    8 => u16,
    9 => u32,
    10 => u32,
    11 => u32,
    12 => u32,
    13 => u32,
    14 => u32,
    15 => u32,
    16 => u32,
    17 => u64,
    18 => u64,
    19 => u64,
    20 => u64,
    21 => u64,
    22 => u64,
    23 => u64,
    24 => u64,
    25 => u64,
    26 => u64,
    27 => u64,
    28 => u64,
    29 => u64,
    30 => u64,
    31 => u64,
    32 => u64,
    33 => u128,
    34 => u128,
    35 => u128,
    36 => u128,
    37 => u128,
    38 => u128,
    39 => u128,
    40 => u128,
    41 => u128,
    42 => u128,
    43 => u128,
    44 => u128,
    45 => u128,
    46 => u128,
    47 => u128,
    48 => u128,
    49 => u128,
    50 => u128,
    51 => u128,
    52 => u128,
    53 => u128,
    54 => u128,
    55 => u128,
    56 => u128,
    57 => u128,
    58 => u128,
    59 => u128,
    60 => u128,
    61 => u128,
    62 => u128,
    63 => u128,
    64 => u128
}

/// A type alias for easily getting the encoded k-mer integer type from
/// `MAX_LEN` and the [`KmerEncoder`].
pub(crate) type MaxLenToType<const MAX_LEN: usize, E> = <KmerLen<MAX_LEN, E> as SupportedKmerLen>::T;
