use std::{
    hash::Hash,
    ops::{Add, AddAssign, BitAnd, BitOr, Not, Shl, ShlAssign, Shr, ShrAssign, Sub},
    simd::SimdElement,
};

pub trait Int:
    SimdElement
    + Ord
    + Hash
    + Add<Output = Self>
    + Sub<Output = Self>
    + AddAssign
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + Not<Output = Self>
    + Shr<usize, Output = Self>
    + Shl<usize, Output = Self>
    + ShrAssign<usize>
    + ShlAssign<usize> {
    fn zero() -> Self;
    fn one() -> Self;

    fn bit_0b100() -> Self;
    fn bit_0b111() -> Self;

    fn as_usize(self) -> usize;
    fn checked_addition(&self, other: Self) -> Option<Self>;
    const MAX: Self;
    const MIN: Self;
    const BITS: u32;
}

pub trait Uint: Int + From<u8> {}

macro_rules! impl_int {
    { $($ty:ty),* } => {
        $(
        impl Int for $ty {
            #[inline]
            fn zero() -> Self { 0 }
            #[inline]
            fn one() -> Self { 1 }

            #[inline]
            fn bit_0b100() -> Self { 0b100 }
            #[inline]
            fn bit_0b111() -> Self { 0b111 }

            // TODO: FIX CLIPPY
            #[allow(renamed_and_removed_lints)]
            #[allow(clippy::cast_sign_loss, cast_possible_truncation)]
            #[inline]
            fn as_usize(self) -> usize { self as usize }

            #[inline]
            fn checked_addition(&self,other: Self) -> Option<Self> {
                self.checked_add(other)
            }
            const MAX: $ty = <$ty>::MAX;
            const MIN: $ty = <$ty>::MIN;
            const BITS: u32 = <$ty>::BITS;
        } )*

     }
}

macro_rules! impl_uint {
    { $($ty:ty),* } => {
        $( impl Uint for $ty {} )*
    }
}

impl_int!(i8, i16, i32, i64, isize, u8, u16, u32, u64, usize);
impl_uint!(u8, u16, u32, u64, usize);
