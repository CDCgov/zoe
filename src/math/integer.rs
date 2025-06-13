use std::{
    fmt::Debug,
    hash::Hash,
    ops::{Add, AddAssign, BitAnd, BitOr, Not, Shl, ShlAssign, Shr, ShrAssign, Sub, SubAssign},
};

use crate::{
    math::cast::{CastAsNumeric, CastFromBool, CastFromNumeric},
    private::Sealed,
};

pub trait AnyInt:
    CastAsNumeric
    + CastFromNumeric
    + CastFromBool
    + Sized
    + Copy
    + Debug
    + Ord
    + Hash
    + Default
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + AddAssign
    + SubAssign
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + Not<Output = Self>
    + Shr<usize, Output = Self>
    + Shl<usize, Output = Self>
    + ShrAssign<usize>
    + ShlAssign<usize>
    + TryInto<usize>
    + Sealed {
    const ZERO: Self;
    const ONE: Self;

    fn from_literal<T: Int>(num: T) -> Self {
        let mut out = Self::ZERO;
        let mut i = T::ZERO;
        while i < num {
            out += Self::ONE;
            i += T::ONE;
        }
        out
    }

    fn checked_add(self, other: Self) -> Option<Self>;
    fn trailing_zeros(self) -> u32;
    fn count_ones(self) -> u32;

    #[must_use]
    fn wrapping_shl(self, rhs: u32) -> Self;

    #[must_use]
    fn rotate_left(self, n: u32) -> Self;

    #[must_use]
    fn rotate_right(self, n: u32) -> Self;

    const MAX: Self;
    const MIN: Self;
    const BITS: u32;
    const SIGNED: bool;
}

pub trait Uint: AnyInt + From<u8> {}
pub trait Int: AnyInt + From<i8> {}

macro_rules! impl_int {
    { $($ty:ty),* } => {
        $(
        impl AnyInt for $ty {
            const ZERO: $ty = 0;
            const ONE: $ty = 1;
            const SIGNED: bool = stringify!($ty).as_bytes()[0] == b'i';

            #[inline]
            fn checked_add(self, other: Self) -> Option<Self> {
                self.checked_add(other)
            }

            #[inline]
            fn trailing_zeros(self) -> u32 {
                self.trailing_zeros()
            }

            #[inline]
            fn count_ones(self) -> u32 {
                self.count_ones()
            }

            #[inline]
            fn wrapping_shl(self, rhs: u32) -> Self {
                self.wrapping_shl(rhs)
            }

            #[inline]
            fn rotate_left(self, n: u32) -> Self {
                self.rotate_left(n)
            }

            #[inline]
            fn rotate_right(self, n: u32) -> Self {
                self.rotate_right(n)
            }

            const MAX: $ty = <$ty>::MAX;
            const MIN: $ty = <$ty>::MIN;
            const BITS: u32 = <$ty>::BITS;
        } )*

     }
}

macro_rules! impl_xint {
    ($tr:ty; $($ty:ty),*) => {
        $( impl $tr for $ty {} )*
    }
}

impl_int!(i8, i16, i32, i64, isize, i128, u8, u16, u32, u64, usize, u128);

impl_xint!(Uint; u8, u16, u32, u64, usize, u128);
impl_xint!(Int; i8, i16, i32, i64, isize, i128);

pub trait FromSameSignedness<T>: From<T> + Sealed {}

macro_rules! impl_from_same_sign {
    ($smaller:ty => $($bigger:ty),*) => {
        $(
            impl FromSameSignedness<$smaller> for $bigger {}
        )*
    };
}

impl_from_same_sign! {u8 => u8, u16, u32, u64, u128}
impl_from_same_sign! {u16 => u16, u32, u64, u128}
impl_from_same_sign! {u32 => u32, u64, u128}
impl_from_same_sign! {u64 => u64, u128}
impl_from_same_sign! {u128 => u128}

impl_from_same_sign! {i8 => i8, i16, i32, i64, i128}
impl_from_same_sign! {i16 => i16, i32, i64, i128}
impl_from_same_sign! {i32 => i32, i64, i128}
impl_from_same_sign! {i64 => i64, i128}
impl_from_same_sign! {i128 => i128}
