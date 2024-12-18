use std::{
    hash::Hash,
    ops::{Add, AddAssign, BitAnd, BitOr, Not, Shl, ShlAssign, Shr, ShrAssign, Sub},
};

pub trait Int:
    Sized
    + Copy
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
            const ZERO: $ty = 0;
            const ONE: $ty = 1;

            #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
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

impl_int!(i8, i16, i32, i64, isize, i128, u8, u16, u32, u64, usize, u128);
impl_uint!(u8, u16, u32, u64, usize, u128);
