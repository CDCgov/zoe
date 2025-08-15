use std::{
    fmt::Debug,
    hash::Hash,
    ops::{Add, AddAssign, BitAnd, BitOr, Not, Shl, ShlAssign, Shr, ShrAssign, Sub, SubAssign},
    simd::SimdElement,
};

use crate::{
    math::cast::{CastAsNumeric, CastFromBool, CastFromNumeric},
    private::Sealed,
};

/// A trait for providing generic functionality over integer numbers.
///
/// These can be signed or unsigned. For more specific traits, consider [`Int`]
/// for signed integers and [`Uint`] for unsigned integers.
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
    /// The constant zero.
    const ZERO: Self;

    /// The constant one.
    const ONE: Self;

    /// The largest value that can be represented by this integer type.
    const MAX: Self;

    /// The smallest value that can be represented by this integer type.
    const MIN: Self;

    /// The size of this integer type in bits.
    const BITS: u32;

    /// Whether or not this integer type is signed.
    const SIGNED: bool;

    /// Converts a literal number to this type via repeated addition.
    ///
    /// The compiler will often optimize this away.
    fn from_literal<T: Int>(num: T) -> Self {
        let mut out = Self::ZERO;
        let mut i = T::ZERO;
        while i < num {
            out += Self::ONE;
            i += T::ONE;
        }
        out
    }

    /// Checked integer addition. Computes `self + rhs`, returning `None` if
    /// overflow occurred.
    fn checked_add(self, other: Self) -> Option<Self>;

    /// Returns the number of trailing zeros in the binary representation of
    /// `self`.
    fn trailing_zeros(self) -> u32;

    /// Returns the number of ones in the binary representation of `self`.
    fn count_ones(self) -> u32;

    /// Panic-free bitwise shift-left; yields `self << mask(rhs)`, where `mask`
    /// removes any high-order bits of `rhs` that would cause the shift to
    /// exceed the bitwidth of the type.
    ///
    /// See [`u64::wrapping_shl`] for more details.
    #[must_use]
    fn wrapping_shl(self, rhs: u32) -> Self;

    /// Shifts the bits to the left by a specified amount, `n`, wrapping the
    /// truncated bits to the end of the resulting integer.
    ///
    /// [`u64::rotate_left`] for more details.
    #[must_use]
    fn rotate_left(self, n: u32) -> Self;

    /// Shifts the bits to the right by a specified amount, n, wrapping the
    /// truncated bits to the beginning of the resulting integer.
    ///
    /// [`u64::rotate_right`] for more details.
    #[must_use]
    fn rotate_right(self, n: u32) -> Self;
}

/// A trait for providing generic functionality over unsigned integer numbers.
///
/// See [`AnyInt`] for more details.
pub trait Uint: AnyInt + From<u8> {}

/// A trait for providing generic functionality over signed integer numbers.
///
/// See [`AnyInt`] for more details.
pub trait Int: AnyInt + From<i8> {}

/// A macro for implementing [`AnyInt`] on primitives.
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

/// A macro for implementing [`Uint`] and [`Int`] on primitives.
macro_rules! impl_xint {
    ($tr:ty; $($ty:ty),*) => {
        $( impl $tr for $ty {} )*
    }
}

impl_int!(i8, i16, i32, i64, isize, i128, u8, u16, u32, u64, usize, u128);

impl_xint!(Uint; u8, u16, u32, u64, usize, u128);
impl_xint!(Int; i8, i16, i32, i64, isize, i128);

/// A marker trait for integer types fo the same sign.
pub trait FromSameSignedness<T>: From<T> + Sealed {}

/// A macro for implementing [`FromSameSignedness`] on primitive types.
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

/// Valid integer widths for sequence alignment scores.
pub trait AlignableIntWidth: AnyInt + SimdElement + Sealed {}

impl AlignableIntWidth for u8 {}
impl AlignableIntWidth for u16 {}
impl AlignableIntWidth for u32 {}
impl AlignableIntWidth for i8 {}
impl AlignableIntWidth for i16 {}
impl AlignableIntWidth for i32 {}
