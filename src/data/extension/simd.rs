use std::{
    ops::{AddAssign, Shl, Shr, SubAssign},
    simd::{LaneCount, SimdElement, SupportedLaneCount, prelude::*},
};

// Do not use in a prelude to avoid conflicts with `std::simd::SimdInt` and `std::simd::SimdUint`.
pub trait SimdAnyInt<T, const N: usize>:
    SimdOrd
    + SimdPartialEq<Mask = Mask<<T as SimdElement>::Mask, N>>
    + AddAssign<Simd<T, N>>
    + SubAssign<Simd<T, N>>
    + Shl<T, Output = Simd<T, N>>
    + Shr<T, Output = Simd<T, N>>
where
    LaneCount<N>: SupportedLaneCount,
    T: SimdElement, {
    #[must_use]
    fn reduce_max(self) -> T;
    #[must_use]
    fn reduce_min(self) -> T;
    #[must_use]
    fn saturating_sub(self, rhs: Self) -> Self;
    #[must_use]
    fn saturating_add(self, rhs: Self) -> Self;
}

macro_rules! impl_simd_any_int_signed {
    ($($t:ty),*) => {
        $(
            impl<const N: usize> SimdAnyInt<$t, N> for Simd<$t, N>
            where
                LaneCount<N>: SupportedLaneCount,
            {
                #[inline]
                fn reduce_max(self) -> $t {
                    <Self as std::simd::num::SimdInt>::reduce_max(self)
                }

                #[inline]
                fn reduce_min(self) -> $t {
                    <Self as std::simd::num::SimdInt>::reduce_min(self)
                }

                #[inline]
                fn saturating_sub(self, rhs: Self) -> Self {
                    <Self as std::simd::num::SimdInt>::saturating_sub(self, rhs)
                }

                #[inline]
                fn saturating_add(self, rhs: Self) -> Self {
                    <Self as std::simd::num::SimdInt>::saturating_add(self, rhs)
                }
            }
        )*
    };
}
impl_simd_any_int_signed!(i8, i16, i32, i64, isize);

macro_rules! impl_simd_any_int_unsigned {
    ($($t:ty),*) => {
        $(
            impl<const N: usize> SimdAnyInt<$t, N> for Simd<$t, N>
            where
                LaneCount<N>: SupportedLaneCount,
            {
                #[inline]
                fn reduce_max(self) -> $t {
                    <Self as std::simd::num::SimdUint>::reduce_max(self)
                }

                #[inline]
                fn reduce_min(self) -> $t {
                    <Self as std::simd::num::SimdUint>::reduce_min(self)
                }

                #[inline]
                fn saturating_sub(self, rhs: Self) -> Self {
                    <Self as std::simd::num::SimdUint>::saturating_sub(self, rhs)
                }

                #[inline]
                fn saturating_add(self, rhs: Self) -> Self {
                    <Self as std::simd::num::SimdUint>::saturating_add(self, rhs)
                }
            }
        )*
    };
}
impl_simd_any_int_unsigned!(u8, u16, u32, u64, usize);

#[allow(dead_code)]
pub(crate) trait SimdByteFunctions<const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    fn is_ascii(&self) -> Mask<i8, N>;

    fn is_ascii_uppercase(&self) -> Mask<i8, N>;
    fn is_ascii_lowercase(&self) -> Mask<i8, N>;
    fn is_ascii_alphabetic(&self) -> Mask<i8, N>;
    fn is_ascii_digit(&self) -> Mask<i8, N>;
    fn is_ascii_whitespace(&self) -> Mask<i8, N>;
    fn is_ascii_graphic(&self) -> Mask<i8, N>;

    fn to_ascii_uppercase(&self) -> Self;
    fn to_ascii_lowercase(&self) -> Self;

    fn make_ascii_uppercase(&mut self);
    fn make_ascii_lowercase(&mut self);

    /// All instances of `find` are mapped to `replace`.
    ///
    /// ## Examples:
    /// ```ignore
    ///     let mut s = Simd::from_array([1,2,4,1]);
    ///     s.if_value_then_replace(1, 4);
    ///
    ///     assert_eq!(s, Simd::from_array([4,2,4,4]) );
    /// ```
    fn if_value_then_replace(&mut self, find: u8, replace: u8);

    /// All instances of `this` or `that` byte are swapped or exchanged.
    ///
    /// ## Examples:
    /// ```ignore
    ///     let mut s = Simd::from_array([1,2,4,1]);
    ///     s.exchange_byte_pairs(1, 4);
    ///
    ///     assert_eq!(s, Simd::from_array([4,2,1,4]) );
    /// ```
    ///
    fn exchange_byte_pairs(&mut self, this: u8, that: u8);
}

impl<const N: usize> SimdByteFunctions<N> for Simd<u8, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn is_ascii(&self) -> Mask<i8, N> {
        self.simd_lt(Simd::splat(128))
    }

    #[inline]
    fn is_ascii_uppercase(&self) -> Mask<i8, N> {
        self.simd_ge(Simd::splat(b'A')) & self.simd_le(Simd::splat(b'Z'))
    }

    #[inline]
    fn is_ascii_lowercase(&self) -> Mask<i8, N> {
        self.simd_ge(Simd::splat(b'a')) & self.simd_le(Simd::splat(b'z'))
    }

    #[inline]
    fn is_ascii_alphabetic(&self) -> Mask<i8, N> {
        self.is_ascii_lowercase() & self.is_ascii_uppercase()
    }

    #[inline]
    fn is_ascii_digit(&self) -> Mask<i8, N> {
        self.simd_ge(Simd::splat(b'0')) & self.simd_le(Simd::splat(b'9'))
    }

    #[inline]
    /// Checks for graphic ASCII characters from `!` (33) to `~` (126).
    fn is_ascii_graphic(&self) -> Mask<i8, N> {
        self.simd_ge(Simd::splat(b'!')) & self.simd_le(Simd::splat(b'~'))
    }

    #[inline]
    /// Tests is bytes are tab, line feed,  vertical tab, form feed, carriage
    /// return or space. This is different than in `std`.
    fn is_ascii_whitespace(&self) -> Mask<i8, N> {
        (self.simd_ge(Simd::splat(b'\t')) & self.simd_le(Simd::splat(b'\r'))) | self.simd_eq(Simd::splat(b' '))
    }

    #[inline]
    fn to_ascii_uppercase(&self) -> Self {
        let mask = self.is_ascii_lowercase();
        mask.select(*self ^ Simd::splat(0b0010_0000), *self)
    }

    #[inline]
    fn to_ascii_lowercase(&self) -> Self {
        let mask = self.is_ascii_uppercase();
        mask.select(*self | Simd::splat(0b0010_0000), *self)
    }

    #[inline]
    fn make_ascii_uppercase(&mut self) {
        *self = self.to_ascii_uppercase();
    }

    #[inline]
    fn make_ascii_lowercase(&mut self) {
        *self = self.to_ascii_lowercase();
    }

    #[inline]
    fn if_value_then_replace(&mut self, find: u8, replace: u8) {
        let mask = self.simd_eq(Simd::splat(find));
        *self = mask.select(Simd::splat(replace), *self);
    }

    #[inline]
    fn exchange_byte_pairs(&mut self, this: u8, that: u8) {
        let splat_this = Simd::splat(this);
        let splat_that = Simd::splat(that);

        let mask_this = self.simd_eq(splat_this);
        let mask_that = self.simd_eq(splat_that);

        let halfway = mask_this.select(splat_that, *self);
        *self = mask_that.select(splat_this, halfway);
    }
}

pub(crate) trait SimdMaskFunctions<const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    #[must_use]
    fn make_selected_ascii_uppercase(&self, bytes: &Simd<u8, N>) -> Simd<u8, N>;
    #[must_use]
    fn make_selected_ascii_lowercase(&self, bytes: &Simd<u8, N>) -> Simd<u8, N>;
    #[must_use]
    fn bitmask_offset(&self) -> usize;
}

impl<const N: usize> SimdMaskFunctions<N> for Mask<i8, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn make_selected_ascii_uppercase(&self, bytes: &Simd<u8, N>) -> Simd<u8, N> {
        self.select(*bytes ^ Simd::splat(0b0010_0000), *bytes)
    }

    #[inline]
    fn make_selected_ascii_lowercase(&self, bytes: &Simd<u8, N>) -> Simd<u8, N> {
        self.select(*bytes | Simd::splat(0b0010_0000), *bytes)
    }

    #[inline]
    fn bitmask_offset(&self) -> usize {
        self.to_bitmask().trailing_zeros() as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simd_upper_lowercase() {
        let test1 = b"aeiouA!?.@XYZqwt".to_owned();
        let simd = Simd::from_array(test1);

        assert_eq!(
            String::from_utf8_lossy(&test1.to_ascii_uppercase()),
            String::from_utf8_lossy(simd.to_ascii_uppercase().as_array())
        );
    }
}
