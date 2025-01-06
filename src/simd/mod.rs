use std::simd::{LaneCount, SimdElement, SupportedLaneCount, Swizzle, prelude::*};

#[allow(dead_code)]
pub trait SimdExt<T> {
    // based on rotate_elements_left and rotate_elements_right
    fn shift_elements_left_z<const OFFSET: usize>(self, padding: T) -> Self;
    fn shift_elements_right_z<const OFFSET: usize>(self, padding: T) -> Self;
}

impl<T, const N: usize> SimdExt<T> for Simd<T, N>
where
    T: SimdElement,
    LaneCount<N>: SupportedLaneCount,
{
    fn shift_elements_left_z<const OFFSET: usize>(self, padding: T) -> Self {
        struct Shift<const OFFSET: usize>;

        impl<const OFFSET: usize, const N: usize> Swizzle<N> for Shift<OFFSET> {
            const INDEX: [usize; N] = const {
                let mut index = [N; N];
                let mut i = 0;
                while i + OFFSET < N {
                    index[i] = i + OFFSET;
                    i += 1;
                }
                index
            };
        }

        Shift::<OFFSET>::concat_swizzle(self, Simd::splat(padding))
    }

    fn shift_elements_right_z<const OFFSET: usize>(self, padding: T) -> Self {
        struct Shift<const OFFSET: usize>;

        impl<const OFFSET: usize, const N: usize> Swizzle<N> for Shift<OFFSET> {
            const INDEX: [usize; N] = const {
                let mut index = [N; N];
                let mut i = OFFSET;
                while i < N {
                    index[i] = i - OFFSET;
                    i += 1;
                }
                index
            };
        }

        Shift::<OFFSET>::concat_swizzle(self, Simd::splat(padding))
    }
}

#[allow(dead_code)]
pub trait SimdByteFunctions<const N: usize>
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
    /// # Examples:
    /// ```ignore
    ///     let mut s = Simd::from_array([1,2,4,1]);
    ///     s.if_value_then_replace(1, 4);
    ///
    ///     assert_eq!(s, Simd::from_array([4,2,4,4]) );
    /// ```
    fn if_value_then_replace(&mut self, find: u8, replace: u8);

    /// All instances of `this` or `that` byte are swapped or exchanged.
    ///
    /// # Examples:
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
    fn make_selected_ascii_uppercase(&self, bytes: &Simd<u8, N>) -> Simd<u8, N>;
    fn make_selected_ascii_lowercase(&self, bytes: &Simd<u8, N>) -> Simd<u8, N>;
    fn bitmask_offset(&self) -> usize;
}

impl<const N: usize> SimdMaskFunctions<N> for Mask<i8, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[must_use]
    #[inline]
    fn make_selected_ascii_uppercase(&self, bytes: &Simd<u8, N>) -> Simd<u8, N> {
        self.select(*bytes ^ Simd::splat(0b0010_0000), *bytes)
    }

    #[must_use]
    #[inline]
    fn make_selected_ascii_lowercase(&self, bytes: &Simd<u8, N>) -> Simd<u8, N> {
        self.select(*bytes | Simd::splat(0b0010_0000), *bytes)
    }

    #[must_use]
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

    #[test]
    fn test_shift_elements_left() {
        let v = Simd::from_array([1, 2, 3, 4, 5, 6, 7, 8]);
        assert_eq!(v.shift_elements_left_z::<1>(0).to_array(), [2, 3, 4, 5, 6, 7, 8, 0]);
        assert_eq!(v.shift_elements_left_z::<3>(0).to_array(), [4, 5, 6, 7, 8, 0, 0, 0]);
        assert_eq!(v.shift_elements_left_z::<7>(0).to_array(), [8, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(v.shift_elements_left_z::<17>(0).to_array(), [0, 0, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_shift_elements_right() {
        let v = Simd::from_array([1, 2, 3, 4, 5, 6, 7, 8]);
        assert_eq!(v.shift_elements_right_z::<1>(0).to_array(), [0, 1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(v.shift_elements_right_z::<3>(0).to_array(), [0, 0, 0, 1, 2, 3, 4, 5]);
        assert_eq!(v.shift_elements_right_z::<7>(0).to_array(), [0, 0, 0, 0, 0, 0, 0, 1]);
        assert_eq!(v.shift_elements_right_z::<17>(0).to_array(), [0, 0, 0, 0, 0, 0, 0, 0]);
    }
}
