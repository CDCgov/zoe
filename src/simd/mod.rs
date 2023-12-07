use std::simd::{prelude::*, LaneCount, SupportedLaneCount};

pub(crate) trait SimdByteFunctions<const N: usize>
where
    LaneCount<N>: SupportedLaneCount, {
    fn is_ascii_uppercase(&self) -> Mask<i8, N>;
    fn is_ascii_lowercase(&self) -> Mask<i8, N>;

    fn to_ascii_uppercase(&self) -> Self;
    fn to_ascii_lowercase(&self) -> Self;

    fn make_ascii_uppercase(&mut self);
    fn make_ascii_lowercase(&mut self);

    fn if_value_then_replace(&mut self, find: u8, replace: u8);
    fn swap_byte_pairs(&mut self, this: u8, that: u8);
}

impl<const N: usize> SimdByteFunctions<N> for Simd<u8, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn is_ascii_uppercase(&self) -> Mask<i8, N>
    where
        LaneCount<N>: SupportedLaneCount, {
        self.simd_ge(Simd::splat(b'A')) & self.simd_le(Simd::splat(b'Z'))
    }

    #[inline]
    fn is_ascii_lowercase(&self) -> Mask<i8, N>
    where
        LaneCount<N>: SupportedLaneCount, {
        self.simd_ge(Simd::splat(b'a')) & self.simd_le(Simd::splat(b'z'))
    }

    #[inline]
    fn to_ascii_uppercase(&self) -> Self
    where
        LaneCount<N>: SupportedLaneCount, {
        let mask = self.is_ascii_lowercase();
        mask.select(*self ^ Simd::splat(0b0010_0000), *self)
    }

    #[inline]
    fn to_ascii_lowercase(&self) -> Self
    where
        LaneCount<N>: SupportedLaneCount, {
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
    fn if_value_then_replace(&mut self, find: u8, replace: u8)
    where
        LaneCount<N>: SupportedLaneCount, {
        let mask = self.simd_eq(Simd::splat(find));
        *self = mask.select(Simd::splat(replace), *self);
    }

    #[inline]
    fn swap_byte_pairs(&mut self, this: u8, that: u8)
    where
        LaneCount<N>: SupportedLaneCount, {
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
