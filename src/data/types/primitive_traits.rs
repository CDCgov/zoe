use std::ops::{Add, AddAssign};
use std::simd::SimdElement;

macro_rules! impl_int {
    { $ident:ident; $($ty:ty),* } => {
        pub trait $ident: SimdElement + Ord + Add + AddAssign {
            fn one() -> Self;
            fn zero() -> Self;
        }

        $(
        impl $ident for $ty {
            fn one() -> Self { 1 }
            fn zero() -> Self { 0 }
        } )*

     }
}

impl_int!(Uint; u8, u16, u32, u64, usize);
impl_int!(Int; u8, u16, u32, u64, usize, i8, i16, i32, i64, isize);
