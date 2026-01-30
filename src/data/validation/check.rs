use crate::simd::SimdByteFunctions;
use std::simd::prelude::*;

/// Provides SIMD-accelerated sequence validation methods.
pub trait CheckSequence {
    /// Checks if all bytes in the sequence are ASCII using SIMD operations
    fn is_ascii_simd<const N: usize>(&self) -> bool;

    /// Checks if all bytes in the sequence are printable ASCII using SIMD operations
    fn is_graphic_simd<const N: usize>(&self) -> bool;
}

impl<T> CheckSequence for T
where
    T: AsRef<[u8]>,
{
    #[inline]
    fn is_ascii_simd<const N: usize>(&self) -> bool {
        let (pre, mid, suffix) = self.as_ref().as_simd::<N>();
        pre.is_ascii() && suffix.is_ascii() && mid.iter().fold(Mask::splat(true), |acc, b| acc & b.is_ascii()).all()
    }

    #[inline]
    fn is_graphic_simd<const N: usize>(&self) -> bool {
        let (pre, mid, suffix) = self.as_ref().as_simd::<N>();
        pre.iter().fold(true, |acc, b| acc & b.is_ascii_graphic())
            && suffix.iter().fold(true, |acc, b| acc & b.is_ascii_graphic())
            && mid.iter().fold(Mask::splat(true), |acc, b| acc & b.is_ascii_graphic()).all()
    }
}

#[cfg(all(test, feature = "rand"))]
mod test {
    use super::CheckSequence;
    use crate::data::alphas::AA_DAIS_WITH_GAPS_X;

    #[test]
    fn is_ascii() {
        let s = crate::generate::rand_sequence(AA_DAIS_WITH_GAPS_X, 151, 42);
        assert_eq!(s.is_ascii(), s.is_ascii_simd::<16>());
    }
}

#[cfg(all(test, feature = "rand"))]
mod bench {
    use super::CheckSequence;
    use crate::data::alphas::AA_DAIS_WITH_GAPS_X;
    use std::sync::LazyLock;
    use test::Bencher;
    extern crate test;

    const N: usize = 151;
    const SEED: u64 = 99;

    static SEQ: LazyLock<Vec<u8>> = LazyLock::new(|| crate::generate::rand_sequence(AA_DAIS_WITH_GAPS_X, N, SEED));

    #[bench]
    fn is_ascii_std(b: &mut Bencher) {
        b.iter(|| SEQ.is_ascii());
    }

    #[bench]
    fn is_ascii_zoe(b: &mut Bencher) {
        let (p, m, s) = SEQ.as_simd::<16>();
        eprintln!("{p} {m} {s}", p = p.len(), m = m.len(), s = s.len());
        b.iter(|| SEQ.is_ascii_simd::<16>());
    }
}
