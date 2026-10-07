use crate::simd::SimdByteFunctions;
use std::{fmt::Display, simd::prelude::*};

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

/// An extension trait for byte data allowing it to be displayed efficiently
/// using SIMD. Specifically, the data is checked using SIMD to see whether it
/// contains ASCII, in which case it is displayed directly. Otherwise, it is
/// lossily converted to a [`String`].
pub(crate) trait SimdDisplay: AsRef<[u8]> {
    /// Displays the bytes efficiently if it solely contains ASCII, otherwise
    /// falls back to an allocating lossy display.
    fn display_ascii_or_lossy(&self) -> AsciiOrLossyDisplay<'_> {
        AsciiOrLossyDisplay(self.as_ref())
    }
}

impl<T: AsRef<[u8]>> SimdDisplay for T {}

/// A wrapper type around a byte slice allowing it to be displayed efficiently
/// using SIMD. Specifically, the data is checked using SIMD to see whether it
/// contains ASCII, in which case it is displayed directly. Otherwise, it is
/// lossily converted to a [`String`].
pub(crate) struct AsciiOrLossyDisplay<'a>(&'a [u8]);

impl Display for AsciiOrLossyDisplay<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.0.is_ascii_simd::<16>() {
            // SAFETY: we just checked it is ASCII using our fast SIMD function.
            // ASCII is valid UTF8.
            f.write_str(unsafe { std::str::from_utf8_unchecked(self.0) })
        } else {
            f.write_str(&String::from_utf8_lossy(self.0))
        }
    }
}

#[cfg(all(test, feature = "rand"))]
mod test {
    use super::CheckSequence;
    use crate::data::validation::SimdDisplay;

    /// DAIS-Ribosome style amino acid codes: IUPAC + gaps + X + partial codons `~`.
    pub(crate) const AA_DAIS_WITH_GAPS_X: &[u8; 45] = b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.~Xx";

    #[test]
    fn is_ascii() {
        let s = crate::generate::rand_sequence(AA_DAIS_WITH_GAPS_X, 151, 42);
        assert_eq!(s.is_ascii(), s.is_ascii_simd::<16>());
    }

    #[test]
    fn display_ascii_or_lossy() {
        let s = crate::generate::rand_sequence(AA_DAIS_WITH_GAPS_X, 151, 42);
        assert_eq!(s.display_ascii_or_lossy().to_string(), str::from_utf8(&s).unwrap());

        let s = vec![255, 255, 255, 255];
        assert_eq!(
            s.display_ascii_or_lossy().to_string(),
            String::from_utf8_lossy(&s).to_string()
        );
    }
}

#[cfg(all(test, feature = "rand"))]
mod bench {
    use super::CheckSequence;
    use crate::data::validation::check::test::AA_DAIS_WITH_GAPS_X;
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
