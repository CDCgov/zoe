#![allow(dead_code, clippy::explicit_iter_loop)]

/// Sorts sequence-like vectors, usually bases, amino acids, or quality scores.
/// Not for use with random bytes or sequences containing UTF or control
/// characters. Only suitable for printable ASCII [32 to 127].
///
/// Performs the prefix sum and output mapping according to the work by Thomas
/// Cormen.
#[inline]
#[must_use]
pub(crate) fn count_sorted_ascii_print(v: &[u8]) -> Vec<u8> {
    if v.len() < 2 {
        v.to_vec()
    } else {
        let mut counts = [0usize; 256];
        let mut output = vec![0u8; v.len()];

        // Create histogram
        for b in v {
            counts[*b as usize] += 1;
        }

        // Perform prefix sum
        for index in 32..=127 {
            counts[index] += counts[index - 1];
        }

        // Map back
        for b in v.iter().rev() {
            let key = *b as usize;

            counts[key] -= 1;
            output[counts[key]] = *b;
        }

        output
    }
}

/// Sorts sequence-like vectors, usually bases, amino acids, or quality scores.
/// Not for use with random bytes or sequences containing UTF or control
/// characters. Only suitable for printable ASCII [32 to 127].
///
/// Uses a nested loop in place of a prefix sum for mapping back onto the
/// original, mutable array. Original algorithm by Thomas Cormen.
#[inline]
pub(crate) fn count_sort_ascii_print(v: &mut [u8]) {
    if v.len() > 1 {
        let mut counts = [0usize; 256];

        // Create histogram
        for b in v.iter() {
            counts[*b as usize] += 1;
        }

        // Re-fill array
        let mut p = 0;
        for i in 32..=127 {
            for _ in 0..counts[i as usize] {
                v[p] = i;
                p += 1;
            }
        }
    }
}

#[cfg(test)]
mod test {
    use crate::generate::rand_sequence;

    #[test]
    fn sort_test() {
        use super::*;

        for _ in 0..1000 {
            let og = rand_sequence(b"ABCDEFGHIJKLMNOPQRSTUVWXYZ", 1000, 42);
            let mut s1 = og.clone();
            let mut s2 = count_sorted_ascii_print(&og);
            s1.sort_unstable();
            assert_eq!(s1, s2);

            s2 = og;
            count_sort_ascii_print(&mut s2);
            assert_eq!(s1, s2);
        }
    }
}

#[allow(clippy::redundant_clone, clippy::stable_sort_primitive)]
#[cfg(test)]
mod bench {
    extern crate test;
    use crate::generate::rand_sequence;
    use test::Bencher;

    const N: usize = 150;
    const ALPHA: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const SEED: u64 = 42;

    #[bench]
    fn clock_rand_clone(b: &mut Bencher) {
        b.iter(|| {
            let v1 = rand_sequence(ALPHA, N, SEED);
            v1.clone()
        });
    }

    #[bench]
    fn sorted_vanilla(b: &mut Bencher) {
        b.iter(|| {
            let v1 = rand_sequence(ALPHA, N, SEED);
            v1.clone().sort();
            v1
        });
    }

    #[bench]
    fn sorted_unstable(b: &mut Bencher) {
        b.iter(|| {
            let v1 = rand_sequence(ALPHA, N, SEED);
            v1.clone().sort_unstable();
        });
    }

    #[bench]
    fn sorted_count(b: &mut Bencher) {
        // Clones within itself
        #[allow(unused)]
        b.iter(|| {
            let v1 = rand_sequence(ALPHA, N, SEED);
            super::count_sorted_ascii_print(&v1);
            v1
        });
    }

    #[bench]
    fn mod_vanilla(b: &mut Bencher) {
        b.iter(|| {
            let mut v1 = rand_sequence(ALPHA, N, SEED);
            v1.sort();
            v1
        });
    }

    #[bench]
    fn mod_unstable(b: &mut Bencher) {
        b.iter(|| {
            let mut v1 = rand_sequence(ALPHA, N, SEED);
            v1.sort_unstable();
            v1
        });
    }

    #[bench]
    fn mod_count(b: &mut Bencher) {
        b.iter(|| {
            let mut v1 = rand_sequence(ALPHA, N, SEED);
            super::count_sort_ascii_print(&mut v1);
            v1
        });
    }
}
