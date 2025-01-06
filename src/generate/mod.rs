#[must_use]
pub fn rand_sequence(alpha: &[u8], length: usize, seed: u64) -> Vec<u8> {
    use rand_xoshiro::{
        Xoshiro256PlusPlus,
        rand_core::{RngCore, SeedableRng},
    };

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    (1..=length).map(|_| alpha[rng.next_u32() as usize % alpha.len()]).collect()
}

#[cfg(test)]
mod test {
    use super::rand_sequence;

    #[test]
    fn rand_test() {
        const LEN: usize = 10_000;

        let random_sequence = rand_sequence(b"ATGC", LEN, 42);
        assert_eq!(LEN, random_sequence.len());

        let (a, c, g, t) = random_sequence.iter().fold((0, 0, 0, 0), |(a, c, g, t), &b| match b {
            b'A' => (a + 1, c, g, t),
            b'C' => (a, c + 1, g, t),
            b'G' => (a, c, g + 1, t),
            b'T' => (a, c, g, t + 1),
            _ => (a, c, g, t),
        });

        assert!(a > 0);
        assert!(c > 0);
        assert!(g > 0);
        assert!(t > 0);
    }
}
