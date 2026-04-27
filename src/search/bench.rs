extern crate test;
use crate::search::{
    ByteSubstring, find_k_repeating, find_k_repeating_scalar, position_by_byte, position_by_byte2, replace_all_bytes,
    replace_all_bytes_simd, substring_match, substring_match_simd,
};
use std::sync::LazyLock;
use test::{Bencher, black_box};

mod find_repeating_at_ends {
    use super::*;

    static SEQ: &[u8] = b"GGGGGGGGGGGGAGCAAGCACAAAACAAGTTAAAGTTACTGGCCATAACAGCCAGAGGAAAATTAACTTAATTATATACAAAAACATATTCCTGTTGGCATAGGCAAATTTTAGAAGACAAATCCATGTAAGGAATAGGGGGGGGGGGGGG";
    static LONG_SEQ: &[u8] = b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGCAAGCACAAAACAAGTTAAAGTTACTGGCCATAACAGCCAGAGGAAAATTAACTTAATTATATACAAAAACATATTCCTGTTGGCATAGGCAAATTTTAGAAGACAAATCCATGTAAGGAATAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

    mod short {
        use super::*;

        #[bench]
        fn bench_starts_with_repeating(b: &mut Bencher) {
            b.iter(|| {
                for _ in 0..10 {
                    let ans = black_box(SEQ).find_repeating_at_start(b'G', 10).is_some()
                        && black_box(SEQ).find_repeating_at_end(b'G', 10).is_some();
                    black_box(ans);
                }
            });
        }

        #[bench]
        fn bench_check_literal(b: &mut Bencher) {
            let needle = vec![b'G'; 10];
            b.iter(|| {
                for _ in 0..10 {
                    let ans = if black_box(SEQ).starts_with(&needle) {
                        let offset = SEQ[10..].iter().take_while(|b| **b == b'G').count();
                        Some(0..10 + offset)
                    } else {
                        None
                    }
                    .is_some()
                        && if black_box(SEQ).ends_with(&needle) {
                            let offset = SEQ[..SEQ.len() - 10].iter().rev().take_while(|b| **b == b'G').count();
                            Some(10 - offset..SEQ.len())
                        } else {
                            None
                        }
                        .is_some();
                    black_box(ans);
                }
            });
        }
    }

    mod long {
        use super::*;

        #[bench]
        fn bench_starts_with_repeating(b: &mut Bencher) {
            b.iter(|| {
                for _ in 0..10 {
                    let ans = black_box(LONG_SEQ).find_repeating_at_start(b'G', 100).is_some()
                        && black_box(LONG_SEQ).find_repeating_at_end(b'G', 100).is_some();
                    black_box(ans);
                }
            });
        }

        #[bench]
        fn bench_check_literal(b: &mut Bencher) {
            let needle = vec![b'G'; 100];
            b.iter(|| {
                for _ in 0..10 {
                    let ans = if black_box(LONG_SEQ).starts_with(&needle) {
                        let offset = LONG_SEQ[100..].iter().take_while(|b| **b == b'G').count();
                        Some(0..100 + offset)
                    } else {
                        None
                    }
                    .is_some()
                        && if black_box(LONG_SEQ).ends_with(&needle) {
                            let offset = LONG_SEQ[..LONG_SEQ.len() - 100]
                                .iter()
                                .rev()
                                .take_while(|b| **b == b'G')
                                .count();
                            Some(100 - offset..LONG_SEQ.len())
                        } else {
                            None
                        }
                        .is_some();
                    black_box(ans);
                }
            });
        }
    }
}

mod k_repeating {
    use super::*;

    extern crate test;

    static WORST_REPEATING: LazyLock<Vec<u8>> = LazyLock::new(|| {
        let mut s = b"abb".to_vec();
        s.extend(b"aabb".repeat(398));
        s.extend(b"aabbb");
        s
    });

    static AVG_REPEATING: LazyLock<Vec<u8>> = LazyLock::new(|| {
        let mut s = b"aaaaaaab".repeat(199);
        s.extend(b"aaaaabbb");
        s
    });

    mod average {
        use super::*;

        #[bench]
        fn using_substring_match(b: &mut Bencher) {
            let needle = vec![b'b'; 3];
            b.iter(|| substring_match_simd::<16>(&AVG_REPEATING, &needle));
        }

        #[bench]
        fn composite(b: &mut Bencher) {
            b.iter(|| find_k_repeating::<16>(&AVG_REPEATING, b'b', 3));
        }

        #[bench]
        fn simd(b: &mut Bencher) {
            b.iter(|| find_k_repeating::<16>(&AVG_REPEATING, b'b', 3));
        }

        #[bench]
        fn scalar(b: &mut Bencher) {
            b.iter(|| find_k_repeating_scalar(&AVG_REPEATING, b'b', 3));
        }
    }

    mod worst {
        use super::*;

        #[bench]
        fn using_substring_match(b: &mut Bencher) {
            let needle = vec![b'b'; 3];
            b.iter(|| substring_match_simd::<16>(&WORST_REPEATING, &needle));
        }

        #[bench]
        fn composite(b: &mut Bencher) {
            b.iter(|| find_k_repeating::<16>(&WORST_REPEATING, b'b', 3));
        }

        #[bench]
        fn simd(b: &mut Bencher) {
            b.iter(|| find_k_repeating::<16>(&WORST_REPEATING, b'b', 3));
        }

        #[bench]
        fn scalar(b: &mut Bencher) {
            b.iter(|| find_k_repeating_scalar(&WORST_REPEATING, b'b', 3));
        }
    }
}

#[cfg(feature = "rand")]
mod replace {
    use super::*;
    use crate::prelude::rand_sequence;

    static LONG: LazyLock<Vec<u8>> = LazyLock::new(|| rand_sequence(b"AGCT", 15_000, 42));
    static SHORT: LazyLock<Vec<u8>> = LazyLock::new(|| rand_sequence(b"AGCT", 150, 42));

    #[bench]
    fn find_replace_long_scalar(b: &mut Bencher) {
        b.iter(|| replace_all_bytes(&mut LONG.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_short_scalar(b: &mut Bencher) {
        b.iter(|| replace_all_bytes(&mut SHORT.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_long_simd32(b: &mut Bencher) {
        b.iter(|| replace_all_bytes_simd::<32>(&mut LONG.clone(), b'A', b'T'));
    }

    #[bench]
    fn find_replace_short_simd32(b: &mut Bencher) {
        b.iter(|| replace_all_bytes_simd::<32>(&mut SHORT.clone(), b'A', b'T'));
    }
}

mod byte_position {
    use super::*;

    static LONG: LazyLock<Vec<u8>> = LazyLock::new(|| b"0".repeat(288).into_iter().chain(std::iter::once(b'1')).collect());
    static SHORT: LazyLock<Vec<u8>> = LazyLock::new(|| b"0".repeat(11).into_iter().chain(std::iter::once(b'1')).collect());

    #[bench]
    fn position_short_scalar(b: &mut Bencher) {
        b.iter(|| SHORT.iter().position(|x| *x == b'1'));
    }

    #[bench]
    fn position_long_scalar(b: &mut Bencher) {
        b.iter(|| LONG.iter().position(|x| *x == b'1'));
    }

    #[bench]
    fn position_short_simd(b: &mut Bencher) {
        b.iter(|| position_by_byte::<16>(&SHORT, b'1'));
    }

    #[bench]
    fn position_long_simd(b: &mut Bencher) {
        b.iter(|| position_by_byte::<16>(&LONG, b'1'));
    }
}

mod byte_position2 {
    use super::*;

    static LONG_INNER: [u8; 289] = {
        let mut out = [b'0'; 289];
        *out.last_mut().unwrap() = b'1';
        out
    };

    static SHORT_INNER: [u8; 12] = {
        let mut out = [b'0'; 12];
        *out.last_mut().unwrap() = b'1';
        out
    };

    static LONG: &[u8] = &LONG_INNER;
    static SHORT: &[u8] = &SHORT_INNER;

    #[bench]
    fn position2_short_scalar(b: &mut Bencher) {
        b.iter(|| {
            black_box(SHORT)
                .iter()
                .position(|x| *x == black_box(b'1') || *x == black_box(b'2'))
        });
    }

    #[bench]
    fn position2_long_scalar(b: &mut Bencher) {
        b.iter(|| {
            black_box(LONG)
                .iter()
                .position(|x| *x == black_box(b'1') || *x == black_box(b'2'))
        });
    }

    #[bench]
    fn position2_short_simd(b: &mut Bencher) {
        b.iter(|| position_by_byte2::<16>(black_box(SHORT), black_box(b'1'), black_box(b'2')));
    }

    #[bench]
    fn position2_long_simd(b: &mut Bencher) {
        b.iter(|| position_by_byte2::<16>(black_box(LONG), black_box(b'1'), black_box(b'2')));
    }
}

mod find_substring {
    use super::*;

    static SEQ: &[u8] = b"GCGGAAATCCGCCACGAATGAGAATGTATTTCCCCGACAATCATAATGGGGCGCTCCTAAGCTTTTCCACTTGGTTGGGCCGGCTAGGCCTCTCTGCCCGGAGTTTCGGCGCACTGCTGCCGACAGCCGGGCATTGTTTTAGGGGCGTTA";
    static LONG_SEQ: &[u8] = b"TTCGAGGGCACTCGGAGCTAACTTGTCGGGACCAGCCGGGGTAGTCATCGGGCTTATACAGCGAAAAGCCCAGGACCCGGCTCCACGCTATGGAACGTCTTTAGCTCCGGCAAGCAATTAAGAACAACGCAAGCATCGCGGATATAAACAGAGAAACGGCCGAATACACCTGTTCGTATCGTATCGGTAAATAGCCTCGCGGAGCCATGTGCCATACTGGTCTGCGGAGCACTCTGGTTATGCATATGGTCCACAGGACACTCGTCGCTTCCGGGTATGCGCTCTATGTGACGGTCTTTAGGCGCACTAATGCTCAGCACCATTTAAACCAGACCGACACCAGATCTGTAAGGTCCGCCACGCAGACGACAGCCCACGGAGATCACCGACCGATCTATCTGATCGGCGACCATTTGTGTGGTACTGGGGCGGAGAGGTAACTACGGTGCCGCTAACAACCCCTCTGTCGTCGCTGACGTTTGTAGTCTAGTCTCATTATGATTGTACGCTATTCAGGGATTGACTGATACCGGAAGACATCTCAGTTGAAGTGGTCTATACGACAGAGACCGTGCACCTACCAAATCTCCTTAGTGTAAGTTCAGACCAATTGGTAGTTTGTCCAGAACTCAGATTTTAACAGCAGAGGACGCATGCTCTACCTTCATGATCCACTGACGTCCCTGAGGCTGCAATACATGCAACGAGGCAGTCTCCGCGGTAAGTCCTAGTGCAATGGCGCTTTTTTACCCTCGTCCTCGAGAAGAGGGGACGCCAGTGCAGATATCTTTAATGTGGTAATTGGGAGGACTCTTGGCCCTCCGCCCTTAGGCAGTGCATACTCTTCCATGCGGAAATCCGCCACGAATGAGAATGTATTTCCCCGACAATCATAATGGGGCGCTCCTAAGCTTTTCCACTTGGTTGGGCCGGCTAGGCCTCTCTGCCCGGAGTTTCGGCGCACTGCTGCCGACAGCCGGGCATTGTTTTAGGGGCGTTA";

    #[bench]
    fn substring_match_scalar_short(b: &mut Bencher) {
        dbg!(substring_match(black_box(SEQ), black_box(b"CGTTA")));
        b.iter(|| substring_match(black_box(SEQ), black_box(b"CGTTA")));
    }

    #[bench]
    fn substring_match_scalar_long(b: &mut Bencher) {
        dbg!(substring_match(black_box(LONG_SEQ), black_box(b"CGTTA")));
        b.iter(|| substring_match(black_box(LONG_SEQ), black_box(b"CGTTA")));
    }

    #[bench]
    fn substring_match_simd_short(b: &mut Bencher) {
        dbg!(substring_match_simd::<32>(black_box(SEQ), black_box(b"CGTTA")));
        b.iter(|| substring_match_simd::<32>(black_box(SEQ), black_box(b"CGTTA")));
    }

    #[bench]
    fn substring_match_simd_long(b: &mut Bencher) {
        dbg!(substring_match_simd::<32>(black_box(LONG_SEQ), black_box(b"CGTTA")));
        b.iter(|| substring_match_simd::<32>(black_box(LONG_SEQ), black_box(b"CGTTA")));
    }
}
