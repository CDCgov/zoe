extern crate test;
use super::{test_data::*, *};
use crate::alignment::profile::StripedProfile;
use test::Bencher;

#[bench]
fn sw_alignment_scalar(b: &mut Bencher) {
    let query_profile = &*SCALAR_PROFILE;
    b.iter(|| sw_scalar_alignment(REFERENCE, query_profile));
}

#[bench]
fn sw_score_scalar(b: &mut Bencher) {
    let query_profile = &*SCALAR_PROFILE;
    b.iter(|| sw_scalar_score(REFERENCE, query_profile));
}

mod int08 {
    use super::*;

    #[bench]
    fn sw_simd_w0064n08i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i8, 8, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0128n16u(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u8, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0128n16i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i8, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0256n32u(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u8, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0256n32i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i8, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0512n64i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i8, 64, 5>(QUERY, &query_profile));
    }
}

mod int16 {
    use super::*;

    #[bench]
    fn sw_simd_w0064n04i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i16, 4, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0128n08i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i16, 8, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0256n16u(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u16, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0256n16i_score(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i16, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0256n16i_ends(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score_ends::<i16, 16, 5>(QUERY, &query_profile));
    }

    #[cfg(feature = "dev-3pass")]
    #[bench]
    fn sw_simd_w0256n16i_ranges(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score_ranges::<i16, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0512n32u(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &BIASED_WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<u16, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0512n32i_ends(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score_ends::<i16, 32, 5>(QUERY, &query_profile));
    }

    #[cfg(feature = "dev-3pass")]
    #[bench]
    fn sw_simd_w0512n32i_ranges(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score_ranges::<i16, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0512n32i_score(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i16, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w1024n64i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i16, 64, 5>(QUERY, &query_profile));
    }
}

mod int32 {
    use super::*;

    #[bench]
    fn sw_simd_w0064n02i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i32, 2, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0128n04i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i32, 4, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0256n08i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i32, 8, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w0512n16i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i32, 16, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w1024n32i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i32, 32, 5>(QUERY, &query_profile));
    }

    #[bench]
    fn sw_simd_w2048n64i(b: &mut Bencher) {
        let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
        b.iter(|| sw_simd_score::<i32, 64, 5>(QUERY, &query_profile));
    }
}

mod simd_aln {
    use super::*;

    mod int08 {
        use super::*;
        #[bench]
        fn sw_aln_simd_w128n16i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i8, 16, 5>(QUERY, &query_profile));
        }

        #[bench]
        fn sw_aln_simd_w256n32i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i8, 32, 5>(QUERY, &query_profile));
        }

        #[bench]
        fn sw_aln_simd_w512n64i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i8, 64, 5>(QUERY, &query_profile));
        }
    }

    mod int16 {
        use super::*;
        #[bench]
        fn sw_aln_simd_w0256n16i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i16, 16, 5>(QUERY, &query_profile));
        }

        #[bench]
        fn sw_aln_simd_w0512n32i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i16, 32, 5>(QUERY, &query_profile));
        }

        #[bench]
        fn sw_aln_simd_w1024n64i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i16, 64, 5>(QUERY, &query_profile));
        }
    }

    mod int32 {
        use super::*;
        #[bench]
        fn sw_aln_simd_w0256n08i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i32, 8, 5>(QUERY, &query_profile));
        }

        #[bench]
        fn sw_aln_simd_w0512n16i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i32, 16, 5>(QUERY, &query_profile));
        }

        #[bench]
        fn sw_aln_simd_w1024n32i(b: &mut Bencher) {
            let query_profile = StripedProfile::new(REFERENCE, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
            b.iter(|| sw_simd_alignment::<i32, 32, 5>(QUERY, &query_profile));
        }
    }
}

mod banded {
    extern crate test;
    use super::*;
    use test::Bencher;

    #[bench]
    fn sw_banded_alignment_w03(b: &mut Bencher) {
        let query_profile = &*SCALAR_PROFILE;
        b.iter(|| sw_banded_alignment(REFERENCE, query_profile, 3));
    }

    #[bench]
    fn sw_banded_alignment_w05(b: &mut Bencher) {
        let query_profile = &*SCALAR_PROFILE;
        b.iter(|| sw_banded_alignment(REFERENCE, query_profile, 5));
    }

    #[bench]
    fn sw_banded_alignment_w10(b: &mut Bencher) {
        let query_profile = &*SCALAR_PROFILE;
        b.iter(|| sw_banded_alignment(REFERENCE, query_profile, 10));
    }

    #[bench]
    fn sw_banded_alignment_w20(b: &mut Bencher) {
        let query_profile = &*SCALAR_PROFILE;
        b.iter(|| sw_banded_alignment(REFERENCE, query_profile, 20));
    }

    #[bench]
    fn sw_banded_alignment_w50(b: &mut Bencher) {
        let query_profile = &*SCALAR_PROFILE;
        b.iter(|| sw_banded_alignment(REFERENCE, query_profile, 50));
    }
}
