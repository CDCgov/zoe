use super::*;
use crate::alignment::AlignmentStates;

#[test]
fn test_expand() {
    let cigar = Cigar::from_slice_unchecked("4S10M2I2D3M4H4P");
    let expanded: ExpandedCigar = "SSSSMMMMMMMMMMIIDDMMMHHHHPPPP".into();
    assert_eq!(cigar.expand_cigar(), expanded);
}

#[test]
fn test_iter() {
    let cigar = Cigar::from_slice_unchecked("1M22I333D4444N55555S6666H777P88X9=");
    let mut cigar = cigar.into_iter();

    assert_eq!(cigar.next(), Some(Ciglet { op: b'M', inc: 1 }));
    assert_eq!(cigar.next(), Some(Ciglet { op: b'I', inc: 22 }));
    assert_eq!(cigar.next(), Some(Ciglet { op: b'D', inc: 333 }));
    assert_eq!(cigar.next(), Some(Ciglet { op: b'N', inc: 4444 }));
    assert_eq!(cigar.next(), Some(Ciglet { op: b'S', inc: 55555 }));
    assert_eq!(cigar.next(), Some(Ciglet { op: b'H', inc: 6666 }));
    assert_eq!(cigar.next(), Some(Ciglet { op: b'P', inc: 777 }));
    assert_eq!(cigar.next(), Some(Ciglet { op: b'X', inc: 88 }));
    assert_eq!(cigar.next(), Some(Ciglet { op: b'=', inc: 9 }));
    assert_eq!(cigar.next(), None);

    // Illegal cigar op
    let cigar = Cigar::from_slice_unchecked("8K");
    assert_eq!(cigar.into_iter().next(), None);

    // Leading zeroes don't matter
    let cigar = Cigar::from_slice_unchecked("000000000000000000000000000000155M");
    assert_eq!(cigar.into_iter().next(), Some(Ciglet { op: b'M', inc: 155 }));

    // Overflows
    let cigar = Cigar::from_slice_unchecked("100000000000000000000000000000155M");
    assert_eq!(cigar.into_iter().next(), None);

    // Bad order
    let cigar = Cigar::from_slice_unchecked("M155M");
    assert_eq!(cigar.into_iter().next(), None);

    // usize == u64
    if USIZE_WIDTH == 20 {
        let cigar = Cigar::from_slice_unchecked("18446744073709551615M");
        assert_eq!(
            cigar.into_iter().next(),
            Some(Ciglet {
                op:  b'M',
                inc: 18_446_744_073_709_551_615,
            })
        );

        let cigar = Cigar::from_slice_unchecked("001234567890123456789M");
        assert_eq!(
            cigar.into_iter().next(),
            Some(Ciglet {
                op:  b'M',
                inc: 1_234_567_890_123_456_789,
            })
        );
    }
}

#[test]
fn test_ref_len_in_alignment() {
    let cigars = [
        ("4S10M2I2D3M4H4P", 15),
        ("", 0),
        ("3M2D1M", 6),
        ("255M", 255),
        ("3M1D4I8X9=4M", 25),
        ("M", 0),
    ];

    for (c, l) in cigars {
        let cigar = Cigar::from_slice_unchecked(c.as_bytes());
        assert_eq!(cigar.ref_len_in_alignment(), l);
    }
}

#[test]
fn test_query_len_in_alignment() {
    let cigars = [
        ("4S10M2I2D3M4H4P", 19),
        ("", 0),
        ("3M2D1M", 4),
        ("255M", 255),
        ("3M1D4I8X9=4M", 28),
        ("M", 0),
    ];

    for (c, l) in cigars {
        let cigar = Cigar::from_slice_unchecked(c.as_bytes());
        assert_eq!(cigar.query_len_in_alignment(), l);
    }
}

#[test]
fn test_condense_cigar() {
    let cigars: [&str; 5] = ["4S10M2I2D3M4H4P", "", "3M2D1M", "255M", "3M1D4I8X9=4M"];

    for c in cigars {
        let cigar = Cigar::from_slice_unchecked(c);
        assert_eq!(cigar.expand_cigar().condense_to_cigar(), cigar);
    }
}

#[test]
fn test_add_state() {
    let cigar = Cigar::from_slice_unchecked("4S10M2I2D3M4H4P");
    let mut states = AlignmentStates::new();
    for Ciglet { inc, op } in &cigar {
        for _ in 0..inc {
            states.add_state(op);
        }
    }

    assert_eq!(cigar, states.to_cigar_unchecked());
}

#[test]
fn test_is_valid() {
    assert!(Cigar::from_slice_unchecked("10M5I20D").is_valid());
    assert!(!Cigar::from_slice_unchecked("M10D").is_valid());
    assert!(!Cigar::from_slice_unchecked("10M5I0D").is_valid());
    assert!(!Cigar::from_slice_unchecked("10MM5I").is_valid());
    assert!(Cigar::try_from("19446744073709551616M").is_err());
}

#[test]
fn test_try_from() {
    assert!(Cigar::try_from("10M5I20D").is_ok());
    assert!(Cigar::try_from("1M").is_ok());
    assert!(Cigar::try_from("1000000M").is_ok());
    assert!(Cigar::try_from("*").is_ok());
    assert_eq!(Cigar::try_from("10M5I20X@"), Err(CigarError::InvalidOperation));
    assert_eq!(Cigar::try_from("10M#I5D"), Err(CigarError::InvalidOperation));
    assert_eq!(Cigar::try_from("10M5I20"), Err(CigarError::MissingOp));
    assert_eq!(Cigar::try_from("M10D"), Err(CigarError::MissingInc));
    assert_eq!(Cigar::try_from("10MI5D"), Err(CigarError::MissingInc));
    assert_eq!(Cigar::try_from("MI"), Err(CigarError::MissingInc));
    assert_eq!(Cigar::try_from("10M0D"), Err(CigarError::IncZero));
    assert_eq!(Cigar::try_from("0M"), Err(CigarError::IncZero));
    assert_eq!(Cigar::try_from("18446744073709551616M"), Err(CigarError::IncOverflow));
    assert_eq!(Cigar::try_from(""), Ok(Cigar::new()));
    assert_eq!(Cigar::try_from("10M5I20*5D"), Err(CigarError::InvalidOperation));
    assert!(Cigar::try_from("10M3I2M4M").is_ok());
    assert!(Cigar::try_from("10S6S").is_ok());
}
