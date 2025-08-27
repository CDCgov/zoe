use crate::{
    alignment::{AlignmentStates, PairwiseSequence, StatesSequence, profile::ScalarProfile, sw::sw_scalar_alignment},
    data::{
        amino_acids::AminoAcids,
        matrices::WeightMatrix,
        types::cigar::{Cigar, Ciglet},
    },
};

#[test]
fn alignment_invert() {
    const WEIGHTS: WeightMatrix<i8, 5> = WeightMatrix::new_dna_matrix(4, -2, Some(b'N'));
    const GAP_OPEN: i8 = -3;
    const GAP_EXTEND: i8 = -1;

    let reference: &[u8] = b"GGCCACAGGATTGAGC";
    let query: &[u8] = b"TCTCAGATTGCAGTTT";

    let profile = ScalarProfile::<5>::new(query, &WEIGHTS, GAP_OPEN, GAP_EXTEND).unwrap();
    let alignment = sw_scalar_alignment(reference, &profile).unwrap();
    let invert_alignment = alignment.invert();

    assert_eq!(alignment.ref_range, 3..15);
    assert_eq!(alignment.query_range, 1..13);
    assert_eq!(alignment.states, Cigar::from_slice_unchecked("1S5M1D4M1I2M3S"));

    assert_eq!(invert_alignment.ref_range, 1..13);
    assert_eq!(invert_alignment.query_range, 3..15);
    assert_eq!(invert_alignment.states, Cigar::from_slice_unchecked("3S5M1I4M1D2M1S"));
}

#[test]
fn states_sequence() {
    // The clipping is not valid here, but we are using it to exercise all
    // options
    let cigar = Cigar::try_from("3H4S4S3H10M9I3M1D5M3H5S10S1H").unwrap();
    let alignment_states = cigar.iter().collect::<AlignmentStates>();
    let mut cigar_slice = alignment_states.as_slice();
    let mut ciglet_iterator = cigar.into_iter();

    assert_eq!(cigar_slice.peek_op(), Some(b'H'));
    assert_eq!(ciglet_iterator.peek_op(), Some(b'H'));
    assert_eq!(cigar_slice.peek_back_op(), Some(b'H'));
    assert_eq!(ciglet_iterator.peek_back_op(), Some(b'H'));

    // Shouldn't have changed, since we didn't advance
    assert_eq!(cigar_slice.peek_op(), Some(b'H'));
    assert_eq!(ciglet_iterator.peek_op(), Some(b'H'));
    assert_eq!(cigar_slice.peek_back_op(), Some(b'H'));
    assert_eq!(ciglet_iterator.peek_back_op(), Some(b'H'));

    assert_eq!(cigar_slice.remove_clipping_front(), 7);
    assert_eq!(ciglet_iterator.remove_clipping_front(), 7);
    assert_eq!(cigar_slice.remove_clipping_back(), 11);
    assert_eq!(ciglet_iterator.remove_clipping_back(), 11);

    assert_eq!(cigar_slice.peek_op(), Some(b'S'));
    assert_eq!(ciglet_iterator.peek_op(), Some(b'S'));
    assert_eq!(cigar_slice.peek_back_op(), Some(b'S'));
    assert_eq!(ciglet_iterator.peek_back_op(), Some(b'S'));

    assert_eq!(cigar_slice.remove_clipping_front(), 4);
    assert_eq!(ciglet_iterator.remove_clipping_front(), 4);
    assert_eq!(cigar_slice.remove_clipping_back(), 5);
    assert_eq!(ciglet_iterator.remove_clipping_back(), 5);

    assert_eq!(cigar_slice.remove_clipping_front(), 3);
    assert_eq!(ciglet_iterator.remove_clipping_front(), 3);
    assert_eq!(cigar_slice.remove_clipping_back(), 3);
    assert_eq!(ciglet_iterator.remove_clipping_back(), 3);

    assert_eq!(cigar_slice.remove_clipping_front(), 0);
    assert_eq!(ciglet_iterator.remove_clipping_front(), 0);
    assert_eq!(cigar_slice.remove_clipping_back(), 0);
    assert_eq!(ciglet_iterator.remove_clipping_back(), 0);

    assert_eq!(cigar_slice.next_if_op(|op| op == b'I'), None);
    assert_eq!(ciglet_iterator.next_if_op(|op| op == b'I'), None);
    assert_eq!(cigar_slice.next_back_if_op(|op| op == b'I'), None);
    assert_eq!(ciglet_iterator.next_back_if_op(|op| op == b'I'), None);

    assert_eq!(cigar_slice.next_if_op(|op| op == b'M'), Some(Ciglet { inc: 10, op: b'M' }));
    assert_eq!(
        ciglet_iterator.next_if_op(|op| op == b'M'),
        Some(Ciglet { inc: 10, op: b'M' })
    );
    assert_eq!(
        cigar_slice.next_back_if_op(|op| op == b'M'),
        Some(Ciglet { inc: 5, op: b'M' })
    );
    assert_eq!(
        ciglet_iterator.next_back_if_op(|op| op == b'M'),
        Some(Ciglet { inc: 5, op: b'M' })
    );

    assert_eq!(cigar_slice.next_ciglet(), Some(Ciglet { inc: 9, op: b'I' }));
    assert_eq!(ciglet_iterator.next_ciglet(), Some(Ciglet { inc: 9, op: b'I' }));
    assert_eq!(cigar_slice.next_ciglet_back(), Some(Ciglet { inc: 1, op: b'D' }));
    assert_eq!(ciglet_iterator.next_ciglet_back(), Some(Ciglet { inc: 1, op: b'D' }));

    assert!(!cigar_slice.is_empty());
    assert!(!ciglet_iterator.is_empty());

    assert_eq!(cigar_slice.next_ciglet(), Some(Ciglet { inc: 3, op: b'M' }));
    assert_eq!(ciglet_iterator.next_ciglet(), Some(Ciglet { inc: 3, op: b'M' }));
    assert_eq!(cigar_slice.next_ciglet_back(), None);
    assert_eq!(ciglet_iterator.next_ciglet_back(), None);
    assert_eq!(cigar_slice.next_ciglet(), None);
    assert_eq!(ciglet_iterator.next_ciglet(), None);

    assert!(cigar_slice.is_empty());
    assert!(ciglet_iterator.is_empty());
}

#[test]
fn align_with_cigar() {
    let data: [(usize, Cigar, [AminoAcids; 4]); 1] = [(
        2,
        Cigar::from_slice_unchecked("4M"),
        [b"PLEASANTLY".into(), b"MEANLY".into(), b"LEAS".into(), b"MEAN".into()],
    )];

    for (rpos, cig, [ref_input, query_input, ref_output, query_ouptut]) in data {
        assert_eq!(
            ref_input.align_and_collect(&query_input, &cig, rpos),
            (ref_output, query_ouptut)
        );
    }
}

#[test]
fn test_eq_states_cigar() {
    let mut states = AlignmentStates::new();
    states.add_ciglet(Ciglet { inc: 10, op: b'M' });

    let cigar = Cigar::from_slice_unchecked(b"10M3D");
    assert_ne!(states, cigar);
    assert_ne!(cigar, states);

    let cigar = Cigar::from_slice_unchecked(b"10M$$");
    assert_ne!(states, cigar);
    assert_ne!(cigar, states);

    let cigar = Cigar::from_slice_unchecked(b"10M");
    assert_eq!(states, cigar);
    assert_eq!(cigar, states);

    let cigar = Cigar::new();
    let states = AlignmentStates::new();
    assert_eq!(states, cigar);
    assert_eq!(cigar, states);
}
