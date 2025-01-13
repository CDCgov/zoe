use zoe::data::records::fasta::*;

#[test]
fn to_reverse() {
    let mut s = FastaSeq {
        name:     "s1".to_string(),
        sequence: b"AtgcnN-".to_vec(),
    };

    s.reverse_complement();
    assert_eq!(String::from_utf8_lossy(&s.sequence), "-NngcaT");
}

#[test]
fn forward_rev_forward() {
    let mut s = FastaSeq {
        name:     "s1".to_string(),
        sequence: b"atgc".to_vec(),
    };

    s.reverse_complement();
    s.reverse_complement();
    assert_eq!(String::from_utf8_lossy(&s.sequence), "atgc");
}
