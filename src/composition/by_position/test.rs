use super::*;

#[test]
fn empty() {
    let max_base = Nucleotides::from(String::new()).to_base_counts::<u64>().plurality_acgtn();
    assert_eq!(max_base, b'A');
}

#[test]
fn no_valid_bases() {
    let max_base = Nucleotides::from("########".to_string())
        .to_base_counts::<u64>()
        .plurality_acgtn();
    assert_eq!(max_base, b'A');
}

#[test]
fn all_gap_or_other() {
    let max_base = Nucleotides::from("--SW---KMB--".to_string())
        .to_base_counts::<u64>()
        .plurality_acgtn();
    assert_eq!(max_base, b'A');
}

#[test]
fn all_gap_other_or_invalid() {
    let max_base = Nucleotides::from("###--S##W---KM###B--##".to_string())
        .to_base_counts::<u64>()
        .plurality_acgtn();
    assert_eq!(max_base, b'A');
}

#[test]
fn all_unknown() {
    let max_base = Nucleotides::from("NNNNN".to_string())
        .to_base_counts::<u64>()
        .plurality_acgtn();
    assert_eq!(max_base, b'N');
}

#[test]
fn all_tied() {
    let max_base = Nucleotides::from("ACGTN".to_string())
        .to_base_counts::<u64>()
        .plurality_acgtn();
    assert_eq!(max_base, b'A');
}

#[test]
fn test1() {
    let max_base = Nucleotides::from("ACCCCCCCCGGGGGTTN".to_string())
        .to_base_counts::<u64>()
        .plurality_acgtn();
    assert_eq!(max_base, b'C');
}

#[test]
fn test2() {
    let max_base = Nucleotides::from("AACGGGTNNN".to_string())
        .to_base_counts::<u64>()
        .plurality_acgtn();
    assert_eq!(max_base, b'G');
}

#[test]
fn test3() {
    let max_base = Nucleotides::from("##CCC##--ANGT--ACAAACCCC####GCTAAGGGGGGG".to_string())
        .to_base_counts::<u64>()
        .plurality_acgtn();
    assert_eq!(max_base, b'C');
}
