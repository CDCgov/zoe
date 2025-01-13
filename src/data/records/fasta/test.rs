use super::*;
use std::io::Cursor;

#[test]
fn sequence_whitespace() {
    let mut reader = FastaReader::new(Cursor::new(">seq1\nATG C   \n>seq2\r\n\n  AC G\r\n T"));

    let record1 = reader.next().unwrap().unwrap();
    assert_eq!(record1.name, "seq1");
    assert_eq!(record1.sequence, b"ATG C   ");

    let record2 = reader.next().unwrap().unwrap();
    assert_eq!(record2.name, "seq2");
    assert_eq!(record2.sequence, b"  AC G T");

    assert!(reader.next().is_none());
}

#[test]
fn empty_file() {
    let mut reader = FastaReader::new(Cursor::new(""));

    assert_eq!(reader.next().unwrap().unwrap_err().to_string(), "No FASTA data found.");

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn whitespace_only() {
    let mut reader = FastaReader::new(Cursor::new("   \r\n \r\t\n   \t"));

    assert_eq!(reader.next().unwrap().unwrap_err().to_string(), "No FASTA data found.");

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn missing_header() {
    let mut reader = FastaReader::new(Cursor::new("ATGC"));

    assert_eq!(
        reader.next().unwrap().unwrap_err().to_string(),
        "The FASTA file must start with a '>' symbol!"
    );

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_header_first_record() {
    let mut reader = FastaReader::new(Cursor::new(">\nATGC"));

    assert_eq!(reader.next().unwrap().unwrap_err().to_string(), "Missing FASTA header!");

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_sequence_first_record() {
    let mut reader = FastaReader::new(Cursor::new(">seq1\n"));

    assert_eq!(reader.next().unwrap().unwrap_err().to_string(), "Missing FASTA sequence!");

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_header_second_record() {
    let mut reader = FastaReader::new(Cursor::new(">seq1\r\nGADGSDHS\r\n\r\nFDSHJF\n>\r\nGATY"));

    let record1 = reader.next().unwrap().unwrap();
    assert_eq!(record1.name, "seq1");
    assert_eq!(record1.sequence, b"GADGSDHSFDSHJF");
    assert_eq!(reader.next().unwrap().unwrap_err().to_string(), "Missing FASTA header!");

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_sequence_second_record() {
    let mut reader = FastaReader::new(Cursor::new(">seq1\nGCAT\n>seq2\n"));

    let record1 = reader.next().unwrap().unwrap();
    assert_eq!(record1.name, "seq1");
    assert_eq!(record1.sequence, b"GCAT");

    assert_eq!(reader.next().unwrap().unwrap_err().to_string(), "Missing FASTA sequence!");

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn invalid_char_header_first_record() {
    let mut reader = FastaReader::new(Cursor::new(">seq1>seq2\r\nGADGSDHS\n>seq2\nGAT"));

    assert_eq!(
        reader.next().unwrap().unwrap_err().to_string(),
        "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a header!"
    );

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn invalid_char_sequence_first_record() {
    let mut reader = FastaReader::new(Cursor::new(">seq1\r\nGADGSDHS>CAT\n>seq2\nGAT"));

    let record1 = reader.next().unwrap().unwrap();
    assert_eq!(record1.name, "seq1");
    assert_eq!(record1.sequence, b"GADGSDHS");
    assert_eq!(
        reader.next().unwrap().unwrap_err().to_string(),
        "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a sequence!"
    );

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn invalid_char_header_second_record() {
    let mut reader = FastaReader::new(Cursor::new(">seq1\r\nGADGSDHS\n>seq2>seq3\nGAT"));

    let record1 = reader.next().unwrap().unwrap();
    assert_eq!(record1.name, "seq1");
    assert_eq!(record1.sequence, b"GADGSDHS");
    assert_eq!(
        reader.next().unwrap().unwrap_err().to_string(),
        "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a header!"
    );

    // Ensure iterator terminates
    assert!(reader.count() < 100);

    let mut reader = FastaReader::new(Cursor::new(">seq1\nGCAT\n>seq2>invalid"));

    let record1 = reader.next().unwrap().unwrap();
    assert_eq!(record1.name, "seq1");
    assert_eq!(record1.sequence, b"GCAT");

    assert_eq!(
        reader.next().unwrap().unwrap_err().to_string(),
        "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a header!"
    );

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn invalid_char_sequence_second_record() {
    let mut reader = FastaReader::new(Cursor::new(">seq1\r\nGADGSDHS\n>seq2\nGAT>CAT"));

    let Some(Ok(FastaSeq { name, sequence })) = reader.next() else {
        panic!("The first record should process without error.");
    };
    assert_eq!(name, "seq1");
    assert_eq!(sequence, b"GADGSDHS");

    let Some(Ok(FastaSeq { name, sequence })) = reader.next() else {
        panic!("The first record should process without error.");
    };
    assert_eq!(name, "seq2");
    assert_eq!(sequence, b"GAT");

    assert_eq!(
        reader.next().unwrap().unwrap_err().to_string(),
        "FASTA records must start with the '>' symbol on a newline, and no other '>' symbols can occur in a sequence!"
    );

    // Ensure iterator terminates
    assert!(reader.count() < 100);
}
