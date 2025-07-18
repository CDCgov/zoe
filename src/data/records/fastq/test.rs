use crate::prelude::*;
use std::io::Cursor;

#[test]
fn empty_file() {
    let mut reader = FastQReader::new(Cursor::new(""));
    assert!(reader.next().is_none());

    assert_eq!(
        FastQReader::from_readable(Cursor::new("")).err().unwrap().to_string(),
        "No FASTQ data was found!"
    );
}

#[test]
fn whitespace_only() {
    let mut reader = FastQReader::new(Cursor::new(" "));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '@' symbol at header line beginning! Ensure that the FASTQ file is not multi-line."
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn missing_at_sign_first_record() {
    let mut reader = FastQReader::new(Cursor::new("a"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '@' symbol at header line beginning! Ensure that the FASTQ file is not multi-line."
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_header_first_record() {
    let mut reader = FastQReader::new(Cursor::new("@\nATGC+\nIIII"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ header!");
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_sequence_first_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\n\n+\nIIII"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ sequence! See header: @seq1");
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn missing_plus_line_first_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n@seq2"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: @seq1"
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);

    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: @seq1"
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_quality_first_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\n"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ quality scores! See header: @seq1");
    // Ensure iterator terminates
    assert!(reader.count() < 100);

    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ quality scores! See header: @seq1");
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn mismatch_lengths_first_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIII"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Sequence and quality score length mismatch (4 ≠ 3)! See: @seq1"
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn missing_at_sign_second_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIIII\na"));
    let Some(Ok(FastQ {
        header,
        sequence,
        quality,
    })) = reader.next()
    else {
        panic!("Should parse correctly")
    };
    assert_eq!(header, "@seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '@' symbol at header line beginning! Ensure that the FASTQ file is not multi-line."
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_header_second_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIIII\n@\nATGC+\nIIII"));
    let Some(Ok(FastQ {
        header,
        sequence,
        quality,
    })) = reader.next()
    else {
        panic!("Should parse correctly")
    };
    assert_eq!(header, "@seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ header!");
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_sequence_second_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIIII\n@seq2\n\n+\nIIII"));
    let Some(Ok(FastQ {
        header,
        sequence,
        quality,
    })) = reader.next()
    else {
        panic!("Should parse correctly")
    };
    assert_eq!(header, "@seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ sequence! See header: @seq2");
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn missing_plus_line_second_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIIII\n@seq2\nATGC\n@seq3"));
    let Some(Ok(FastQ {
        header,
        sequence,
        quality,
    })) = reader.next()
    else {
        panic!("Should parse correctly")
    };
    assert_eq!(header, "@seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: @seq2"
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);

    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIIII\n@seq2\nATGC"));
    let Some(Ok(FastQ {
        header,
        sequence,
        quality,
    })) = reader.next()
    else {
        panic!("Should parse correctly")
    };
    assert_eq!(header, "@seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: @seq2"
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn empty_quality_second_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIIII\n@seq2\nATGC\n+\n"));
    let Some(Ok(FastQ {
        header,
        sequence,
        quality,
    })) = reader.next()
    else {
        panic!("Should parse correctly")
    };
    assert_eq!(header, "@seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ quality scores! See header: @seq2");
    // Ensure iterator terminates
    assert!(reader.count() < 100);

    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIIII\n@seq2\nATGC\n+"));
    let Some(Ok(FastQ {
        header,
        sequence,
        quality,
    })) = reader.next()
    else {
        panic!("Should parse correctly")
    };
    assert_eq!(header, "@seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ quality scores! See header: @seq2");
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn mismatch_lengths_second_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIIII\n@seq2\nATGC\n+\nIII"));
    let Some(Ok(FastQ {
        header,
        sequence,
        quality,
    })) = reader.next()
    else {
        panic!("Should parse correctly")
    };
    assert_eq!(header, "@seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Sequence and quality score length mismatch (4 ≠ 3)! See: @seq2"
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}
