use crate::prelude::*;
use std::io::{BufReader, Cursor, ErrorKind, Read};

struct InterruptedOnce<R> {
    inner:       R,
    interrupted: bool,
}

impl<R: Read> Read for InterruptedOnce<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.interrupted {
            self.inner.read(buf)
        } else {
            self.interrupted = true;
            Err(ErrorKind::Interrupted.into())
        }
    }
}

#[test]
fn interrupted_header_read_is_retried() {
    let input = InterruptedOnce {
        inner:       Cursor::new(b"@seq1\nA\n+\n!\n"),
        interrupted: false,
    };
    let mut reader = FastQReader::new(input);

    let record = reader.next().unwrap().unwrap();
    assert_eq!(record.header, "seq1");
    assert_eq!(record.sequence, Nucleotides::from_vec_unchecked(b"A".into()));
    assert_eq!(record.quality, QualityScores::try_from(b"!".to_vec()).unwrap());
    assert!(reader.next().is_none());
}

#[test]
fn crlf_with_one_byte_input_buffer() {
    let input = b"@seq1\r\nATGC\r\n+seq1\r\nIIII\r\n";
    let inner = BufReader::with_capacity(1, Cursor::new(input));
    let mut reader = FastQReader::from_bufreader(inner).unwrap();

    let record = reader.next().unwrap().unwrap();
    assert_eq!(record.header, "seq1");
    assert_eq!(record.sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(record.quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());
    assert!(reader.next().is_none());
}

#[test]
fn malformed_header_line_is_consumed() {
    let input = b"not-a-header\n@seq2\nA\n+\n!\n";
    let inner = BufReader::with_capacity(1, Cursor::new(input));
    let mut reader = FastQReader::from_bufreader(inner).unwrap();

    let error = reader.next().unwrap().unwrap_err();
    assert_eq!(error.kind(), ErrorKind::InvalidData);
    assert_eq!(
        error.to_string(),
        "Missing '@' symbol at header line beginning! Ensure that the FASTQ file is not multi-line."
    );

    let record = reader.next().unwrap().unwrap();
    assert_eq!(record.header, "seq2");
    assert_eq!(record.sequence, Nucleotides::from_vec_unchecked(b"A".into()));
    assert_eq!(record.quality, QualityScores::try_from(b"!".to_vec()).unwrap());
    assert!(reader.next().is_none());
}

#[test]
fn invalid_header_utf8_retains_details() {
    let input = b"@\xff\nA\n+\n!\n";
    let mut reader = FastQReader::new(Cursor::new(input));
    let error = reader.next().unwrap().unwrap_err();
    let expected = String::from_utf8(vec![0xff]).unwrap_err();

    assert_eq!(error.kind(), ErrorKind::InvalidData);
    assert_eq!(error.to_string(), expected.to_string());
}

#[test]
fn invalid_quality_is_rejected() {
    let input = b"@seq1\nA\n+\n \n";
    let mut reader = FastQReader::new(Cursor::new(input));
    let error = reader.next().unwrap().unwrap_err();

    assert_eq!(error.kind(), ErrorKind::InvalidData);
    assert_eq!(error.to_string(), "Quality scores contain invalid state!");
}

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
    assert_eq!(e.to_string(), "Missing FASTQ sequence! See header: seq1");
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
        "Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: seq1"
    );
    // Ensure iterator terminates
    assert!(reader.count() < 100);

    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: seq1"
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
    assert_eq!(e.to_string(), "Missing FASTQ quality scores! See header: seq1");
    // Ensure iterator terminates
    assert!(reader.count() < 100);

    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ quality scores! See header: seq1");
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}

#[test]
fn mismatch_lengths_first_record() {
    let mut reader = FastQReader::new(Cursor::new("@seq1\nATGC\n+\nIII"));
    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Sequence and quality score length mismatch (4 ≠ 3)! See: seq1");
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
    assert_eq!(header, "seq1");
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
    assert_eq!(header, "seq1");
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
    assert_eq!(header, "seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ sequence! See header: seq2");
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
    assert_eq!(header, "seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: seq2"
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
    assert_eq!(header, "seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(
        e.to_string(),
        "Missing '+' line! Ensure that the FASTQ file is not multi-line. See header: seq2"
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
    assert_eq!(header, "seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ quality scores! See header: seq2");
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
    assert_eq!(header, "seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Missing FASTQ quality scores! See header: seq2");
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
    assert_eq!(header, "seq1");
    assert_eq!(sequence, Nucleotides::from_vec_unchecked(b"ATGC".into()));
    assert_eq!(quality, QualityScores::try_from(b"IIII".to_vec()).unwrap());

    let Some(Err(e)) = reader.next() else {
        panic!("Should throw error")
    };
    assert_eq!(e.to_string(), "Sequence and quality score length mismatch (4 ≠ 3)! See: seq2");
    // Ensure iterator terminates
    assert!(reader.count() < 100);
}
