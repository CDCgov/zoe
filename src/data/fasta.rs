use super::types::nucleotides::*;
use std::io::BufRead;
use std::fs::File;

#[derive(Debug)]
pub struct FastaSeq {
    pub name: Vec<u8>,
    pub sequence: Vec<u8>,
}

pub struct FastaReader<R: std::io::Read> {
    pub reader: std::io::BufReader<R>,
    pub r_buff: Vec<u8>,
}
pub struct FastaNT {
    pub name: String,
    pub sequence: Nucleotides,
}

// FastaAA & AminoAcids

use std::fmt;
impl fmt::Display for FastaSeq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,">{}\n{}\n", String::from_utf8_lossy(&self.name), String::from_utf8_lossy(&self.sequence))
    }
}

impl FastaSeq {
    pub fn reverse_complement(&mut self) {
        use super::matrices::REV_COMP;
        let sequence: Vec<u8> = self.sequence
            .iter()
            .rev()
            .copied()
            .map(|x| REV_COMP[x as usize])
            .collect();
        self.sequence = sequence;
    }
}

impl<R: std::io::Read> FastaReader<R> {
    pub fn new(inner: R) -> Self {
        FastaReader {
            reader: std::io::BufReader::new(inner),
            r_buff: Vec::new(),
        }
    }
}

impl FastaReader<std::fs::File> {
    pub fn from_filename(filename: &str) -> Result<FastaReader<File>,String> {
        match File::open(filename) {
            Err(why) => Err(format!("Couldn't open fasta file '{filename}': {why}")),
            Ok(file) =>  Ok(FastaReader::new(file)),
        }       
    }
}

impl<R: std::io::Read> Iterator for FastaReader<R> {
    type Item = FastaSeq;

    fn next(&mut self) -> Option<Self::Item> {  
       self.r_buff.clear();

       let bytes = match self.reader.read_until(b'>', &mut self.r_buff) {
            Err(_)            => 0,
            Ok(bytes)   => bytes,
        };

        match bytes {
            0 => None,
            1 => self.next(),
            _ => {
                if self.r_buff.ends_with(b">") {
                    self.r_buff.pop();
                }

                let mut lines = self.r_buff.split(|x| *x == b'\n' || *x == b'\r');
                let name = match lines.next() {
                    Some(t) => t.to_vec(),
                    None => b"UNKNOWN".to_vec(),
                };
    
                let sequence: Vec<u8> = lines.flatten().copied().collect();
    
                if sequence.len() > 0 {
                    return Some(FastaSeq { name, sequence });
                } else {
                    return None
                }
            }
        }
    }
}