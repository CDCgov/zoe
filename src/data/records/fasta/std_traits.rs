use crate::data::{
    fasta::{FastaAA, FastaNT, FastaSeq},
    validation::SimdDisplay,
};
use std::fmt::Display;

impl Display for FastaSeq {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, self.sequence.display_ascii_or_lossy())
    }
}

impl Display for FastaNT {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, self.sequence)
    }
}

impl Display for FastaAA {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, self.sequence)
    }
}
