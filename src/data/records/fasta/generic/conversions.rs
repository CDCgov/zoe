use crate::data::{
    amino_acids::AminoAcids,
    fasta::{FastaAA, FastaNT, FastaSeq, generic::Fasta},
    nucleotides::Nucleotides,
};

impl From<FastaSeq> for Fasta {
    fn from(value: FastaSeq) -> Self {
        Fasta {
            header:   value.name,
            sequence: value.sequence,
        }
    }
}

impl From<FastaNT> for Fasta<Nucleotides> {
    fn from(value: FastaNT) -> Self {
        Fasta {
            header:   value.name,
            sequence: value.sequence,
        }
    }
}

impl From<FastaAA> for Fasta<AminoAcids> {
    fn from(value: FastaAA) -> Self {
        Fasta {
            header:   value.name,
            sequence: value.sequence,
        }
    }
}

impl From<Fasta> for FastaSeq {
    fn from(value: Fasta) -> Self {
        FastaSeq {
            name:     value.header,
            sequence: value.sequence,
        }
    }
}

impl From<Fasta<Nucleotides>> for FastaNT {
    fn from(value: Fasta<Nucleotides>) -> Self {
        FastaNT {
            name:     value.header,
            sequence: value.sequence,
        }
    }
}

impl From<Fasta<AminoAcids>> for FastaAA {
    fn from(value: Fasta<AminoAcids>) -> Self {
        FastaAA {
            name:     value.header,
            sequence: value.sequence,
        }
    }
}
