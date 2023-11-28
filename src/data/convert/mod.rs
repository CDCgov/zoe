use super::types::nucleotides::Nucleotides;

pub trait ToDNA {
    fn filter_to_dna(&self) -> Nucleotides;
}

impl ToDNA for String {
    fn filter_to_dna(&self) -> Nucleotides {
        let mut n = Nucleotides(self.as_bytes().to_vec());
        n.retain_dna_uc();
        n
    }
}

impl ToDNA for Vec<u8> {
    fn filter_to_dna(&self) -> Nucleotides {
        let mut n = Nucleotides(self.clone());
        n.retain_dna_uc();
        n
    }
}

impl ToDNA for &[u8] {
    fn filter_to_dna(&self) -> Nucleotides {
        let mut n = Nucleotides(self.to_vec());
        n.retain_dna_uc();
        n
    }
}
