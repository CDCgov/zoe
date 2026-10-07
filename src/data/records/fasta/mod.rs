use crate::data::{
    id_types::FastaIDs,
    types::{
        amino_acids::AminoAcids,
        nucleotides::{self, Nucleotides, ToDNA, Translate},
    },
    validation::SimdDisplay,
};

mod reader;

pub use reader::*;

#[cfg(feature = "dev-generic-fasta")]
pub mod generic;

#[cfg(test)]
mod test;

/// Provides a container struct for data from a generic
/// [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaSeq {
    pub name:     String,
    pub sequence: Vec<u8>,
}

/// Similar to [`FastaSeq`] but assumes that the `sequence` contains valid
/// [`Nucleotides`].
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaNT {
    pub name:     String,
    pub sequence: Nucleotides,
}

/// Similar to [`FastaSeq`] but assumes that the `sequence` contains valid
/// [`AminoAcids`].
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaAA {
    pub name:     String,
    pub sequence: AminoAcids,
}

impl FastaSeq {
    /// Reverse complements the sequence stored in the struct using a new
    /// buffer.
    pub fn reverse_complement(&mut self) {
        self.sequence = nucleotides::reverse_complement(&self.sequence);
    }

    /// Recodes to uppercase IUPAC DNA with corrected gaps, otherwise
    /// mapping to `N`. Returns [`FastaNT`].
    #[inline]
    #[must_use]
    pub fn recode_to_dna(self) -> FastaNT {
        FastaNT {
            name:     self.name,
            sequence: self.sequence.recode_to_dna(),
        }
    }

    /// Filters and recodes to uppercase IUPAC DNA with corrected gaps. Returns
    /// [`FastaNT`].
    #[inline]
    #[must_use]
    pub fn filter_to_dna(self) -> FastaNT {
        FastaNT {
            name:     self.name,
            sequence: self.sequence.filter_to_dna(),
        }
    }

    /// Filters and recodes to uppercase IUPAC DNA without gaps. Returns
    /// [`FastaNT`].
    #[inline]
    #[must_use]
    pub fn filter_to_dna_unaligned(self) -> FastaNT {
        FastaNT {
            name:     self.name,
            sequence: self.sequence.filter_to_dna_unaligned(),
        }
    }

    /// For an annotated `FASTA` with format `id{annotation}` returns a tuple
    /// of the id and annotated taxon.
    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, &str)> {
        self.name.get_id_taxon()
    }

    /// Translates the stored [`Vec<u8>`] to [`AminoAcids`] using a new buffer.
    ///
    /// See [`translate_sequence`] for more details.
    ///
    /// [`translate_sequence`]: nucleotides::translate_sequence
    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     self.name,
            sequence: AminoAcids(nucleotides::translate_sequence(&self.sequence)),
        }
    }
}

impl FastaNT {
    /// Reverse complements the sequence stored in the struct using a new buffer.
    #[inline]
    pub fn reverse_complement(&mut self) {
        self.sequence.make_reverse_complement();
    }

    /// Translates the stored [`Nucleotides`] to [`AminoAcids`] using a new buffer.
    #[must_use]
    pub fn translate(self) -> FastaAA {
        FastaAA {
            name:     self.name,
            sequence: self.sequence.translate(),
        }
    }

    /// For an annotated FASTA file with format `id{annotation}` returns a tuple
    /// of the id and annotated taxon.
    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, &str)> {
        self.name.get_id_taxon()
    }
}

/// Allows converting from [`FastaSeq`] to [`FastaNT`] without checks or
/// filtering.
impl From<FastaSeq> for FastaNT {
    fn from(record: FastaSeq) -> Self {
        FastaNT {
            name:     record.name,
            sequence: record.sequence.into(),
        }
    }
}

impl FastaAA {
    /// For an annotated `FASTA` with format `id{annotation}` returns a tuple
    /// of the id and annotated taxon.
    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, &str)> {
        self.name.get_id_taxon()
    }
}

impl std::fmt::Display for FastaSeq {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, self.sequence.display_ascii_or_lossy())
    }
}

impl std::fmt::Display for FastaNT {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, self.sequence)
    }
}

impl std::fmt::Display for FastaAA {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, ">{}\n{}\n", self.name, self.sequence)
    }
}
