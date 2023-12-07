use crate::{
    composition::NucleotideCounts,
    data::{
        mappings::{
            GENETIC_CODE, IS_IUPAC_BASE, IS_UNALIGNED_IUPAC_BASE, TO_DNA_UC, TO_REVERSE_COMPLEMENT, TO_UNALIGNED_DNA_UC,
        },
        vec_types::ValidateSequence,
    },
    simd::{SimdByteFunctions, SimdMaskFunctions},
};
use std::simd::{LaneCount, SupportedLaneCount};
use std::slice::ChunksExact;

use super::{amino_acids::AminoAcids, Uint};

#[derive(Debug, Clone, Default, PartialEq)]
#[repr(transparent)]
pub struct Nucleotides(pub(crate) Vec<u8>);

impl Nucleotides {
    // Standard functions
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Nucleotides(Vec::new())
    }

    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_slice()
    }

    #[inline]
    #[must_use]
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        self.0.as_mut_slice()
    }

    #[inline]
    #[must_use]
    pub fn as_vec(&self) -> &Vec<u8> {
        &self.0
    }

    #[inline]
    #[must_use]
    pub fn as_mut_vec(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }

    #[inline]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    // Manipulation
    #[inline]
    pub fn find_and_replace(&mut self, needle: u8, replacement: u8) {
        crate::search::find_and_replace(&mut self.0, needle, replacement);
    }

    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
    }

    #[inline]
    #[must_use]
    pub fn gc_content(&self) -> usize {
        crate::composition::gc_content_simd::<32>(&self.0)
    }

    // Validation

    /// Only retains valid [IUPAC bases](https://www.bioinformatics.org/sms/iupac.html).
    #[inline]
    pub fn retain_iupac_bases(&mut self) {
        self.0.retain_by_validation(IS_IUPAC_BASE);
    }

    /// Only retains valid [IUPAC bases](https://www.bioinformatics.org/sms/iupac.html).
    #[inline]
    pub fn retain_unaligned_bases(&mut self) {
        self.0.retain_by_validation(IS_UNALIGNED_IUPAC_BASE);
    }

    /// Only retains valid [IUPAC](https://www.bioinformatics.org/sms/iupac.html) DNA and converts to uppercase.
    #[inline]
    pub fn retain_dna_uc(&mut self) {
        self.0.retain_by_recoding(TO_DNA_UC);
    }

    /// Only retains valid, unaligned [IUPAC](https://www.bioinformatics.org/sms/iupac.html) DNA and converts to uppercase.
    #[inline]
    pub fn retain_unaligned_dna_uc(&mut self) {
        self.0.retain_by_recoding(TO_UNALIGNED_DNA_UC);
    }

    // Domain functions
    #[inline]
    #[must_use]
    pub fn reverse_complement(&self) -> Self {
        Self(reverse_complement(&self.0))
    }

    #[inline]
    #[must_use]
    pub fn translate(&self) -> AminoAcids {
        AminoAcids(translate_sequence(&self.0))
    }

    #[must_use]
    pub fn into_aa_iter(&self) -> TranslatedNucleotidesIter {
        TranslatedNucleotidesIter {
            codons:        self.0.chunks_exact(3),
            has_remainder: true,
        }
    }

    #[must_use]
    pub fn into_base_counts<T: Uint>(&self) -> NucleotideCounts<T> {
        self.0.iter().fold(NucleotideCounts::new(), |acc, &b| acc + b)
    }

    /// # Distance
    ///
    /// Calculates hamming distance between `self`and another sequence.
    ///
    /// # Example
    /// ```
    /// use zoe::data::types::nucleotides::Nucleotides;
    ///
    /// let s1: Nucleotides = b"ATGCATCGATCGATCGATCGATCGATCGATGC".into();
    /// let s2: Nucleotides = b"ATGCATnGATCGATCGATCGAnCGATCGATnC".into();
    ///
    /// assert!(3 == s1.distance_hamming(&s2));
    ///
    /// let s3: &[u8] = b"ATGCATnGATCGATCGATCGAnCGATCGATnC";
    /// assert!(3 == s1.distance_hamming(&s3));
    /// ```
    ///
    #[inline]
    #[must_use]
    pub fn distance_hamming<T: AsRef<[u8]> + MaybeNucleic>(&self, other_sequence: &T) -> usize {
        crate::distance::hamming_simd::<16>(&self.0, other_sequence.as_ref())
    }

    // Associated functions

    /// Generate a random DNA sequence of given `length` and using a random `seed`.  Canonical DNA only contains A, C, G, or T.
    #[must_use]
    pub fn generate_random_dna(length: usize, seed: u64) -> Self {
        Nucleotides(crate::generate::rand_sequence(b"AGCT", length, seed))
    }
}

/// Marker trait to restrict usage to types that might possible contain nucleotides.
pub trait MaybeNucleic {}
impl MaybeNucleic for Nucleotides {}
impl MaybeNucleic for Vec<u8> {}
impl MaybeNucleic for &[u8] {}

#[inline]
#[must_use]
pub fn translate_sequence(s: &[u8]) -> Vec<u8> {
    let mut codons = s.chunks_exact(3);
    let mut aa_sequence = Vec::with_capacity(s.len() / 3 + 1);
    let mut codon_copy = [0; 3];

    for codon in codons.by_ref() {
        aa_sequence.push(if is_partial_codon(codon) {
            b'~'
        } else {
            codon_copy.copy_from_slice(codon);
            codon_copy.make_ascii_uppercase();
            *GENETIC_CODE.get(&codon_copy).unwrap_or(&b'X')
        });
    }

    let tail = codons.remainder();
    if is_partial_codon(tail) {
        aa_sequence.push(b'~');
    }

    aa_sequence
}

pub struct TranslatedNucleotidesIter<'a> {
    codons:        ChunksExact<'a, u8>,
    has_remainder: bool,
}

impl<'a> Iterator for TranslatedNucleotidesIter<'a> {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(codon) = self.codons.next() {
            if is_partial_codon(codon) {
                Some(b'~')
            } else {
                let mut codon_copy = [0; 3];
                codon_copy.copy_from_slice(codon);
                codon_copy.make_ascii_uppercase();
                Some(*GENETIC_CODE.get(&codon_copy).unwrap_or(&b'X'))
            }
        } else if self.has_remainder && is_partial_codon(self.codons.remainder()) {
            self.has_remainder = false;
            Some(b'~')
        } else {
            None
        }
    }
}

#[inline]
#[must_use]
fn is_partial_codon(codon: &[u8]) -> bool {
    if codon.is_empty() {
        false
    } else if codon.len() < 3 {
        true
    } else {
        let count: u8 = codon.iter().map(|b| u8::from(*b == b'-' || *b == b'~' || *b == b'.')).sum();
        count == 1 || count == 2
    }
}

#[inline]
#[must_use]
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .copied()
        .map(|x| TO_REVERSE_COMPLEMENT[x as usize])
        .collect()
}

/// Reverse complement of a nucleotide sequences using explicit SIMD instructions.
///
/// # Note
///
/// Works well for Haswell and above architectures. Recommend 32 lanes for x86.
/// Both `swizzle_dyn` and `gather_or` are too slow to be relevant versus the
/// scalar implementation.
#[inline]
#[must_use]
#[allow(dead_code)]
pub fn reverse_complement_simd<const N: usize>(bases: &[u8]) -> Vec<u8>
where
    LaneCount<N>: SupportedLaneCount, {
    let (pre, mid, sfx) = bases.as_simd::<N>();
    let mut reverse_complement: Vec<u8> = Vec::with_capacity(bases.len());

    reverse_complement.extend(sfx.iter().rev().copied().map(|x| TO_REVERSE_COMPLEMENT[x as usize]));

    reverse_complement.extend(
        mid.iter()
            .map(|&v| {
                let mut rev = v.reverse();
                let lowercase = rev.is_ascii_lowercase();
                rev = lowercase.make_selected_ascii_uppercase(&rev);

                rev.swap_byte_pairs(b'T', b'A');
                rev.swap_byte_pairs(b'G', b'C');
                rev.swap_byte_pairs(b'R', b'Y');
                rev.swap_byte_pairs(b'K', b'M');
                rev.swap_byte_pairs(b'B', b'V');
                rev.swap_byte_pairs(b'H', b'D');
                rev.if_value_then_replace(b'U', b'A');

                rev = lowercase.make_selected_ascii_lowercase(&rev);
                rev.to_array()
            })
            .rev()
            .flatten(),
    );

    reverse_complement.extend(pre.iter().rev().copied().map(|x| TO_REVERSE_COMPLEMENT[x as usize]));

    reverse_complement
}

#[cfg(test)]
mod bench;
#[cfg(test)]
mod test;

mod std_traits;

#[allow(unused_imports)]
pub use std_traits::*;
