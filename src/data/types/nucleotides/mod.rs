use crate::{
    alignment::DNAProfileIndices,
    composition::NucleotideCounts,
    data::{
        mappings::{
            CodonConvert, ANY_TO_DNA_CANONICAL_UPPER, GENETIC_CODE, IS_IUPAC_BASE, IS_UNALIGNED_IUPAC_BASE,
            IUPAC_TO_DNA_CANONICAL, IUPAC_TO_DNA_CANONICAL_UPPER, TO_DNA_UC, TO_REVERSE_COMPLEMENT, TO_UNALIGNED_DNA_UC,
        },
        types::{amino_acids::AminoAcids, Uint},
        vec_types::ValidateSequence,
    },
    simd::{SimdByteFunctions, SimdMaskFunctions},
};
use std::{
    ops::Range,
    simd::{LaneCount, SupportedLaneCount},
    slice::ChunksExact,
};

/// [`Nucleotides`] is a transparent, new-type wrapper around [`Vec<u8>`]
/// that provides DNA-specific functionality and semantics. It may contain either
/// aligned or unaligned valid IUPAC letters.
///
/// *NB: No type-state guarantees have been decided on at this time.*
#[derive(Debug, Clone, Default, PartialEq)]
#[repr(transparent)]
pub struct Nucleotides(pub(crate) Vec<u8>);

impl Nucleotides {
    // std

    /// Create a new `Nucleotides` empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Nucleotides(Vec::new())
    }

    /// The length of the stored sequence.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Is the sequence empty?
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Obtains the bytes as a slice.
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
    /// Create an iterator over the nucleotides as `u8`.
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }

    /// Gets the base at the zero-based index, returning an `Option`.
    #[inline]
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: std::slice::SliceIndex<[u8]>, {
        self.0.get(index)
    }

    /// The sequence is re-encoded as `DNAProfileIndices` in-place.
    #[must_use]
    pub fn into_dna_profile_indices(mut self) -> DNAProfileIndices {
        for base in &mut self.0 {
            *base = crate::data::mappings::TO_DNA_PROFILE_INDEX[*base as usize];
        }
        DNAProfileIndices(self.0)
    }

    /// Replaces the provided [`Range`] with the specified byte. If the range
    /// does not exist, the function does nothing.
    #[inline]
    pub fn mask_if_exists(&mut self, range: Range<usize>, replacement: u8) {
        if let Some(slice) = self.0.get_mut(range) {
            for b in slice {
                *b = replacement;
            }
        }
    }

    /// Truncates the length of the sequence to the specified `new_length`. This
    /// is equivalent to 3' trimming up to and including the index.
    #[inline]
    pub fn shorten_to(&mut self, new_length: usize) {
        self.0.truncate(new_length);
    }

    /// Cuts the 5' end of the [`Nucleotides`] just prior to the new starting
    /// index (0-based). Be aware that this clones the internal buffer!
    #[inline]
    pub fn cut_to_start(&mut self, new_start: usize) {
        *self = Nucleotides(self.0.drain(new_start..).collect());
    }

    /// Provides the count of G and C bases.
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

    /// Recodes the stored sequence to an uppercase canonical (ACTG + N) one.
    /// Any non-canonical base becomes N
    #[inline]
    pub fn recode_any_to_actgn_uc(&mut self) {
        self.0.recode(ANY_TO_DNA_CANONICAL_UPPER);
    }

    /// Recodes the stored sequence of valid IUPAC codes to a canonical (ACTG + N)
    /// sequence. Ambiguous bases become N while non-IUPAC bytes are left
    /// unchanged..
    #[inline]
    pub fn recode_iupac_to_actgn(&mut self) {
        self.0.recode(IUPAC_TO_DNA_CANONICAL);
    }

    /// Recodes the stored sequence of valid IUPAC codes to an uppercase
    /// canonical (ACTG + N) sequence. Ambiguous bases become N while non-IUPAC
    /// bytes are left unchanged.
    #[inline]
    pub fn recode_iupac_to_actgn_uc(&mut self) {
        self.0.recode(IUPAC_TO_DNA_CANONICAL_UPPER);
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

    /// Translate the stored [Nucleotides] to [`AminoAcids`].
    #[inline]
    #[must_use]
    pub fn translate(&self) -> AminoAcids {
        AminoAcids(translate_sequence(&self.0))
    }

    /// Creates an iterator for [`AminoAcids`] translation.
    #[inline]
    #[must_use]
    pub fn into_aa_iter(&self) -> TranslatedNucleotidesIter {
        TranslatedNucleotidesIter {
            codons:        self.0.chunks_exact(3),
            has_remainder: true,
        }
    }

    /// Creates [`NucleotideCounts`] (ACGT + N + -) statistics using the specified `const` [Uint].
    #[must_use]
    pub fn into_base_counts<T: Uint>(&self) -> NucleotideCounts<T> {
        self.0.iter().fold(NucleotideCounts::new(), |acc, &b| acc + b)
    }

    /// # Distance
    ///
    /// Calculates hamming distance between [self] and another sequence.
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

    /// Generate a random DNA sequence of given `length` and using a random
    /// `seed`.  Canonical DNA only contains A, C, G, or T.
    #[must_use]
    pub fn generate_random_dna(length: usize, seed: u64) -> Self {
        Nucleotides(crate::generate::rand_sequence(b"AGCT", length, seed))
    }
}

/// Marker trait to restrict usage to types that might possible contain
/// nucleotides.
pub trait MaybeNucleic {}
impl MaybeNucleic for Nucleotides {}
impl MaybeNucleic for Vec<u8> {}
impl MaybeNucleic for &[u8] {}

/// Translates a byte slice into an amino acid byte vector.
#[inline]
#[must_use]
pub fn translate_sequence(s: &[u8]) -> Vec<u8> {
    let mut codons = s.chunks_exact(3);
    let mut aa_sequence = Vec::with_capacity(s.len() / 3 + 1);

    for codon in codons.by_ref() {
        aa_sequence.push(if is_partial_codon(codon) {
            b'~'
        } else {
            *GENETIC_CODE.get(&codon.to_upper_u32()).unwrap_or(&b'X')
        });
    }

    let tail = codons.remainder();
    if is_partial_codon(tail) {
        aa_sequence.push(b'~');
    }

    aa_sequence
}

/// Iterator for translating [Nucleotides] into [`AminoAcids`].
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
                Some(*GENETIC_CODE.get(&codon.to_upper_u32()).unwrap_or(&b'X'))
            }
        } else if self.has_remainder && is_partial_codon(self.codons.remainder()) {
            self.has_remainder = false;
            Some(b'~')
        } else {
            None
        }
    }
}

/// A codon is considered *partial* if it has fewer than 3 IUPAC bases.
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

/// Performs the DNA reverse complement of the byte slice into a new vector.
/// Assumes ASCII input.
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
/// Both [`swizzle_dyn`](std::simd::prelude::Simd::swizzle_dyn) and
/// [`gather_or`](std::simd::prelude::Simd::gather_or) are too slow to be
/// relevant versus the scalar implementation.
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

                rev.exchange_byte_pairs(b'T', b'A');
                rev.exchange_byte_pairs(b'G', b'C');
                rev.exchange_byte_pairs(b'R', b'Y');
                rev.exchange_byte_pairs(b'K', b'M');
                rev.exchange_byte_pairs(b'B', b'V');
                rev.exchange_byte_pairs(b'H', b'D');
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
