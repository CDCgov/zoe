#![allow(clippy::doc_markdown)]
use crate::{
    DEFAULT_SIMD_LANES, data::types::nucleotides::NucleotidesReadable, distance::hamming_simd, math::MapFloat,
    private::Sealed, simd::SimdByteFunctions,
};
use std::simd::prelude::*;

#[cfg(test)]
mod bench;
/// Helper functions for pairwise frequency tabulation.
mod tabulation;

#[cfg(test)]
pub(crate) mod test;

pub use tabulation::dna_substitution_matrix;
pub(crate) use tabulation::{hamming_dist_from_sub_matrix, total_and_frequencies};

/// Calculates the p-distance (proportion of differing sites) between two
/// nucleotide sequences. The p-distance considers only valid, canonical
/// nucleotide pairs (A, C, G, T/U). Casing is ignored and only shared sequence
/// length--the smaller of the two--are compared. The algorithm uses SIMD
/// operations for improved performance.
///
/// ## Returns
///
/// - `Some(f64)` - The p-distance between sequences if valid positions are found.
/// - `None` - If no valid positions are found to compare.
///
/// ## Type Parameters
///
/// - `N` - SIMD lane count, must be a supported lane count
///
/// ## Example
///
/// ```
/// # use zoe::distance::dna::p_distance_acgt;
///
/// let seq1 = b"ATGCATGC";
/// let seq2 = b"atgtttnc";
/// let distance = p_distance_acgt::<16>(seq1, seq2);
/// assert!((distance.unwrap() - 0.2857142).abs() < 0.1);
/// ```
#[must_use]
#[allow(clippy::cast_precision_loss)]
#[cfg_attr(feature = "multiversion", multiversion::multiversion(targets = "simd"))]
pub fn p_distance_acgt<const N: usize>(x: &[u8], y: &[u8]) -> Option<f64> {
    let alpha_v = [Simd::splat(b'A'), Simd::splat(b'C'), Simd::splat(b'G'), Simd::splat(b'T')];

    let (x, y) = if x.len() < y.len() {
        (x, &y[..x.len()])
    } else {
        (&x[..y.len()], y)
    };

    let mut mismatches: usize = 0;
    let mut valid_length = 0;

    let mut x = x.chunks_exact(N * 255);
    let mut y = y.chunks_exact(N * 255);

    for (c1, c2) in x.by_ref().zip(y.by_ref()) {
        let mut accum: Simd<u8, N> = Simd::splat(0);

        let mut c1 = c1.chunks_exact(N);
        let mut c2 = c2.chunks_exact(N);

        for (v1, v2) in c1.by_ref().zip(c2.by_ref()) {
            let mut v1: Simd<u8, N> = Simd::from_slice(v1).to_ascii_uppercase();
            let mut v2: Simd<u8, N> = Simd::from_slice(v2).to_ascii_uppercase();
            v1.if_value_then_replace(b'U', b'T');
            v2.if_value_then_replace(b'U', b'T');

            let mut valid = Mask::from_array([false; N]);
            for a in alpha_v {
                valid |= a.simd_eq(v1);
                valid |= a.simd_eq(v2);
            }

            valid_length += valid.to_bitmask().count_ones() as usize;

            let m = (v1.simd_ne(v2) & valid).to_simd();
            // True => -1, so - -1 => +1
            accum -= m.cast();
        }

        let accum2: Simd<u16, N> = accum.cast();
        mismatches += accum2.reduce_sum() as usize;
    }

    let x = x.remainder();
    let y = y.remainder();
    let mut accum: Simd<u8, N> = Simd::splat(0);
    let mut c1 = x.chunks_exact(N);
    let mut c2 = y.chunks_exact(N);

    for (v1, v2) in c1.by_ref().zip(c2.by_ref()) {
        let mut v1: Simd<u8, N> = Simd::from_slice(v1).to_ascii_uppercase();
        let mut v2: Simd<u8, N> = Simd::from_slice(v2).to_ascii_uppercase();
        v1.if_value_then_replace(b'U', b'T');
        v2.if_value_then_replace(b'U', b'T');
        let mut valid1 = Mask::from_array([false; N]);
        let mut valid2 = Mask::from_array([false; N]);
        for a in alpha_v {
            valid1 |= a.simd_eq(v1);
            valid2 |= a.simd_eq(v2);
        }
        let valid = valid1 & valid2;

        valid_length += valid.to_bitmask().count_ones() as usize;

        let m = (v1.simd_ne(v2) & valid).to_simd();
        // True => -1, so - -1 => +1
        accum -= m.cast();
    }
    let accum2: Simd<u16, N> = accum.cast();
    mismatches += accum2.reduce_sum() as usize;

    let r1 = c1.remainder();
    let r2 = c2.remainder();

    mismatches += r1
        .iter()
        .zip(r2)
        .filter(|(a, b)| {
            let mut a = a.to_ascii_uppercase();
            let mut b = b.to_ascii_uppercase();
            if a == b'U' {
                a = b'T';
            }

            if b == b'U' {
                b = b'T';
            }
            let valid = matches!(a, b'A' | b'G' | b'T' | b'C') && matches!(b, b'A' | b'G' | b'T' | b'C');
            valid_length += usize::from(valid);
            a != b && valid
        })
        .count();

    if valid_length > 0 {
        Some(mismatches as f64 / valid_length as f64)
    } else {
        None
    }
}

/// An extension trait to add DNA distance methods to data that looks like [`Nucleotides`].
///
/// [`Nucleotides`]: crate::data::types::nucleotides::Nucleotides
pub trait NucleotidesDistance: NucleotidesReadable + Sealed {
    /// Calculates hamming distance between [`self`] and another sequence.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     data::types::nucleotides::{Nucleotides, NucleotidesView},
    /// #     distance::dna::NucleotidesDistance
    /// # };
    ///
    /// let s1: Nucleotides = b"ATGCATCGATCGATCGATCGATCGATCGATGC".into();
    /// let s2: Nucleotides = b"ATGCATnGATCGATCGATCGAnCGATCGATnC".into();
    ///
    /// assert!(3 == s1.distance_hamming(&s2));
    ///
    /// let s3 = NucleotidesView::from_bytes_unchecked(b"ATGCATnGATCGATCGATCGAnCGATCGATnC");
    /// assert!(3 == s1.distance_hamming(&s3));
    /// ```
    #[inline]
    #[must_use]
    fn distance_hamming<T: NucleotidesReadable>(&self, other_sequence: &T) -> usize {
        hamming_simd::<16>(self.nucleotide_bytes(), other_sequence.nucleotide_bytes())
    }

    /// Computes the JC69 distance between [`self`] and another sequence. See
    /// [`jukes_cantor_69`] for more details.
    #[inline]
    #[must_use]
    fn distance_jc69<T: NucleotidesReadable>(&self, other_sequence: &T) -> Option<f64> {
        jukes_cantor_69(self.nucleotide_bytes(), other_sequence.nucleotide_bytes())
    }

    /// Computes the K80 distance between [`self`] and another sequence. See
    /// [`kimura_80`] for more details.
    #[inline]
    #[must_use]
    fn distance_k80<T: NucleotidesReadable>(&self, other_sequence: &T) -> Option<f64> {
        kimura_80(self.nucleotide_bytes(), other_sequence.nucleotide_bytes())
    }

    /// Computes the K81 distance between [`self`] and another sequence. See
    /// [`kimura_81`] for more details.
    #[inline]
    #[must_use]
    fn distance_k81<T: NucleotidesReadable>(&self, other_sequence: &T) -> Option<f64> {
        kimura_81(self.nucleotide_bytes(), other_sequence.nucleotide_bytes())
    }

    /// Computes the F81 distance between [`self`] and another sequence. See
    /// [`felsenstein_81`] for more details.
    #[inline]
    #[must_use]
    fn distance_f81<T: NucleotidesReadable>(&self, other_sequence: &T) -> Option<f64> {
        felsenstein_81(self.nucleotide_bytes(), other_sequence.nucleotide_bytes())
    }

    /// Computes the TN93 distance between [`self`] and another sequence. See
    /// [`tamura_nei_93`] for more details.
    #[inline]
    #[must_use]
    fn distance_tn93<T: NucleotidesReadable>(&self, other_sequence: &T) -> Option<f64> {
        tamura_nei_93(self.nucleotide_bytes(), other_sequence.nucleotide_bytes())
    }

    #[inline]
    #[must_use]
    fn to_dna_substitution_matrix<T: NucleotidesReadable>(&self, other_sequence: &T) -> [[u32; 4]; 4] {
        dna_substitution_matrix(self.nucleotide_bytes(), other_sequence.nucleotide_bytes())
    }
}

impl<T: NucleotidesReadable + Sealed> NucleotidesDistance for T {}

/// ## Jukes-Cantor (JC-69) nucleotide substitution model.
///
/// This model estimates the number of substitutions per site between two sequences
/// under the assumption of equal base frequencies and equal substitution rates among
/// all possible nucleotide changes.
///
/// See [Assumptions](super::dna)
///
/// The formula used is:
///
/// $$ d = - \frac{3}{4} \ln \left(1 - \frac{4}{3} p \right) $$
///
/// where $d$ is the JC-69 evolutionary distance and $p$ is the proportion of differing sites
/// between the two sequences. (p-distance)
///
/// ## Citations
///
/// - Jukes, T., and Cantor, C. (1969). "Evolution of Protein Molecules."
///   Mammalian Protein Metabolism, New York: Academic Press, III(3), 21–132.
#[inline]
#[must_use]
pub fn jukes_cantor_69(seq1: &[u8], seq2: &[u8]) -> Option<f64> {
    if seq1.is_empty() || seq2.is_empty() {
        return None;
    }

    // Specialize the JC
    if let Some(p) = p_distance_acgt::<{ DEFAULT_SIMD_LANES }>(seq1, seq2) {
        (-0.75 * (1.0 - p * 4.0 / 3.0).ln()).into_option()
    } else {
        None
    }
}

/// ## Kimura 2-Parameter (K-80) nucleotide substitution model.
///
/// This model accounts for different rates of **transitions** (A ↔ G or C ↔ T) and
/// **transversions** (A ↔ C, A ↔ T, C ↔ G, G ↔ T) between sequences.
///
/// See [Assumptions](super::dna)
///
/// The formula used is:
///
/// $$ d = -\frac{1}{2} \ln\left(1 - 2p - q\right) - \frac{1}{4} \ln\left(1 - 2q\right) $$
///
/// where $d$ is the evolutionary distance, $p$ is the proportion of transitions, $\left( p = \frac{\text{transitions}}{\text{total bases}}\right)$
/// and $q$ is the proportion of transversions. $\left( q = \frac{\text{transversions}}{\text{total bases}}\right)$
///
///
/// ## Citations
///
/// - Kimura, M. (1980). "A simple method for estimating evolutionary rates of
///   base substitutions through comparative studies of nucleotide sequences."
///   Journal of Molecular Evolution. 16, 111-120.
#[inline]
#[must_use]
pub fn kimura_80(seq1: &[u8], seq2: &[u8]) -> Option<f64> {
    if seq1.is_empty() || seq2.is_empty() {
        return None;
    }
    let sub_matrix = dna_substitution_matrix(seq1, seq2);
    sub_matrix.k80_distance()
}

/// ## Kimura 3-Parameter (K-81) nucleotide substitution model.
///
/// Calculates evolutionary distance between two sequences, accounting for different rates
/// of **transitions** (A ↔ G or C ↔ T) and two types of
/// transversions, **type 1 transversions** (A ↔ T and C ↔ G), and **type 2 transversions**
/// (A ↔ C and G ↔ T) between sequences. Takes a [substitution matrix](dna_substitution_matrix).
///
/// See [Assumptions](super::dna).
///
/// The formula used is:
///
/// $$ d = -\frac{1}{4} \ln(a_1 \cdot a_2 \cdot a_3) $$
///
/// where $d$ is the evolutionary distance, and $a_1$, $a_2$, and $a_3$ are
/// intermediate variables defined as such:
/// $$ a_1 = 1 - 2p - 2q $$
/// $$ a_2 = 1 - 2p - 2r $$
/// $$ a_3 = 1 - 2q - 2r $$
/// where $p$ is the proportion of transitions, $\left( p = \frac{\text{transitions}}{\text{total bases}}\right)$
/// $q$ is the proportion of type 1 transversions, $\left( q = \frac{\text{transv. 1}}{\text{total bases}}\right)$
/// and $r$ is the proportion of type 2 transversions. $\left( r = \frac{\text{transv. 2}}{\text{total bases}}\right)$
///
/// ## Citations
///
/// - Kimura, M.  1981. "Estimation of evolutionary distances between homologous
///   nucleotide sequences."  Proc. Natl. Acad. Sci. U.S.A. 78, 454–458.
#[inline]
#[must_use]
pub fn kimura_81(seq1: &[u8], seq2: &[u8]) -> Option<f64> {
    if seq1.is_empty() || seq2.is_empty() {
        return None;
    }
    let sub_matrix = dna_substitution_matrix(seq1, seq2);
    sub_matrix.k81_distance()
}

/// ## Felsenstein (F-81) nucleotide substitution model.
///
/// Calculates evolutionary distance between two sequences, accounting for different
/// **base frequencies** of each nucleotide. Takes a [substitution matrix](dna_substitution_matrix).
///
/// See [Assumptions](super::dna)
///
/// The formula used is:
///
/// $$ d = -e \ln \left(1-\frac{p}{e} \right) $$
///
/// where $d$ is the evolutionary distance, e is an intermediate variable that represents the
/// effective number of possible substitutions:
///
/// $$ e = 1 - \sum{i=1}^{4} f_i^2 $$
///
/// and $p$ is the proportion of transitions, $\left( p = \frac{\text{transitions}}{\text{total bases}}\right)$.
///
/// ## Citations
///
/// - Felsenstein, J.  1981. "Evolutionary trees from DNA sequences: a maximum
///   likelihood approach."  J. Mol. Evol. 17, 368–376.
#[inline]
#[must_use]
pub fn felsenstein_81(seq1: &[u8], seq2: &[u8]) -> Option<f64> {
    if seq1.is_empty() || seq2.is_empty() {
        return None;
    }
    let sub_matrix = dna_substitution_matrix(seq1, seq2);
    sub_matrix.f81_distance()
}

/// ## Tamura-Nei (TN-93) nucleotide substitution model.
///
/// Calculates evolutionary distance between two sequences, accounting for different
/// **base frequencies** of each nucleotide, as well as different rates for different
/// substitution types: **transitions**, (A ↔ G or C ↔ T) **type 1 transversions**,
/// (A ↔ T and C ↔ G) and **type 2 transversions**. (A ↔ C and G ↔ T) Takes a [substitution matrix](dna_substitution_matrix).
///
/// See [Assumptions](super::dna)
///
/// The formula used is:
///
/// $$ d = -k_1 \ln(w_1) - k_2 \ln(w_2) - k_3 \ln(w_3) $$
///
/// where $d$ is the evolutionary distance, $k_1$, $k_2$, and $k_3$ are intermediate variables representing
/// nucleotide frequencies:
/// $$ k_1 = \frac{2 \cdot f_A \cdot f_G}{\text{Purines}} $$
/// $$ k_2 = \frac{2 \cdot f_C \cdot f_T}{\text{Pyrimidines}} $$
/// $$ k_3 = 2 \cdot ( \text{Purines}\cdot\text{Pyrimidines} - k_1 \cdot \text{Pyrimidines} - k_2 \cdot \text{Purines}) $$,
///
/// and $w_1$, $w_2$, and $w_3$ are intermediate variables that represent substitution rates *and*
/// substitution frequencies:
/// $$ w_1 = 1 - \frac{q}{k_1} - \frac{p}{2 \cdot \text{Purines}} $$
/// $$ w_2 = 1 - \frac{r}{k_2} - \frac{p}{2 \cdot \text{Pyrimidines}} $$
/// $$ w_3 = 1 - \frac{p}{2 \cdot \text{Purines} \cdot \text{Pyrimidines}} $$
///
/// where $p$ is the proportion of transitions, $\left( p = \frac{\text{transitions}}{\text{total bases}}\right)$
/// $q$ is the proportion of type 1 transversions, $\left( q = \frac{\text{transv. 1}}{\text{total bases}}\right)$
/// and $r$ is the proportion of type 2 transversions. $\left( r = \frac{\text{transv. 2}}{\text{total bases}}\right)$
///
/// ## Citations
///
/// - Tamura, K., and M. Nei.  1993. "Estimation of the number of nucleotide
///   substitutions in the control region of mitochondrial DNA in humans and
///   chimpanzees."  Mol. Biol. Evol. 10, 512–526.
#[inline]
#[must_use]
pub fn tamura_nei_93(seq1: &[u8], seq2: &[u8]) -> Option<f64> {
    if seq1.is_empty() || seq2.is_empty() {
        return None;
    }

    let sub_matrix = dna_substitution_matrix(seq1, seq2);
    sub_matrix.tn93_distance()
}

/// A trait for calculating evolutionary distances using the same substitution
/// matrix.
///
/// Substitution matrices are built from two `&[u8]` slices using
/// [`dna_substitution_matrix`]
pub trait DistanceFromMatrix: Sealed {
    /// Computes the JC69 distance from the given matrix. See
    /// [`jukes_cantor_69`] for more details.
    #[must_use]
    fn jc69_distance(&self) -> Option<f64>;

    /// Computes the K80 distance from the given matrix. See [`kimura_80`] for
    /// more details.
    #[must_use]
    fn k80_distance(&self) -> Option<f64>;

    /// Computes the K81 distance from the given matrix. See [`kimura_81`] for
    /// more details.
    #[must_use]
    fn k81_distance(&self) -> Option<f64>;

    /// Computes the F81 distance from the given matrix. See [`felsenstein_81`]
    /// for more details.
    #[must_use]
    fn f81_distance(&self) -> Option<f64>;

    /// Computes the TN93 distance from the given matrix. See [`tamura_nei_93`]
    /// for more details.
    #[must_use]
    fn tn93_distance(&self) -> Option<f64>;
}

impl DistanceFromMatrix for [[u32; 4]; 4] {
    #[inline]
    fn jc69_distance(&self) -> Option<f64> {
        let total_bases: u32 = self.iter().flatten().sum();
        let hamming = hamming_dist_from_sub_matrix(self);

        let p: f64 = f64::from(hamming) / f64::from(total_bases);
        (-0.75 * (1.0 - p * 4.0 / 3.0).ln()).into_option()
    }

    #[inline]
    fn k80_distance(&self) -> Option<f64> {
        let total_bases = self.iter().flatten().sum::<u32>();
        let synonymous = self[0][0] + self[1][1] + self[2][2] + self[3][3];
        let transitions = self[0][2] + self[2][0] + self[1][3] + self[3][1];
        let transversions = total_bases - synonymous - transitions;

        let p = f64::from(transitions) / f64::from(total_bases);
        let q = f64::from(transversions) / f64::from(total_bases);

        (-0.5 * (1.0 - 2.0 * p - q).ln() - 0.25 * (1.0 - 2.0 * q).ln()).into_option()
    }

    #[inline]
    fn k81_distance(&self) -> Option<f64> {
        let total_bases = f64::from(self.iter().flatten().sum::<u32>());
        let transv1 = f64::from(self[0][3] + self[3][0] + self[1][2] + self[2][1]);
        let transv2 = f64::from(self[0][1] + self[1][0] + self[2][3] + self[3][2]);
        let transitions = f64::from(self[1][3] + self[3][1] + self[0][2] + self[2][0]);

        let p = transitions / total_bases;
        let q = transv1 / total_bases;
        let r = transv2 / total_bases;

        let a1 = 1.0 - 2.0 * p - 2.0 * q;
        let a2 = 1.0 - 2.0 * p - 2.0 * r;
        let a3 = 1.0 - 2.0 * q - 2.0 * r;

        (-0.25 * (a1 * a2 * a3).ln()).into_option()
    }

    #[inline]
    fn f81_distance(&self) -> Option<f64> {
        let (total_bases, bf) = total_and_frequencies(self);

        let hamming = f64::from(hamming_dist_from_sub_matrix(self));

        let e = 1.0 - bf.iter().map(|&freq| freq.powi(2)).sum::<f64>();

        let p = hamming / f64::from(total_bases);

        (-e * (1.0 - p / e).ln()).into_option()
    }

    #[inline]
    fn tn93_distance(&self) -> Option<f64> {
        let (total_bases, bf) = total_and_frequencies(self);
        let total_bases = f64::from(total_bases);

        let purines = bf[0] + bf[2];
        let pyrimidines = bf[1] + bf[3];

        let k1 = 2.0 * bf[0] * bf[2] / purines;
        let k2 = 2.0 * bf[1] * bf[3] / pyrimidines;
        let k3 =
            2.0 * (purines * pyrimidines - bf[0] * bf[2] * pyrimidines / purines - bf[1] * bf[3] * purines / pyrimidines);

        let hamming = f64::from(hamming_dist_from_sub_matrix(self));

        let transv1 = f64::from(self[0][2] + self[2][0]);
        let transv2 = f64::from(self[1][3] + self[3][1]);

        let p1 = transv1 / total_bases;
        let p2 = transv2 / total_bases;
        let q = (hamming - transv1 - transv2) / total_bases;

        let w1 = 1.0 - p1 / k1 - q / (2.0 * purines);
        let w2 = 1.0 - p2 / k2 - q / (2.0 * pyrimidines);
        let w3 = 1.0 - q / (2.0 * purines * pyrimidines);

        (-k1 * w1.ln() - k2 * w2.ln() - k3 * w3.ln()).into_option()
    }
}
