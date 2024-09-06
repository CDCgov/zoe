use crate::data::constants::mappings::to_dna_index;

#[must_use]
/// ## Generates a 4x4 substitution matrix from two aligned sequences of the type `&[u8]`.
///
/// Each column refers to nucleotide bases in the first sequence (A, C, G, T, respectively),
/// and each row refers to nucleotide bases in the second sequence.
///
/// The substitution matrix can then be indexed to find counts of each type of substitution
/// between the sequences, or used to calculate evolutionary distances.
///
/// ### Example
/// ```
/// # use zoe::distance::dna_substitution_matrix;
/// let seq1: &[u8] = b"GATCAGATTTGCATTGGTT";
/// let seq2: &[u8] = b"GATCATATTAGCATTGCTT";
///
/// let sub_matrix: [[u32; 4]; 4] = dna_substitution_matrix(seq1, seq2);
/// ```
/// This should return the matrix (bases provided here for context)
/// $$ \begin{pmatrix}  & \text{A} & \text{C} & \text{G} & \text{T} \cr \text{A} & 4 & 0 & 0 & 1 \cr \text{C} & 0 & 2 & 1 & 0 \cr \text{G} & 0 & 0 & 3 & 0 \cr \text{T} & 0 & 0 & 1 & 7 \end{pmatrix} $$
pub fn dna_substitution_matrix(seq1: &[u8], seq2: &[u8]) -> [[u32; 4]; 4] {
    let mut sub_matrix = [[0u32; 4]; 4];
    std::iter::zip(seq1.iter().copied().map(to_dna_index), seq2.iter().copied().map(to_dna_index))
        .filter(|(a, b)| *a < 4 && *b < 4)
        .for_each(|(a, b)| sub_matrix[a as usize][b as usize] += 1);
    sub_matrix
}

#[allow(clippy::needless_range_loop)]
#[inline]
#[must_use]
/// ## Calculates the Hamming distance from a nucleotide [substitution matrix](dna_substitution_matrix).
///
/// The Hamming distance is the number of substitutions between
/// aligned sequences, or sum of non-diagonal values in the substitution matrix.
pub fn hamming_dist_from_sub_matrix(sub_matrix: &[[u32; 4]; 4]) -> u32 {
    let mut hamming = 0;

    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                hamming += sub_matrix[i][j];
            }
        }
    }
    hamming
}

#[must_use]
/// ## Finds total number of bases and nucleotide base frequencies from a [substitution matrix](dna_substitution_matrix).
///
/// Returns a tuple of the total number of nucleotides in the two sequences,
/// and an array of base frequencies, corresponding to the
/// frequencies of A, C, G, and T bases in the sequences, respectively.
pub fn total_and_frequencies(sub_matrix: &[[u32; 4]; 4]) -> (u32, [f64; 4]) {
    let mut bf = [0f64; 4];
    for i in 0..4 {
        // sum of ith column and row to get base frequencies
        bf[i] = f64::from(sub_matrix[i].iter().sum::<u32>() + sub_matrix.iter().map(|row| row[i]).sum::<u32>());
    }
    let total_bases = sub_matrix.iter().flatten().sum::<u32>();
    // counts of nucleotides are divided by 2 to remove double counting
    bf.iter_mut().for_each(|freq| *freq /= 2.0 * f64::from(total_bases));

    (total_bases, bf)
}
