use crate::data::err::{DistanceError, DistanceError::*};

/// Calculates a "physiochemical" distance measure using the euclidean distances
/// over physiochemical factors. Only valid amino acids are permitted in the
/// denominator.
///
/// # Example
/// ```
/// use zoe::distance::physiochemical;
///
/// let s1 = b"MANATEE";
/// let s2 = b"MANGAEE";
///
/// assert!( (physiochemical(s1,s2).unwrap()  - 1.1464952).abs() < 0.01 )
/// ```
///
/// # Citation
/// For factor analysis used by the function:
///
/// 1. Atchley et al. 2008. "Solving the protein sequence metric problem." Proc
///    Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Published 2005 Apr 25.
///
///
/// # Errors
///
/// If either argument is empty or the sequence characters are invalid, an error
/// is thrown. See [`DistanceError`].
//
// TO-DO: Could make a generic function that depends on distance matrix.
#[allow(clippy::cast_precision_loss)]
pub fn physiochemical(seq1: &[u8], seq2: &[u8]) -> Result<f32, DistanceError> {
    use crate::data::matrices::PHYSIOCHEMICAL_FACTORS as dm;

    if seq1.is_empty() || seq2.is_empty() {
        return Err(NoData);
    }

    if seq1 == seq2 {
        return Ok(0.0);
    }

    let (number_valid, mut pcd_distance) = seq1
        .iter()
        .zip(seq2)
        .filter_map(|(a1, a2)| dm[*a1 as usize][*a2 as usize])
        .map(|val| (1, val))
        .fold((0u32, 0.0), |(a_valid, a_value), (i_valid, i_value)| {
            (a_valid + i_valid, a_value + i_value)
        });

    if number_valid > 0 {
        pcd_distance /= number_valid as f32;
        Ok(pcd_distance)
    } else {
        Err(NotComparable)
    }
}
