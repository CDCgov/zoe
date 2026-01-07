use crate::{
    DEFAULT_SIMD_LANES,
    data::types::amino_acids::AminoAcidsReadable,
    distance::{
        DistanceError::{self, NoData, NotComparable},
        hamming_simd,
    },
    prelude::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
    private::Sealed,
};

pub trait AminoAcidsDistance: AminoAcidsReadable + Sealed {
    /// ## Distance
    ///
    /// Calculates hamming distance between `self` and another sequence.
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     data::types::amino_acids::{AminoAcids, AminoAcidsView},
    /// #     distance::aa::AminoAcidsDistance
    /// # };
    ///
    /// let s1: AminoAcids = b"MANATEEMANATEEMANATEE".into();
    /// let s2: AminoAcids = b"MANGAEEMANATEEMANGAEE".into();
    ///
    /// assert!(4 == s1.distance_hamming(&s2));
    ///
    /// let s3 = AminoAcidsView::from_bytes_unchecked(b"MANGAEEMANATEEMANGAEE");
    /// assert!(4 == s1.distance_hamming(&s3));
    /// ```
    #[inline]
    #[must_use]
    fn distance_hamming<T: AminoAcidsReadable>(&self, other_sequence: &T) -> usize {
        hamming_simd::<{ DEFAULT_SIMD_LANES }>(self.amino_acids_bytes(), other_sequence.amino_acids_bytes())
    }

    /// Calculates physiochemical distance between `self`and another protein
    /// sequence. See: [`crate::distance::aa::physiochemical`]
    ///
    /// ## Citation
    ///
    /// For factor analysis used by the function:
    ///
    /// > Atchley et al. 2008. "Solving the protein sequence metric problem."
    /// > Proc Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Published 2005
    /// > Apr 25.
    ///
    ///
    /// ## Errors
    ///
    /// If either argument is empty or the sequence characters are invalid, an
    /// error is thrown. See [`DistanceError`].
    ///
    /// ## Example
    ///
    /// ```
    /// # use zoe::{
    /// #     data::types::amino_acids::{AminoAcids, AminoAcidsView},
    /// #     distance::aa::AminoAcidsDistance
    /// # };
    ///
    /// let s1: AminoAcids = b"MANATEEMANATEEMANATEE".into();
    /// let s2: AminoAcids = b"MANGAEEMANATEEMANGAEE".into();
    ///
    /// assert!( (s1.distance_physiochemical(&s2).unwrap() - 0.7643302).abs() < 0.01 );
    ///
    /// let s3 = AminoAcidsView::from_bytes_unchecked(b"MANGAEEMANATEEMANGAEE");
    /// assert!( (s1.distance_physiochemical(&s3).unwrap() - 0.7643302).abs() < 0.01 );
    /// ```
    #[inline]
    fn distance_physiochemical<T: AminoAcidsReadable>(&self, other_sequence: &T) -> Result<f32, DistanceError> {
        physiochemical(self.amino_acids_bytes(), other_sequence.amino_acids_bytes())
    }
}

impl AminoAcidsDistance for AminoAcids {}
impl AminoAcidsDistance for AminoAcidsView<'_> {}
impl AminoAcidsDistance for AminoAcidsViewMut<'_> {}

/// Calculates a "physiochemical" distance measure using the euclidean distances
/// over physiochemical factors. Only valid amino acids are permitted in the
/// denominator.
///
/// ## Example
///
/// ```
/// # use zoe::distance::aa::physiochemical;
///
/// let s1 = b"MANATEE";
/// let s2 = b"MANGAEE";
///
/// assert!( (physiochemical(s1,s2).unwrap()  - 1.1464952).abs() < 0.01 )
/// ```
///
/// ## Citation
///
/// For factor analysis used by the function:
///
/// 1. Atchley et al. 2008. "Solving the protein sequence metric problem." Proc
///    Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Published 2005 Apr 25.
///
///
/// ## Errors
///
/// If either argument is empty or the sequence characters are invalid, an error
/// is thrown. See [`DistanceError`].
//
// TODO: Could make a generic function that depends on distance matrix.
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
