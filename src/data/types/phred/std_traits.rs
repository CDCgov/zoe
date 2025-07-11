use crate::{
    data::{types::phred::EncodedQS, validation::CheckSequence},
    prelude::*,
};
use std::{
    io::{Error as IOError, ErrorKind},
    ops::Index,
    slice::SliceIndex,
};

impl AsRef<[u8]> for QualityScores {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsRef<[u8]> for QualityScoresView<'_> {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl AsRef<[u8]> for QualityScoresViewMut<'_> {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl AsMut<[u8]> for QualityScores {
    #[inline]
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl AsMut<[u8]> for QualityScoresViewMut<'_> {
    #[inline]
    fn as_mut(&mut self) -> &mut [u8] {
        self.0
    }
}

impl AsRef<Vec<u8>> for QualityScores {
    #[inline]
    fn as_ref(&self) -> &Vec<u8> {
        &self.0
    }
}

impl AsMut<Vec<u8>> for QualityScores {
    #[inline]
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }
}

impl TryFrom<String> for QualityScores {
    type Error = IOError;
    #[inline]
    fn try_from(encoded_scores: String) -> Result<Self, Self::Error> {
        QualityScores::try_from(encoded_scores.into_bytes())
    }
}

impl TryFrom<Vec<EncodedQS>> for QualityScores {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: Vec<EncodedQS>) -> Result<Self, Self::Error> {
        if encoded_scores.is_graphic_simd::<16>() {
            Ok(QualityScores(encoded_scores))
        } else {
            Err(IOError::new(ErrorKind::InvalidData, "Quality scores contain invalid state!"))
        }
    }
}

impl<const N: usize> TryFrom<[EncodedQS; N]> for QualityScores {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: [EncodedQS; N]) -> Result<Self, Self::Error> {
        QualityScores::try_from(encoded_scores.to_vec())
    }
}

impl TryFrom<&[EncodedQS]> for QualityScores {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &[EncodedQS]) -> Result<Self, Self::Error> {
        QualityScores::try_from(encoded_scores.to_vec())
    }
}

impl<'a> TryFrom<&'a [EncodedQS]> for QualityScoresView<'a> {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &'a [EncodedQS]) -> Result<Self, Self::Error> {
        if encoded_scores.is_graphic_simd::<16>() {
            Ok(QualityScoresView(encoded_scores))
        } else {
            Err(IOError::new(ErrorKind::InvalidData, "Quality scores contain invalid state!"))
        }
    }
}

impl TryFrom<&mut [EncodedQS]> for QualityScores {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &mut [EncodedQS]) -> Result<Self, Self::Error> {
        QualityScores::try_from(encoded_scores.to_vec())
    }
}

impl<'a> TryFrom<&'a mut [EncodedQS]> for QualityScoresView<'a> {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &'a mut [EncodedQS]) -> Result<Self, Self::Error> {
        QualityScoresView::try_from(&*encoded_scores)
    }
}

impl<'a> TryFrom<&'a mut [EncodedQS]> for QualityScoresViewMut<'a> {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &'a mut [EncodedQS]) -> Result<Self, Self::Error> {
        if encoded_scores.is_graphic_simd::<16>() {
            Ok(QualityScoresViewMut(encoded_scores))
        } else {
            Err(IOError::new(ErrorKind::InvalidData, "Quality scores contain invalid state!"))
        }
    }
}

impl<const N: usize> TryFrom<&[EncodedQS; N]> for QualityScores {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &[EncodedQS; N]) -> Result<Self, Self::Error> {
        QualityScores::try_from(encoded_scores.to_vec())
    }
}

impl<'a, const N: usize> TryFrom<&'a [EncodedQS; N]> for QualityScoresView<'a> {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &'a [EncodedQS; N]) -> Result<Self, Self::Error> {
        QualityScoresView::try_from(encoded_scores.as_slice())
    }
}

impl<const N: usize> TryFrom<&mut [EncodedQS; N]> for QualityScores {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &mut [EncodedQS; N]) -> Result<Self, Self::Error> {
        QualityScores::try_from(encoded_scores.to_vec())
    }
}

impl<'a, const N: usize> TryFrom<&'a mut [EncodedQS; N]> for QualityScoresView<'a> {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &'a mut [EncodedQS; N]) -> Result<Self, Self::Error> {
        QualityScoresView::try_from(encoded_scores.as_slice())
    }
}

impl<'a, const N: usize> TryFrom<&'a mut [EncodedQS; N]> for QualityScoresViewMut<'a> {
    type Error = IOError;

    #[inline]
    fn try_from(encoded_scores: &'a mut [EncodedQS; N]) -> Result<Self, Self::Error> {
        QualityScoresViewMut::try_from(encoded_scores.as_mut_slice())
    }
}

impl IntoIterator for QualityScores {
    type Item = u8;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for QualityScoresView<'a> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for QualityScoresViewMut<'a> {
    type Item = &'a mut u8;
    type IntoIter = std::slice::IterMut<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<'a> IntoIterator for &'a QualityScores {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a QualityScoresView<'_> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a QualityScoresViewMut<'_> {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for QualityScores {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for QualityScoresView<'_> {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<I: SliceIndex<[u8]>> Index<I> for QualityScoresViewMut<'_> {
    type Output = <I as SliceIndex<[u8]>>::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl std::fmt::Display for QualityScores {
    // Safety: `QualityScores` is guaranteed to be utf8 because we checked that
    // it is "graphic ascii", which will be valid utf8.
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(unsafe { std::str::from_utf8_unchecked(&self.0) })
    }
}

impl std::fmt::Display for QualityScoresView<'_> {
    // Safety: `QualityScores` is guaranteed to be utf8 because we checked that
    // it is "graphic ascii", which will be valid utf8.
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(unsafe { std::str::from_utf8_unchecked(self.0) })
    }
}

impl std::fmt::Display for QualityScoresViewMut<'_> {
    // Safety: `QualityScores` is guaranteed to be utf8 because we checked that
    // it is "graphic ascii", which will be valid utf8.
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(unsafe { std::str::from_utf8_unchecked(self.0) })
    }
}

impl std::fmt::Debug for QualityScores {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl std::fmt::Debug for QualityScoresView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

impl std::fmt::Debug for QualityScoresViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self}")
    }
}
