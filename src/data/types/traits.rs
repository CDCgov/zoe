use crate::data::{
    amino_acids::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
    nucleotides::{Nucleotides, NucleotidesView, NucleotidesViewMut},
    views::{AssocOwnedType, ToView},
};
use std::ops::Deref;

/// A trait unifying different byte sequence types, like `Vec<u8>`, `&[u8]`,
/// [`Nucleotides`], [`AminoAcids`], and view types.
///
/// Specifically, this is a type which can be viewed as `&[u8]` and from which
/// an owned version can be created from a `Vec<u8>`. It also can be converted
/// into either `&[u8]` (for borrowed data) or `Vec<u8>` (for owned data) with
/// [`into_bytes`], and can be converted into an iterator of `u8` with
/// [`into_bytes_iter`].
///
/// [`Nucleotides`]: crate::data::types::nucleotides::Nucleotides
/// [`AminoAcids`]: crate::data::types::amino_acids::AminoAcids
/// [`into_bytes`]: ByteSeq::into_bytes
/// [`into_bytes_iter`]: ByteSeq::into_bytes_iter
pub trait ByteSeq: AsRef<[u8]> + AssocOwnedType<Owned: From<Vec<u8>>> {
    /// The contained bytes within the [`ByteSeq`], either `&[u8]` or `Vec<u8>`.
    type Bytes: Deref<Target = [u8]>;
    /// An iterator over the bytes in [`ByteSeq`].
    type IntoIter: Iterator<Item = u8>;

    /// Extracts the bytes from the [`ByteSeq`] without cloning, returning
    /// either `&[u8]` for borrowed data or `Vec<u8>` for owned data.
    fn into_bytes(self) -> Self::Bytes;

    /// Returns an iterator of the bytes within the [`ByteSeq`].
    fn into_bytes_iter(self) -> Self::IntoIter;
}

impl<'a, T: ByteSeq> ByteSeq for &'a T {
    type Bytes = &'a [u8];
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, u8>>;

    fn into_bytes(self) -> Self::Bytes {
        self.as_ref()
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        self.as_ref().iter().copied()
    }
}

impl<'a, T: ByteSeq> ByteSeq for &'a mut T {
    type Bytes = &'a [u8];
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, u8>>;

    fn into_bytes(self) -> Self::Bytes {
        (&*self).into_bytes()
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        (&*self).into_bytes_iter()
    }
}

impl ByteSeq for Vec<u8> {
    type Bytes = Vec<u8>;
    type IntoIter = std::vec::IntoIter<u8>;

    fn into_bytes(self) -> Self::Bytes {
        self
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        self.into_iter()
    }
}

impl<'a> ByteSeq for &'a [u8] {
    type Bytes = &'a [u8];
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, u8>>;

    fn into_bytes(self) -> Self::Bytes {
        self
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        self.iter().copied()
    }
}

impl<'a> ByteSeq for &'a mut [u8] {
    type Bytes = &'a [u8];
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, u8>>;

    fn into_bytes(self) -> Self::Bytes {
        self
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        self.iter().copied()
    }
}

impl ByteSeq for Nucleotides {
    type Bytes = Vec<u8>;
    type IntoIter = std::vec::IntoIter<u8>;

    fn into_bytes(self) -> Self::Bytes {
        self.into()
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        self.into_iter()
    }
}

impl<'a> ByteSeq for NucleotidesView<'a> {
    type Bytes = &'a [u8];
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, u8>>;

    fn into_bytes(self) -> Self::Bytes {
        self.into()
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        self.iter().copied()
    }
}

impl<'a> ByteSeq for NucleotidesViewMut<'a> {
    type Bytes = &'a [u8];
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, u8>>;

    fn into_bytes(self) -> Self::Bytes {
        self.to_view().into()
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        <&mut [_]>::from(self).iter().copied()
    }
}

impl ByteSeq for AminoAcids {
    type Bytes = Vec<u8>;
    type IntoIter = std::vec::IntoIter<u8>;

    fn into_bytes(self) -> Self::Bytes {
        self.into()
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        self.into_iter()
    }
}

impl<'a> ByteSeq for AminoAcidsView<'a> {
    type Bytes = &'a [u8];
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, u8>>;

    fn into_bytes(self) -> Self::Bytes {
        self.into()
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        self.iter().copied()
    }
}

impl<'a> ByteSeq for AminoAcidsViewMut<'a> {
    type Bytes = &'a [u8];
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, u8>>;

    fn into_bytes(self) -> Self::Bytes {
        self.to_view().into()
    }

    fn into_bytes_iter(self) -> Self::IntoIter {
        <&mut [_]>::from(self).iter().copied()
    }
}

/// A macro for implementing standard library traits on sequence types and
/// views.
macro_rules! impl_std_traits_for_sequence {
    ($owned:ident, $view:ident, $viewmut:ident) => {
        impl AsRef<[u8]> for $owned {
            #[inline]
            fn as_ref(&self) -> &[u8] {
                &self.0
            }
        }

        impl AsRef<Vec<u8>> for $owned {
            #[inline]
            fn as_ref(&self) -> &Vec<u8> {
                &self.0
            }
        }

        impl AsRef<[u8]> for $view<'_> {
            #[inline]
            fn as_ref(&self) -> &[u8] {
                self.0
            }
        }

        impl AsRef<[u8]> for $viewmut<'_> {
            #[inline]
            fn as_ref(&self) -> &[u8] {
                self.0
            }
        }

        impl AsMut<[u8]> for $owned {
            #[inline]
            fn as_mut(&mut self) -> &mut [u8] {
                &mut self.0
            }
        }

        impl AsMut<Vec<u8>> for $owned {
            #[inline]
            fn as_mut(&mut self) -> &mut Vec<u8> {
                &mut self.0
            }
        }

        impl AsMut<[u8]> for $viewmut<'_> {
            #[inline]
            fn as_mut(&mut self) -> &mut [u8] {
                self.0
            }
        }

        impl From<Vec<u8>> for $owned {
            #[inline]
            fn from(vec: Vec<u8>) -> Self {
                Self(vec)
            }
        }

        impl From<&Vec<u8>> for $owned {
            #[inline]
            fn from(vec: &Vec<u8>) -> Self {
                Self::from(vec.clone())
            }
        }

        impl From<&mut Vec<u8>> for $owned {
            #[inline]
            fn from(vec: &mut Vec<u8>) -> Self {
                Self::from(vec.clone())
            }
        }

        impl From<&[u8]> for $owned {
            #[inline]
            fn from(bytes: &[u8]) -> Self {
                Self::from(bytes.to_vec())
            }
        }

        impl From<&mut [u8]> for $owned {
            #[inline]
            fn from(bytes: &mut [u8]) -> Self {
                Self::from(bytes.to_vec())
            }
        }

        impl<const N: usize> From<[u8; N]> for $owned {
            #[inline]
            fn from(bytes: [u8; N]) -> Self {
                Self::from(bytes.to_vec())
            }
        }

        impl<const N: usize> From<&[u8; N]> for $owned {
            #[inline]
            fn from(bytes: &[u8; N]) -> Self {
                Self::from(bytes.to_vec())
            }
        }

        impl<const N: usize> From<&mut [u8; N]> for $owned {
            #[inline]
            fn from(bytes: &mut [u8; N]) -> Self {
                Self::from(bytes.to_vec())
            }
        }

        impl From<String> for $owned {
            #[inline]
            fn from(s: String) -> Self {
                Self::from(s.into_bytes())
            }
        }

        impl From<&String> for $owned {
            #[inline]
            fn from(s: &String) -> Self {
                Self::from(s.as_bytes())
            }
        }

        impl From<&mut String> for $owned {
            #[inline]
            fn from(s: &mut String) -> Self {
                Self::from(s.as_bytes())
            }
        }

        impl From<&str> for $owned {
            #[inline]
            fn from(s: &str) -> Self {
                Self::from(s.as_bytes())
            }
        }

        impl From<&mut str> for $owned {
            #[inline]
            fn from(s: &mut str) -> Self {
                Self::from(s.as_bytes())
            }
        }

        impl<'a> From<&'a Vec<u8>> for $view<'a> {
            #[inline]
            fn from(bytes: &'a Vec<u8>) -> Self {
                Self(bytes)
            }
        }

        impl<'a> From<&'a mut Vec<u8>> for $view<'a> {
            #[inline]
            fn from(bytes: &'a mut Vec<u8>) -> Self {
                Self(bytes)
            }
        }

        impl<'a> From<&'a [u8]> for $view<'a> {
            #[inline]
            fn from(bytes: &'a [u8]) -> Self {
                Self(bytes)
            }
        }

        impl<'a> From<&'a mut [u8]> for $view<'a> {
            #[inline]
            fn from(bytes: &'a mut [u8]) -> Self {
                Self(bytes)
            }
        }

        impl<'a, const N: usize> From<&'a [u8; N]> for $view<'a> {
            #[inline]
            fn from(bytes: &'a [u8; N]) -> Self {
                Self(bytes)
            }
        }

        impl<'a, const N: usize> From<&'a mut [u8; N]> for $view<'a> {
            #[inline]
            fn from(bytes: &'a mut [u8; N]) -> Self {
                Self(bytes)
            }
        }

        impl<'a> From<&'a str> for $view<'a> {
            #[inline]
            fn from(s: &'a str) -> Self {
                Self(s.as_bytes())
            }
        }

        impl<'a> From<&'a mut str> for $view<'a> {
            #[inline]
            fn from(s: &'a mut str) -> Self {
                Self(s.as_bytes())
            }
        }

        impl<'a> From<&'a mut Vec<u8>> for $viewmut<'a> {
            #[inline]
            fn from(bytes: &'a mut Vec<u8>) -> Self {
                Self(bytes)
            }
        }

        impl<'a> From<&'a mut [u8]> for $viewmut<'a> {
            #[inline]
            fn from(bytes: &'a mut [u8]) -> Self {
                Self(bytes)
            }
        }

        impl<'a, const N: usize> From<&'a mut [u8; N]> for $viewmut<'a> {
            #[inline]
            fn from(bytes: &'a mut [u8; N]) -> Self {
                Self(bytes)
            }
        }

        impl From<$owned> for Vec<u8> {
            #[inline]
            fn from(value: $owned) -> Vec<u8> {
                value.0
            }
        }

        impl<'a> From<&'a $owned> for &'a [u8] {
            #[inline]
            fn from(value: &'a $owned) -> &'a [u8] {
                &value.0
            }
        }

        impl<'a> From<$view<'a>> for &'a [u8] {
            #[inline]
            fn from(value: $view<'a>) -> &'a [u8] {
                value.0
            }
        }

        impl<'a> From<$viewmut<'a>> for &'a mut [u8] {
            #[inline]
            fn from(value: $viewmut<'a>) -> &'a mut [u8] {
                value.0
            }
        }

        impl<'a> From<$owned> for ::std::borrow::Cow<'a, [u8]> {
            #[inline]
            fn from(value: $owned) -> ::std::borrow::Cow<'a, [u8]> {
                value.0.into()
            }
        }

        impl<'a> From<&'a $owned> for ::std::borrow::Cow<'a, [u8]> {
            #[inline]
            fn from(value: &'a $owned) -> ::std::borrow::Cow<'a, [u8]> {
                (&value.0).into()
            }
        }

        impl<'a> From<$view<'a>> for ::std::borrow::Cow<'a, [u8]> {
            #[inline]
            fn from(value: $view<'a>) -> ::std::borrow::Cow<'a, [u8]> {
                value.0.into()
            }
        }

        impl IntoIterator for $owned {
            type Item = u8;
            type IntoIter = std::vec::IntoIter<Self::Item>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                self.0.into_iter()
            }
        }

        impl<'a> IntoIterator for &'a $owned {
            type Item = &'a u8;
            type IntoIter = std::slice::Iter<'a, u8>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                self.0.iter()
            }
        }

        impl<'a> IntoIterator for &'a mut $owned {
            type Item = &'a mut u8;
            type IntoIter = std::slice::IterMut<'a, u8>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                self.0.iter_mut()
            }
        }

        impl<'a> IntoIterator for $view<'a> {
            type Item = &'a u8;
            type IntoIter = std::slice::Iter<'a, u8>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                self.0.iter()
            }
        }

        impl<'a> IntoIterator for &'a $view<'_> {
            type Item = &'a u8;
            type IntoIter = std::slice::Iter<'a, u8>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                self.0.iter()
            }
        }

        impl<'a> IntoIterator for $viewmut<'a> {
            type Item = &'a mut u8;
            type IntoIter = std::slice::IterMut<'a, u8>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                self.0.iter_mut()
            }
        }

        impl<'a> IntoIterator for &'a $viewmut<'_> {
            type Item = &'a u8;
            type IntoIter = std::slice::Iter<'a, u8>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                self.0.iter()
            }
        }

        impl<'a> IntoIterator for &'a mut $viewmut<'_> {
            type Item = &'a mut u8;
            type IntoIter = std::slice::IterMut<'a, u8>;

            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                self.0.iter_mut()
            }
        }

        impl FromIterator<u8> for $owned {
            #[inline]
            fn from_iter<T: IntoIterator<Item = u8>>(iterable: T) -> Self {
                Self(iterable.into_iter().collect())
            }
        }

        impl<I: ::std::slice::SliceIndex<[u8]>> ::std::ops::Index<I> for $owned {
            type Output = <I as ::std::slice::SliceIndex<[u8]>>::Output;

            #[inline]
            fn index(&self, index: I) -> &Self::Output {
                &self.0[index]
            }
        }

        impl<I: ::std::slice::SliceIndex<[u8]>> ::std::ops::Index<I> for $view<'_> {
            type Output = <I as ::std::slice::SliceIndex<[u8]>>::Output;

            #[inline]
            fn index(&self, index: I) -> &Self::Output {
                &self.0[index]
            }
        }

        impl<I: ::std::slice::SliceIndex<[u8]>> ::std::ops::Index<I> for $viewmut<'_> {
            type Output = <I as ::std::slice::SliceIndex<[u8]>>::Output;

            #[inline]
            fn index(&self, index: I) -> &Self::Output {
                &self.0[index]
            }
        }

        impl<I: ::std::slice::SliceIndex<[u8]>> ::std::ops::IndexMut<I> for $owned {
            #[inline]
            fn index_mut(&mut self, index: I) -> &mut Self::Output {
                &mut self.0[index]
            }
        }

        impl<I: ::std::slice::SliceIndex<[u8]>> ::std::ops::IndexMut<I> for $viewmut<'_> {
            #[inline]
            fn index_mut(&mut self, index: I) -> &mut Self::Output {
                &mut self.0[index]
            }
        }

        impl std::fmt::Display for $owned {
            #[inline]
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(
                    f,
                    "{}",
                    $crate::data::validation::check::SimdDisplay::display_ascii_or_lossy(&self.0)
                )
            }
        }

        impl std::fmt::Display for $view<'_> {
            #[inline]
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(
                    f,
                    "{}",
                    $crate::data::validation::check::SimdDisplay::display_ascii_or_lossy(&self.0)
                )
            }
        }

        impl std::fmt::Display for $viewmut<'_> {
            #[inline]
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(
                    f,
                    "{}",
                    $crate::data::validation::check::SimdDisplay::display_ascii_or_lossy(&self.0)
                )
            }
        }
    };
}

pub(crate) use impl_std_traits_for_sequence;
