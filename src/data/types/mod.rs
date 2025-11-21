/// Structs and types for maniuplating amino acid sequences.
pub mod amino_acids;
/// Structs and types for CIGAR strings.
pub mod cigar;
/// Structs and types for nucleotide sequence manipulation.
pub mod nucleotides;
/// Structs and types for processing Phred quality scores.
pub mod phred;

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
                if $crate::data::validation::check::CheckSequence::is_ascii_simd::<16>(&self.0) {
                    // SAFETY: we just checked it is ASCII using our fast SIMD function.
                    // ASCII is valid UTF8.
                    f.write_str(unsafe { std::str::from_utf8_unchecked(&self.0) })
                } else {
                    f.write_str(&String::from_utf8_lossy(&self.0))
                }
            }
        }

        impl std::fmt::Display for $view<'_> {
            #[inline]
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                if $crate::data::validation::check::CheckSequence::is_ascii_simd::<16>(&self.0) {
                    // SAFETY: we just checked it is ASCII using our fast SIMD function.
                    // ASCII is valid UTF8.
                    f.write_str(unsafe { std::str::from_utf8_unchecked(self.0) })
                } else {
                    f.write_str(&String::from_utf8_lossy(self.0))
                }
            }
        }

        impl std::fmt::Display for $viewmut<'_> {
            #[inline]
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                if $crate::data::validation::check::CheckSequence::is_ascii_simd::<16>(&self.0) {
                    // SAFETY: we just checked it is ASCII using our fast SIMD function.
                    // ASCII is valid UTF8.
                    f.write_str(unsafe { std::str::from_utf8_unchecked(self.0) })
                } else {
                    f.write_str(&String::from_utf8_lossy(self.0))
                }
            }
        }
    };
}

pub(crate) use impl_std_traits_for_sequence;
