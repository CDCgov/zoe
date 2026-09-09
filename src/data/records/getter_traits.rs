//! Getter traits for fields in record types, such as headers, annotations, and
//! sequences.

use std::ops::Bound;

use crate::{
    alignment::{LocalProfiles, SharedProfiles},
    data::{
        amino_acids::{AminoAcids, AminoAcidsView, AminoAcidsViewMut},
        fasta::{FastaAA, FastaNT, FastaNTAnnot, FastaSeq},
        fastq::{FastQ, FastQView, FastQViewMut},
        nucleotides::{Nucleotides, NucleotidesView, NucleotidesViewMut},
    },
    search::ToStrRangeSearch,
};

/// Getter trait for structures providing read access to a header/name.
pub trait HeaderReadable {
    /// Gets the header from the record.
    #[must_use]
    fn header(&self) -> &str;
}

impl HeaderReadable for FastQ {
    #[inline]
    fn header(&self) -> &str {
        &self.header
    }
}

impl HeaderReadable for FastQView<'_> {
    #[inline]
    fn header(&self) -> &str {
        self.header
    }
}

impl HeaderReadable for FastQViewMut<'_> {
    #[inline]
    fn header(&self) -> &str {
        self.header
    }
}

impl HeaderReadable for FastaSeq {
    #[inline]
    fn header(&self) -> &str {
        &self.name
    }
}

impl HeaderReadable for FastaNT {
    #[inline]
    fn header(&self) -> &str {
        &self.name
    }
}

impl HeaderReadable for FastaAA {
    #[inline]
    fn header(&self) -> &str {
        &self.name
    }
}

impl HeaderReadable for FastaNTAnnot {
    #[inline]
    fn header(&self) -> &str {
        &self.name
    }
}

/// Getter trait for structures providing mutable access to a header/name.
pub trait HeaderMutable {
    /// Gets the header from the record.
    #[must_use]
    fn header_mut(&mut self) -> &mut String;
}

impl HeaderMutable for FastQ {
    #[inline]
    fn header_mut(&mut self) -> &mut String {
        &mut self.header
    }
}

impl HeaderMutable for FastQViewMut<'_> {
    #[inline]
    fn header_mut(&mut self) -> &mut String {
        self.header
    }
}

impl HeaderMutable for FastaSeq {
    #[inline]
    fn header_mut(&mut self) -> &mut String {
        &mut self.name
    }
}

impl HeaderMutable for FastaNT {
    #[inline]
    fn header_mut(&mut self) -> &mut String {
        &mut self.name
    }
}

impl HeaderMutable for FastaAA {
    #[inline]
    fn header_mut(&mut self) -> &mut String {
        &mut self.name
    }
}

impl HeaderMutable for FastaNTAnnot {
    #[inline]
    fn header_mut(&mut self) -> &mut String {
        &mut self.name
    }
}

/// Getter trait for structures providing read access to a sequence.
///
/// The sequence can be either nucleotides or amino acids, and is returned as a
/// byte slice.
pub trait SequenceReadable {
    /// Get the sequence from the struct as a byte slice.
    #[must_use]
    fn sequence_bytes(&self) -> &[u8];
}

impl SequenceReadable for Nucleotides {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for NucleotidesView<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for NucleotidesViewMut<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for AminoAcids {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for AminoAcidsView<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for AminoAcidsViewMut<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.as_ref()
    }
}

impl SequenceReadable for FastQ {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastQView<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastQViewMut<'_> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastaSeq {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastaAA {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastaNT {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl SequenceReadable for FastaNTAnnot {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> SequenceReadable for LocalProfiles<'_, M, N, O, S> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        &self.seq
    }
}

impl<const M: usize, const N: usize, const O: usize, const S: usize> SequenceReadable for SharedProfiles<'_, M, N, O, S> {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        &self.seq
    }
}

/// A trait providing access to annotations within a header, as well as the ID
/// for the header (the text appearing before the first annotation).
pub trait GetAnnotation: HeaderReadable {
    /// Returns an iterator of the annotations present within the header, from
    /// left to right. Annotations are text surrounded by braces `{}`.
    ///
    /// Empty annotations are included as empty strings. Text not enclosed in
    /// braces between annotations is excluded.
    ///
    /// ## Errors
    ///
    /// The opening braces and closing braces must correspond to each other and
    /// not be nested.
    fn get_annotations(&self) -> std::io::Result<AnnotationIter<'_>> {
        self.split_annotations().map(|v| v.1)
    }

    /// Splits a header into the ID (before the first annotation) and an
    /// iterator of the annotations, from left to right. Annotations are text
    /// surrounded by braces `{}`.
    ///
    /// Empty annotations are included as empty strings. Text not enclosed in
    /// braces between annotations is excluded. If there are no annotations, the
    /// entire header is returned as the ID.
    ///
    /// ## Errors
    ///
    /// The opening braces and closing braces must correspond to each other and
    /// not be nested.
    fn split_annotations(&self) -> std::io::Result<(&str, AnnotationIter<'_>)> {
        let header = self.header();

        // Find leftmost brace (beginning of annotations)
        let Some(leftmost_brace_idx) = header.find(['}', '{']) else {
            // No annotation is present: return empty iterator
            return Ok((header, AnnotationIter { inner: None }));
        };

        // Confirm the leftmost brace is {, or error
        if header[leftmost_brace_idx..].starts_with('}') {
            return Err(std::io::Error::other(
                "Found a closing brace before an opening brace within the header",
            ));
        }

        // Find rightmost brace (end of annotations)
        let Some(rightmost_brace_idx) = header.rfind(['}', '{']) else {
            return Err(std::io::Error::other(
                "Found an opening brace with no closing brace in the header",
            ));
        };

        // Confirm the rightmost brace is }, or error
        if header[rightmost_brace_idx..].starts_with('{') {
            return Err(std::io::Error::other(
                "Found an opening brace with no closing brace in the header",
            ));
        }

        // Extract annotations
        let annotations = &header[(Bound::Excluded(leftmost_brace_idx), Bound::Excluded(rightmost_brace_idx))];

        // Extract ID
        let id = &header[..leftmost_brace_idx];

        Ok((
            id,
            AnnotationIter {
                inner: Some(annotations),
            },
        ))
    }
}

impl<T: HeaderReadable> GetAnnotation for T {}

/// A trait providing access to annotations within a header, with the ability to
/// mutate the underlying header to include only the ID appearing before the
/// first annotation.
pub trait GetAnnotationOwned: HeaderMutable {
    /// Truncates a header to exclude the annotations, returning an iterator of
    /// the annotations from left to right. Annotations are text surrounded by
    /// braces `{}`.
    ///
    /// Empty annotations are included as empty strings. Text not enclosed in
    /// braces between annotations is excluded. If there are no annotations, no
    /// truncation occurs.
    ///
    /// ## Limitations
    ///
    /// This requires an allocation for each annotation. Consider using
    /// [`split_annotations`] if references are sufficient.
    ///
    /// ## Errors
    ///
    /// The opening braces and closing braces must correspond to each other and
    /// not be nested. Assuming the leftmost brace is opening and the rightmost
    /// brace is closing, any other errors appear when using the iterator. The
    /// header is truncated even if the iterator returns an error.
    ///
    /// [`split_annotations`]: GetAnnotation::split_annotations
    fn split_off_annotations(&mut self) -> std::io::Result<AnnotationIntoIter> {
        let header = self.header_mut();

        // Find leftmost brace (beginning of annotations)
        let Some(leftmost_brace_idx) = header.find(['}', '{']) else {
            // No annotation is present: No mutation and return empty iterator
            return Ok(AnnotationIntoIter { inner: None });
        };

        // Confirm the leftmost brace is {, or error
        if header[leftmost_brace_idx..].starts_with('}') {
            return Err(std::io::Error::other(
                "Found a closing brace before an opening brace within the header",
            ));
        }

        // Find rightmost brace (end of annotations)
        let Some(rightmost_brace_idx) = header.rfind(['}', '{']) else {
            return Err(std::io::Error::other(
                "Found an opening brace with no closing brace in the header",
            ));
        };

        // Confirm the rightmost brace is }, or error
        if header[rightmost_brace_idx..].starts_with('{') {
            return Err(std::io::Error::other(
                "Found an opening brace with no closing brace in the header",
            ));
        }

        // Extract annotations
        let annotations = header[(Bound::Excluded(leftmost_brace_idx), Bound::Excluded(rightmost_brace_idx))].to_string();

        // Truncate header to before leftmost brace
        header.truncate(leftmost_brace_idx);

        Ok(AnnotationIntoIter {
            inner: Some(annotations),
        })
    }
}

impl<T: HeaderMutable> GetAnnotationOwned for T {}

/// An iterator over annotations within a header, as string slices.
pub struct AnnotationIter<'a> {
    /// The slice of the header being processed, excluding the leftmost opening
    /// brace and the rightmost closing brace.
    ///
    /// This is `None` is there are no more annotations. This is `Some("")` if
    /// there is one final annotation which happens to be empty.
    inner: Option<&'a str>,
}

impl<'a> Iterator for AnnotationIter<'a> {
    type Item = std::io::Result<&'a str>;

    fn next(&mut self) -> Option<Self::Item> {
        // Take the inner string (clearing the field to None) or abort. This
        // ensures errors abort the iterator
        let inner = std::mem::take(&mut self.inner)?;

        let Some(closing_brace_idx) = inner.find(['}', '{']) else {
            // Happy path: a single annotation and only one set of braces. This
            // will abort the iterator since self.inner is now None
            return Some(Ok(inner));
        };

        // Confirm the brace is }, or error
        if inner[closing_brace_idx..].starts_with('{') {
            return Some(Err(std::io::Error::other(
                "Found a second opening brace before a closing brace was found in the header",
            )));
        }

        let Some(opening_brace_idx) = inner.str_search_in(closing_brace_idx + 1..).find(['}', '{']) else {
            return Some(Err(std::io::Error::other(
                "Found a closing brace without a preceding opening brace in the header",
            )));
        };

        // Confirm the brace is {, or error
        if inner[opening_brace_idx..].starts_with('}') {
            return Some(Err(std::io::Error::other(
                "Found a closing brace without a preceding opening brace in the header",
            )));
        }

        // Reset inner so that the iterator continues
        self.inner = Some(&inner[opening_brace_idx + 1..]);

        Some(Ok(&inner[..closing_brace_idx]))
    }
}

/// An iterator over annotations within a header, as owned strings.
pub struct AnnotationIntoIter {
    /// The slice of the header being processed, excluding the leftmost opening
    /// brace and the rightmost closing brace.
    ///
    /// This is `None` is there are no more annotations. This is `Some("")` if
    /// there is one final annotation which happens to be empty. The opening and
    /// closing braces are excluded for efficiency, especially when there is
    /// only a single annotation. Otherwise, two allocations would likely be
    /// needed (one for the String with the braces, and one without which gets
    /// yielded).
    inner: Option<String>,
}

impl Iterator for AnnotationIntoIter {
    type Item = std::io::Result<String>;

    fn next(&mut self) -> Option<Self::Item> {
        // Take the inner string (clearing the field to None) or abort. This
        // ensures errors abort the iterator
        let mut inner = std::mem::take(&mut self.inner)?;

        let Some(closing_brace_idx) = inner.find(['}', '{']) else {
            // Happy path: a single annotation and only one set of braces. The
            // bare minimum number of allocations are performed. This will abort
            // the iterator since self.inner is now None
            return Some(Ok(inner));
        };

        // Confirm the brace is }, or error
        if inner[closing_brace_idx..].starts_with('{') {
            return Some(Err(std::io::Error::other(
                "Found a second opening brace before a closing brace was found in the header",
            )));
        }

        let Some(opening_brace_idx) = inner.str_search_in(closing_brace_idx + 1..).find(['}', '{']) else {
            return Some(Err(std::io::Error::other(
                "Found a closing brace without a preceding opening brace in the header",
            )));
        };

        // Confirm the brace is {, or error
        if inner[opening_brace_idx..].starts_with('}') {
            return Some(Err(std::io::Error::other(
                "Found a closing brace without a preceding opening brace in the header",
            )));
        }

        let annotation = inner[..closing_brace_idx].to_string();

        // Reset inner so that the iterator continues. Use drain to avoid
        // another allocation
        inner.drain(..=opening_brace_idx);
        self.inner = Some(inner);

        Some(Ok(annotation))
    }
}
