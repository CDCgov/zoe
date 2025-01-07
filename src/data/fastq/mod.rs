use crate::prelude::*;

mod reader;
mod std_traits;
mod view_traits;

pub use reader::*;

/// Holds data for a single [Read](https://en.wikipedia.org/wiki/Read_(biology))
/// or [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) record.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastQ {
    pub header:   String,
    pub sequence: Nucleotides,
    pub quality:  QualityScores,
}

/// The corresponding immutable view type for [`FastQ`]. See
/// [Views](crate::data#header) for more details.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastQView<'a> {
    pub header:   &'a str,
    pub sequence: NucleotidesView<'a>,
    pub quality:  QualityScoresView<'a>,
}

/// The corresponding mutable view type for [`FastQ`]. See
/// [Views](crate::data#header) for more details.
#[derive(Eq, PartialEq, Hash, Debug)]
pub struct FastQViewMut<'a> {
    pub header:   &'a mut String,
    pub sequence: NucleotidesViewMut<'a>,
    pub quality:  QualityScoresViewMut<'a>,
}

impl FastQ {
    // Conversions and indexing

    /// Creates a new [`FastQ`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        FastQ {
            header:   String::new(),
            sequence: Nucleotides::new(),
            quality:  QualityScores::new(),
        }
    }

    // FastQ-specific methods

    /// Returns the reverse complement of the record.
    #[inline]
    #[must_use]
    pub fn to_reverse_complement(&self) -> FastQ {
        FastQ {
            header:   self.header.clone(),
            sequence: self.sequence.to_reverse_complement(),
            quality:  self.quality.to_reverse(),
        }
    }

    /// Computes the reverse complement of the record in-place.
    #[inline]
    pub fn make_reverse_complement(&mut self) {
        self.sequence.make_reverse_complement();
        self.quality.make_reverse();
    }
}

impl FastQView<'_> {
    // Conversions and indexing

    /// Creates a new [`FastQView`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        FastQView {
            header:   "",
            sequence: NucleotidesView::new(),
            quality:  QualityScoresView::new(),
        }
    }

    // FastQ-specific methods

    /// Returns the reverse complement of the record.
    #[inline]
    #[must_use]
    pub fn to_reverse_complement(&self) -> FastQ {
        FastQ {
            header:   self.header.to_string(),
            sequence: self.sequence.to_reverse_complement(),
            quality:  self.quality.to_reverse(),
        }
    }
}

impl FastQViewMut<'_> {
    // No new() exists for FastQViewMut because it is not possible to create an
    // empty &mut String

    // FastQ-specific methods

    /// Returns the reverse complement of the record.
    #[inline]
    #[must_use]
    pub fn to_reverse_complement(&self) -> FastQ {
        FastQ {
            header:   (*self.header).to_string(),
            sequence: self.sequence.to_reverse_complement(),
            quality:  self.quality.to_reverse(),
        }
    }

    /// Computes the reverse complement of the record in-place.
    #[inline]
    pub fn make_reverse_complement(&mut self) {
        self.sequence.make_reverse_complement();
        self.quality.make_reverse();
    }
}

#[cfg(test)]
mod test;
