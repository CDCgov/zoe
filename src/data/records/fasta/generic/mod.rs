use crate::{
    data::{
        types::{
            amino_acids::AminoAcids,
            nucleotides::{Nucleotides, ToDNA, Translate},
        },
        views::{AssocViewMutType, AssocViewType},
    },
    prelude::{AminoAcidsView, NucleotidesView},
};

mod reader;
mod std_traits;
mod view_traits;

pub use reader::*;

/// A [FASTA](https://en.wikipedia.org/wiki/FASTA_format) record containing a
/// header and a sequence.
///
/// ## Parameters
///
/// By default, [`Fasta`] holds `Vec<u8>` data, but it can also hold
/// [`Nucleotides`] or [`AminoAcids`] by specifying this for `S`.
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct Fasta<S = Vec<u8>> {
    pub header:   String,
    pub sequence: S,
}

/// The corresponding immutable view type for [`Fasta`].
///
/// See [Views](crate::data#views) for more details.
///
/// ## Parameters
///
/// By default, [`FastaView`] holds arbitrary `&[u8]` data, but it can also hold
/// [`NucleotidesView`] or [`AminoAcidsView`]. To specify these, set `S` to
/// [`Nucleotides`] or [`AminoAcids`] (the owned data type).
///
/// [`NucleotidesView`]: crate::data::types::nucleotides::NucleotidesView
/// [`AminoAcidsView`]: crate::data::types::amino_acids::AminoAcidsView
#[derive(Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaView<'a, S = Vec<u8>>
where
    S: AssocViewType, {
    pub header:   &'a str,
    pub sequence: S::View<'a>,
}

/// The corresponding mutable view type for [`Fasta`].
///
/// See [Views](crate::data#views) for more details.
///
/// ## Parameters
///
/// By default, [`FastaView`] holds arbitrary `&mut [u8]` data, but it can also
/// hold [`NucleotidesViewMut`] or [`AminoAcidsViewMut`]. To specify these, set
/// `S` to [`Nucleotides`] or [`AminoAcids`] (the owned data type).
///
/// [`NucleotidesViewMut`]: crate::data::types::nucleotides::NucleotidesViewMut
/// [`AminoAcidsViewMut`]: crate::data::types::amino_acids::AminoAcidsViewMut
/// [`as_view_mut`]: crate::data::views::AsViewMut::as_view_mut
#[derive(Eq, PartialEq, Hash, Debug)]
pub struct FastaViewMut<'a, S = Vec<u8>>
where
    S: AssocViewMutType, {
    pub header:   &'a mut String,
    pub sequence: S::ViewMut<'a>,
}

impl Fasta<Vec<u8>> {
    /// Creates a new [`Fasta`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            header:   String::new(),
            sequence: Vec::<u8>::new(),
        }
    }
}

impl Fasta<Nucleotides> {
    /// Creates a new [`Fasta`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            header:   String::new(),
            sequence: Nucleotides::new(),
        }
    }
}

impl Fasta<AminoAcids> {
    /// Creates a new [`Fasta`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            header:   String::new(),
            sequence: AminoAcids::new(),
        }
    }
}

impl FastaView<'_, Vec<u8>> {
    /// Creates a new [`FastaView`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            header:   "",
            sequence: &[],
        }
    }
}

impl FastaView<'_, Nucleotides> {
    /// Creates a new [`FastaView`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            header:   "",
            sequence: NucleotidesView::new(),
        }
    }
}

impl FastaView<'_, AminoAcids> {
    /// Creates a new [`FastaView`] empty object.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self {
            header:   "",
            sequence: AminoAcidsView::new(),
        }
    }
}

impl<S> Fasta<S> {
    /// Transforms the header in a [`Fasta`] record using a closure `f`.
    #[inline]
    #[must_use]
    pub fn map_header<F>(self, f: F) -> Self
    where
        F: FnOnce(String) -> String, {
        Self {
            header:   f(self.header),
            sequence: self.sequence,
        }
    }

    /// Transforms the sequence in a [`Fasta`] record using a closure `f`.
    ///
    /// This may change the type of the sequence.
    #[inline]
    #[must_use]
    pub fn map_sequence<U, F>(self, f: F) -> Fasta<U>
    where
        F: FnOnce(S) -> U, {
        Fasta {
            header:   self.header,
            sequence: f(self.sequence),
        }
    }
}

impl Fasta<Vec<u8>> {
    /// Recodes to uppercase IUPAC DNA with corrected gaps, otherwise mapping to
    /// `N`. Returns the sequence as [`Nucleotides`].
    #[inline]
    #[must_use]
    pub fn recode_to_dna(self) -> Fasta<Nucleotides> {
        self.map_sequence(ToDNA::recode_to_dna)
    }

    /// Filters and recodes to uppercase IUPAC DNA with corrected gaps. Returns
    /// the sequence as [`Nucleotides`].
    #[inline]
    #[must_use]
    pub fn filter_to_dna(self) -> Fasta<Nucleotides> {
        self.map_sequence(ToDNA::filter_to_dna)
    }

    /// Filters and recodes to uppercase IUPAC DNA without gaps. Returns
    /// the sequence as [`Nucleotides`].
    #[inline]
    #[must_use]
    pub fn filter_to_dna_unaligned(self) -> Fasta<Nucleotides> {
        self.map_sequence(ToDNA::filter_to_dna_unaligned)
    }

    /// Converts the sequence type to [`Nucleotides`] without any checking.
    #[inline]
    #[must_use]
    pub fn into_dna(self) -> Fasta<Nucleotides> {
        self.map_sequence(Into::into)
    }

    /// Converts the sequence type to [`AminoAcids`] without any checking.
    #[inline]
    #[must_use]
    pub fn into_aa(self) -> Fasta<AminoAcids> {
        self.map_sequence(Into::into)
    }
}

impl Fasta<Nucleotides> {
    /// Computes the reverse complement of the sequence in-place.
    #[inline]
    pub fn make_reverse_complement(&mut self) {
        self.sequence.make_reverse_complement();
    }

    /// Translates the DNA sequence to [`AminoAcids`].
    ///
    /// ## Limitations
    ///
    /// This uses a new buffer for the sequence.
    #[must_use]
    pub fn translate(self) -> Fasta<AminoAcids> {
        self.map_sequence(|x| Translate::translate(&x))
    }
}
