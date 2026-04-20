use crate::{
    data::{
        id_types::FastaIDs,
        types::{
            amino_acids::AminoAcids,
            nucleotides::{Nucleotides, ToDNA, Translate},
        },
        views::ViewAssocTypes,
    },
    prelude::{AminoAcidsView, DataOwned, NucleotidesView},
};

mod reader;
mod std_traits;
mod taxon;
mod view_traits;

pub use reader::*;
pub use taxon::*;

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
    S: ViewAssocTypes, {
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
/// ## Taxon Annotation
///
/// Unlike [`Fasta`] and [`FastaView`], [`FastaViewMut`] does not have a
/// `annotate_taxon` method due to ownership semantics. We recommend annotating
/// an owned [`Fasta`] with [`annotate_taxon`] and then calling [`as_view_mut`].
///
/// [`NucleotidesViewMut`]: crate::data::types::nucleotides::NucleotidesViewMut
/// [`AminoAcidsViewMut`]: crate::data::types::amino_acids::AminoAcidsViewMut
/// [`annotate_taxon`]: Fasta::annotate_taxon
/// [`as_view_mut`]: DataOwned::as_view_mut
#[derive(Eq, PartialEq, Hash, Debug)]
pub struct FastaViewMut<'a, S = Vec<u8>>
where
    S: ViewAssocTypes, {
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
    /// Constructs an annotated record from the provided metadata.
    #[inline]
    #[must_use]
    pub fn with_annot<M>(self, annot: M) -> FastaAnnot<M, S> {
        FastaAnnot {
            header: self.header,
            sequence: self.sequence,
            annot,
        }
    }

    /// Splits the taxon off of the header and include it as an annotation in a
    /// [`FastaAnnot`] object.
    ///
    /// ## Limitations
    ///
    /// This will involve two copies of strings (one for the header, and one for
    /// the taxon). Consider calling [`as_view`] first to avoid this.
    ///
    /// ## Errors
    ///
    /// If no taxon is found, the original input is returned as an `Err` for
    /// proper handling.
    ///
    /// [`as_view`]: DataOwned::as_view
    pub fn annotate_taxon(self) -> Result<FastaAnnot<Taxon, S>, Fasta<S>> {
        if let Some((id, taxon)) = self.header.get_id_taxon() {
            Ok(FastaAnnot {
                header:   id.to_string(),
                sequence: self.sequence,
                annot:    Taxon::from(taxon),
            })
        } else {
            Err(self)
        }
    }

    /// Extracts the ID and taxon from the header, if the header is of the
    /// format `id{taxon}`.
    ///
    /// The braces are not included in the outputs. Anything after the closing
    /// brace is ignored.
    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, TaxonView<'_>)> {
        let (id, taxon) = self.header.get_id_taxon()?;
        Some((id, TaxonView::from(taxon)))
    }
}

impl<'a, S> FastaView<'a, S>
where
    S: DataOwned,
{
    /// Constructs an annotated record from the provided metadata.
    #[inline]
    #[must_use]
    pub fn with_annot<M>(self, annot: M::View<'a>) -> FastaAnnotView<'a, M, S>
    where
        M: DataOwned, {
        FastaAnnotView {
            header: self.header,
            sequence: self.sequence,
            annot,
        }
    }

    /// Splits the taxon off of the header and include it as an annotation in a
    /// [`FastaAnnot`] object.
    ///
    /// ## Errors
    ///
    /// If no taxon is found, the original input is returned as an `Err` for
    /// proper handling.
    pub fn annotate_taxon(self) -> Result<FastaAnnotView<'a, Taxon, S>, FastaView<'a, S>> {
        if let Some((id, taxon)) = self.header.get_id_taxon() {
            Ok(FastaAnnotView {
                header:   id,
                sequence: self.sequence,
                annot:    TaxonView::from(taxon),
            })
        } else {
            Err(self)
        }
    }

    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, TaxonView<'_>)> {
        let (id, taxon) = self.header.get_id_taxon()?;
        Some((id, TaxonView::from(taxon)))
    }
}

impl<'a, S> FastaViewMut<'a, S>
where
    S: DataOwned,
{
    /// Constructs an annotated record from the provided metadata.
    #[inline]
    #[must_use]
    pub fn with_annot<M>(self, annot: M::ViewMut<'a>) -> FastaAnnotViewMut<'a, M, S>
    where
        M: DataOwned, {
        FastaAnnotViewMut {
            header: self.header,
            sequence: self.sequence,
            annot,
        }
    }

    #[inline]
    #[must_use]
    pub fn get_id_taxon(&self) -> Option<(&str, TaxonView<'_>)> {
        let (id, taxon) = self.header.get_id_taxon()?;
        Some((id, TaxonView::from(taxon)))
    }
}

/// A [FASTA](https://en.wikipedia.org/wiki/FASTA_format) record containing a
/// header and a sequence, along with an annotation of type `M`.
///
/// ## Parameters
///
/// - `M`: The type for the annotation/metadata
/// - `S`: The type of the sequence, which defaults to `Vec<u8>` but can also be
///   [`Nucleotides`] or [`AminoAcids`]
#[derive(Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct FastaAnnot<M, S = Vec<u8>> {
    pub header:   String,
    pub sequence: S,
    pub annot:    M,
}

/// The corresponding immutable view type for [`FastaAnnot`].
///
/// See [Views](crate::data#views) for more details.
///
/// ## Parameters
///
/// - `M`: The owned type corresponding to the annotation/metadata
/// - `S`: The owned type corresponding to the sequence, which defaults to
///   `Vec<u8>` but can also be [`Nucleotides`] or [`AminoAcids`]
#[derive(Eq, PartialEq, Hash, Debug)]
pub struct FastaAnnotView<'a, M, S = Vec<u8>>
where
    S: ViewAssocTypes,
    M: ViewAssocTypes, {
    pub header:   &'a str,
    pub sequence: S::View<'a>,
    pub annot:    M::View<'a>,
}

/// The corresponding mutable view type for [`FastaAnnot`].
///
/// See [Views](crate::data#views) for more details.
///
/// ## Parameters
///
/// - `M`: The owned type corresponding to the annotation/metadata
/// - `S`: The owned type corresponding to the sequence, which defaults to
///   `Vec<u8>` but can also be [`Nucleotides`] or [`AminoAcids`]
///
/// [`NucleotidesViewMut`]: crate::data::types::nucleotides::NucleotidesViewMut
/// [`AminoAcidsViewMut`]: crate::data::types::amino_acids::AminoAcidsViewMut
#[derive(Eq, PartialEq, Hash, Debug)]
pub struct FastaAnnotViewMut<'a, M, S = Vec<u8>>
where
    S: ViewAssocTypes,
    M: ViewAssocTypes, {
    pub header:   &'a mut String,
    pub sequence: S::ViewMut<'a>,
    pub annot:    M::ViewMut<'a>,
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

impl<M, S> FastaAnnot<M, S> {
    /// Transforms the header in a [`FastaAnnot`] record using a closure `f`.
    #[inline]
    #[must_use]
    pub fn map_header<F>(self, f: F) -> Self
    where
        F: FnOnce(String) -> String, {
        Self {
            header:   f(self.header),
            sequence: self.sequence,
            annot:    self.annot,
        }
    }

    /// Transforms the sequence in a [`FastaAnnot`] record using a closure `f`.
    ///
    /// This may change the type of the sequence.
    #[inline]
    #[must_use]
    pub fn map_sequence<U, F>(self, f: F) -> FastaAnnot<M, U>
    where
        F: FnOnce(S) -> U, {
        FastaAnnot {
            header:   self.header,
            sequence: f(self.sequence),
            annot:    self.annot,
        }
    }

    /// Transforms the annotation in a [`FastaAnnot`] record using a closure
    /// `f`.
    ///
    /// This may change the type of the annotation.
    #[inline]
    #[must_use]
    pub fn map_annot<U, F>(self, f: F) -> FastaAnnot<U, S>
    where
        F: FnOnce(M) -> U, {
        FastaAnnot {
            header:   self.header,
            sequence: self.sequence,
            annot:    f(self.annot),
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

impl<M> FastaAnnot<M, Vec<u8>> {
    /// Recodes to uppercase IUPAC DNA with corrected gaps, otherwise mapping to
    /// `N`. Returns the sequence as [`Nucleotides`].
    #[inline]
    #[must_use]
    pub fn recode_to_dna(self) -> FastaAnnot<M, Nucleotides> {
        self.map_sequence(ToDNA::recode_to_dna)
    }

    /// Filters and recodes to uppercase IUPAC DNA with corrected gaps. Returns
    /// the sequence as [`Nucleotides`].
    #[inline]
    #[must_use]
    pub fn filter_to_dna(self) -> FastaAnnot<M, Nucleotides> {
        self.map_sequence(ToDNA::filter_to_dna)
    }

    /// Filters and recodes to uppercase IUPAC DNA without gaps. Returns
    /// the sequence as [`Nucleotides`].
    #[inline]
    #[must_use]
    pub fn filter_to_dna_unaligned(self) -> FastaAnnot<M, Nucleotides> {
        self.map_sequence(ToDNA::filter_to_dna_unaligned)
    }

    /// Converts the sequence type to [`Nucleotides`] without any checking.
    #[inline]
    #[must_use]
    pub fn into_dna(self) -> FastaAnnot<M, Nucleotides> {
        self.map_sequence(Into::into)
    }

    /// Converts the sequence type to [`AminoAcids`] without any checking.
    #[inline]
    #[must_use]
    pub fn into_aa(self) -> FastaAnnot<M, AminoAcids> {
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

impl<M> FastaAnnot<M, Nucleotides> {
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
    pub fn translate(self) -> FastaAnnot<M, AminoAcids> {
        self.map_sequence(|x| Translate::translate(&x))
    }
}
