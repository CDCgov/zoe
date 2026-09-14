use crate::data::{
    cigar::ToCigletIterator,
    nucleotides::NucleotidesView,
    phred::QualityScoresView,
    sam::{SamData, SamDataView, ToOptFieldsIterator, is_missing_sam_field},
    views::AsView,
};
use std::fmt::Display;

/// A trait providing getters for SAM-like data allowing it to be used in
/// generic contexts.
///
/// This trait provides individual getters for each field, in addition to
/// allowing flexible types to be returned for CIGAR and the optional fields. It
/// also provides a [`display`] method for displaying the data in SAM format.
///
/// All the getters return `Option`, which is to make it easier to implement the
/// trait for custom structs that do not hold some SAM fields. Returning `Some`
/// does not mean the field is present though. For example, `Some("*")` could be
/// returned for the `cigar` field of an unmapped record.
pub trait GetSamFields: Sized {
    /// Returns the QNAME field if it is contained in the struct.
    fn qname(&self) -> Option<&str>;

    /// Returns the FLAG field if it is contained in the struct.
    fn flag(&self) -> Option<u16>;

    /// Returns the RNAME field if it is contained in the struct.
    fn rname(&self) -> Option<&str>;

    /// Returns the POS field if it is contained in the struct.
    fn pos(&self) -> Option<usize>;

    /// Returns the MAPQ field if it is contained in the struct.
    fn mapq(&self) -> Option<u8>;

    /// Returns the CIGAR field if it is present in the struct.
    ///
    /// The returned type must be able to be iterated over (via the
    /// [`ToCigletIterator`] trait) as well as displayed. Examples include
    /// [`Cigar`], [`CigarView`], and [`AlignmentStates`].
    ///
    /// [`Cigar`]: crate::data::types::cigar::Cigar
    /// [`CigarView`]: crate::data::types::cigar::CigarView
    /// [`AlignmentStates`]: crate::alignment::AlignmentStates
    fn cigar(&self) -> Option<impl ToCigletIterator + Display>;

    /// Returns the SEQ field if it is present in the struct.
    fn seq(&self) -> Option<NucleotidesView<'_>>;

    /// Returns the QUAL field if it is present in the struct.
    fn qual(&self) -> Option<QualityScoresView<'_>>;

    /// Returns any optional SAM fields stored in the struct.
    ///
    /// The returned type must be able to be iterated over (via the
    /// [`ToOptFieldsIterator`] trait) as well as displayed. Examples include
    /// [`SamOptRaw`] and [`SamOptRawView`]. Custom types can be used that
    /// implement these traits as well.
    ///
    /// [`SamOptRaw`]: crate::data::records::sam::SamOptRaw
    /// [`SamOptRawView`]: crate::data::records::sam::SamOptRawView
    fn opt_fields(&self) -> Option<impl ToOptFieldsIterator + Display>;

    /// Returns a displayable representation of the data in SAM format.
    ///
    /// See [`SamFieldsDisplay`] for more details.
    fn display(&self) -> SamFieldsDisplay<&Self> {
        SamFieldsDisplay(self)
    }
}

impl<'a, T: GetSamFields> GetSamFields for &'a T
where
    T: 'a,
{
    fn qname(&self) -> Option<&str> {
        (*self).qname()
    }

    fn flag(&self) -> Option<u16> {
        (*self).flag()
    }

    fn rname(&self) -> Option<&str> {
        (*self).rname()
    }

    fn pos(&self) -> Option<usize> {
        (*self).pos()
    }

    fn mapq(&self) -> Option<u8> {
        (*self).mapq()
    }

    fn cigar(&self) -> Option<impl ToCigletIterator + Display> {
        (*self).cigar()
    }

    fn seq(&self) -> Option<NucleotidesView<'_>> {
        (*self).seq()
    }

    fn qual(&self) -> Option<QualityScoresView<'_>> {
        (*self).qual()
    }

    fn opt_fields(&self) -> Option<impl ToOptFieldsIterator + Display> {
        (*self).opt_fields()
    }
}

impl GetSamFields for SamData {
    fn qname(&self) -> Option<&str> {
        Some(&self.qname)
    }

    fn flag(&self) -> Option<u16> {
        Some(self.flag)
    }

    fn rname(&self) -> Option<&str> {
        Some(&self.rname)
    }

    fn pos(&self) -> Option<usize> {
        Some(self.pos)
    }

    fn mapq(&self) -> Option<u8> {
        Some(self.mapq)
    }

    fn cigar(&self) -> Option<impl ToCigletIterator + Display> {
        Some(self.cigar.as_view())
    }

    fn seq(&self) -> Option<NucleotidesView<'_>> {
        Some(self.seq.as_view())
    }

    fn qual(&self) -> Option<QualityScoresView<'_>> {
        Some(self.qual.as_view())
    }

    fn opt_fields(&self) -> Option<impl ToOptFieldsIterator + Display> {
        Some(&self.opt_fields)
    }
}

impl GetSamFields for SamDataView<'_> {
    fn qname(&self) -> Option<&str> {
        Some(self.qname)
    }

    fn flag(&self) -> Option<u16> {
        Some(self.flag)
    }

    fn rname(&self) -> Option<&str> {
        Some(self.rname)
    }

    fn pos(&self) -> Option<usize> {
        Some(self.pos)
    }

    fn mapq(&self) -> Option<u8> {
        Some(self.mapq)
    }

    fn cigar(&self) -> Option<impl ToCigletIterator + Display> {
        Some(self.cigar)
    }

    fn seq(&self) -> Option<NucleotidesView<'_>> {
        Some(self.seq)
    }

    fn qual(&self) -> Option<QualityScoresView<'_>> {
        Some(self.qual)
    }

    fn opt_fields(&self) -> Option<impl ToOptFieldsIterator + Display> {
        Some(self.opt_fields)
    }
}

/// A display representation for [`GetSamFields`] in the SAM file format.
///
/// Any fields that are not stored (i.e., return `None` from their getter
/// methods) will be replaced with default values per the SAM specs:
///
/// - QNAME: `*`
/// - FLAG: 0
/// - RNAME: `*`
/// - POS: 0
/// - MAPQ: 255
/// - CIGAR: `*`
/// - RNEXT: `*`
/// - PNEXT: 0
/// - TLEN: 0
/// - SEQ: `*`
/// - QUAL: `*`
///
/// Empty QNAME, RNAME, CIGAR, SEQ, and QUAL also use the above defaults.
pub struct SamFieldsDisplay<T>(T);

impl<T: GetSamFields> Display for SamFieldsDisplay<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let qname = self.0.qname().filter(|q| !is_missing_sam_field(q)).unwrap_or("*");
        let flag = self.0.flag().unwrap_or(0);
        let rname = self.0.rname().filter(|r| !is_missing_sam_field(r)).unwrap_or("*");
        let pos = self.0.pos().unwrap_or(0);
        let mapq = self.0.mapq().unwrap_or(255);
        let cigar = CigarDisplay(self.0.cigar());
        let rnext = '*';
        let pnext = 0;
        let tlen = 0;
        let seq = self
            .0
            .seq()
            .filter(|s| !is_missing_sam_field(s))
            .unwrap_or(NucleotidesView::from(b"*"));
        // Safety: `b"*"` is graphic ASCII, which satisfies
        // `QualityScoresView`'s byte invariant.
        let qual = self
            .0
            .qual()
            .filter(|q| !is_missing_sam_field(q))
            .unwrap_or(unsafe { QualityScoresView::from_bytes_unchecked(b"*") });

        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
        )?;

        if let Some(opt_fields) = self.0.opt_fields()
            && opt_fields.to_field_iter().next().is_some()
        {
            write!(f, "\t{opt_fields}")?;
        }

        Ok(())
    }
}

/// An internal display wrapper for an optional CIGAR string, where `None` and
/// an empty iterator both display `*`.
struct CigarDisplay<C>(Option<C>);

impl<C: ToCigletIterator + Display> Display for CigarDisplay<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(cigar) = &self.0
            && cigar.to_ciglet_iterator_checked().next().is_some()
        {
            write!(f, "{cigar}")
        } else {
            write!(f, "*")
        }
    }
}
