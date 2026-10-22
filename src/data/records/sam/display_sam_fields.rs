use crate::data::sam::GetSamFields;
use std::fmt::{Display, Write};

/// An extension of [`GetSamFields`] that provides a
/// [`display`](DisplaySamFields::display) method for automatically displaying
/// the struct in the SAM format, with fields separated by tabs.
///
/// For this to be implemented, [`GetSamFields`] must be implemented with all
/// associated types implementing [`Display`].
pub trait DisplaySamFields:
    GetSamFields<
        Qname: Display,
        Flag: Display,
        Rname: Display,
        Pos: Display,
        Mapq: Display,
        Cigar: Display,
        Seq: Display,
        Qual: Display,
        OptFields: Display,
    > {
    /// Returns a displayable representation of the data in SAM format.
    ///
    /// See [`SamFieldsDisplay`] for more details.
    fn display(&self) -> SamFieldsDisplay<&Self> {
        SamFieldsDisplay(self)
    }
}

impl<R: GetSamFields> DisplaySamFields for R
where
    R::Qname: Display,
    R::Flag: Display,
    R::Rname: Display,
    R::Pos: Display,
    R::Mapq: Display,
    R::Cigar: Display,
    R::Seq: Display,
    R::Qual: Display,
    R::OptFields: Display,
{
}

/// A display representation for [`DisplaySamFields`] in the SAM file format.
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

impl<T: DisplaySamFields> Display for SamFieldsDisplay<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let qname = self.0.qname().display_or_star();
        let flag = self.0.flag().display_or("0");
        let rname = self.0.rname().display_or_star();
        let pos = self.0.pos().display_or("0");
        let mapq = self.0.mapq().display_or("255");
        let cigar = self.0.cigar().display_or_star();
        let rnext = '*';
        let pnext = 0;
        let tlen = 0;
        let seq = self.0.seq().display_or_star();
        let qual = self.0.qual().display_or_star();
        let opt_fields = self.0.opt_fields().display_leading_tab();

        // opt_fields adds its own tab if present
        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}{opt_fields}"
        )
    }
}

/// An extension trait for displaying data in an [`Option`].
trait DisplayOrExt: Sized {
    /// The inner type within the [`Option`].
    type Inner;

    /// Returns a type displaying the contained value if `Some` or a default
    /// static string if `None`.
    fn display_or(self, default: &'static str) -> DisplayOr<Self::Inner>;

    /// Returns a type displaying the contained value if `Some`, or `*` if
    /// `None` or the display implementation of the contained value displays no
    /// bytes.
    fn display_or_star(self) -> DisplayOrStar<Self::Inner>;

    /// Returns a type displaying the contained value with a leading tab if
    /// `Some`, or displaying nothing if `None`.
    fn display_leading_tab(self) -> DisplayLeadingTab<Self::Inner>;
}

impl<T> DisplayOrExt for Option<T> {
    type Inner = T;

    fn display_or(self, default: &'static str) -> DisplayOr<Self::Inner> {
        DisplayOr { value: self, default }
    }

    fn display_or_star(self) -> DisplayOrStar<Self::Inner> {
        DisplayOrStar { value: self }
    }

    fn display_leading_tab(self) -> DisplayLeadingTab<Self::Inner> {
        DisplayLeadingTab { value: self }
    }
}

/// A wraper type for displaying the contained value in an [`Option`] or a
/// default static string if `None`.
struct DisplayOr<T> {
    value:   Option<T>,
    default: &'static str,
}

impl<T: Display> Display for DisplayOr<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.value {
            Some(value) => write!(f, "{value}"),
            None => f.write_str(self.default),
        }
    }
}

/// A wrapper type for displaying the contained value in an [`Option`], or a
/// star if `None` or the display implementation does not write any bytes.
struct DisplayOrStar<T> {
    value: Option<T>,
}

impl<T: Display> Display for DisplayOrStar<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Write star if None
        let Some(value) = &self.value else {
            return f.write_char('*');
        };

        let mut f = TrackedFmt {
            f,
            anything_written: false,
        };

        write!(&mut f, "{value}")?;

        let TrackedFmt { f, anything_written } = f;

        // Write star if nothing written
        if anything_written { Ok(()) } else { f.write_char('*') }
    }
}

/// A wrapper type for displaying the contained value with a leading tab if
/// `Some`, or displaying nothing if `None`.
struct DisplayLeadingTab<T> {
    value: Option<T>,
}

impl<T: Display> Display for DisplayLeadingTab<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Some(value) = &self.value else { return Ok(()) };

        let mut f = LeadingTabFmt { f, tab_written: false };
        write!(f, "{value}")
    }
}

/// A helper formatter that tracks whether anything was written, while also
/// writing to the wrapped formatter as well.
struct TrackedFmt<'a, 'b> {
    f:                &'a mut std::fmt::Formatter<'b>,
    anything_written: bool,
}

impl std::fmt::Write for TrackedFmt<'_, '_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        if !s.is_empty() {
            self.anything_written = true;
        }
        self.f.write_str(s)
    }

    fn write_char(&mut self, c: char) -> std::fmt::Result {
        self.anything_written = true;
        self.f.write_char(c)
    }
}

/// A helper formatter that displays a leading tab before the first written
/// byte, but excludes it if no bytes are written.
struct LeadingTabFmt<'a, 'b> {
    f:           &'a mut std::fmt::Formatter<'b>,
    tab_written: bool,
}

impl std::fmt::Write for LeadingTabFmt<'_, '_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        if !self.tab_written && !s.is_empty() {
            self.tab_written = true;
            self.f.write_char('\t')?;
        }

        self.f.write_str(s)
    }

    fn write_char(&mut self, c: char) -> std::fmt::Result {
        if !self.tab_written {
            self.tab_written = true;
            self.f.write_char('\t')?;
        }

        self.f.write_char(c)
    }
}
