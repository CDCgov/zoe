use super::*;

impl std::fmt::Display for SamData {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let SamData {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
            opt_fields,
        } = self;

        let seq = if is_missing_sam_field(seq) {
            NucleotidesView::from(b"*")
        } else {
            seq.as_view()
        };
        let qual = if is_missing_sam_field(qual) {
            // Safety: `b"*"` is graphic ASCII, which satisfies
            // `QualityScoresView`'s byte invariant.
            unsafe { QualityScoresView::from_bytes_unchecked(b"*") }
        } else {
            qual.as_view()
        };

        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
        )?;

        for opt_field in &opt_fields.0 {
            write!(f, "\t{opt_field}")?;
        }
        Ok(())
    }
}

impl std::fmt::Display for SamDataView<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let SamDataView {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
        } = self;

        let seq = if is_missing_sam_field(seq) {
            NucleotidesView::from(b"*")
        } else {
            *seq
        };
        let qual = if is_missing_sam_field(qual) {
            // Safety: `b"*"` is graphic ASCII, which satisfies
            // `QualityScoresView`'s byte invariant.
            unsafe { QualityScoresView::from_bytes_unchecked(b"*") }
        } else {
            *qual
        };

        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
        )
    }
}

impl std::fmt::Display for SamDataViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let SamDataViewMut {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
        } = self;

        let seq = if is_missing_sam_field(seq) {
            NucleotidesView::from(b"*")
        } else {
            seq.as_view()
        };
        let qual = if is_missing_sam_field(qual) {
            // Safety: `b"*"` is graphic ASCII, which satisfies
            // `QualityScoresView`'s byte invariant.
            unsafe { QualityScoresView::from_bytes_unchecked(b"*") }
        } else {
            qual.as_view()
        };

        write!(
            f,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
        )
    }
}
