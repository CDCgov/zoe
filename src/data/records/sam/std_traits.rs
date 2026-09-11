use super::*;

impl std::fmt::Display for SamData {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.as_view().fmt(f)
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
            opt_fields,
        } = self;

        let qname = if is_missing_sam_field(qname) { "*" } else { qname };
        let rname = if is_missing_sam_field(rname) { "*" } else { rname };

        let cigar = if is_missing_sam_field(cigar.as_bytes()) {
            CigarView::from_slice_unchecked(b"*")
        } else {
            *cigar
        };

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
        )?;

        for opt_field in opt_fields.iter_raw() {
            write!(f, "\t{opt_field}")?;
        }
        Ok(())
    }
}

#[allow(deprecated)]
impl std::fmt::Display for SamDataViewMut<'_> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.as_view().fmt(f)
    }
}

impl Display for SamOptField {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}:{}", self.tag[0] as char, self.tag[1] as char, self.value)
    }
}

impl Display for SamOptValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SamOptValue::Char(val) => write!(f, "A:{val}", val = *val as char),
            SamOptValue::Int(val) => write!(f, "i:{val}"),
            SamOptValue::Float(val) => write!(f, "f:{val}"),
            SamOptValue::String(val) => write!(f, "Z:{val}"),
            SamOptValue::Hex(val) => write!(f, "H:{val}"),
            SamOptValue::Array(val) => match val {
                OptArray::I8(vec) => OptArray::fmt_opt_array(f, 'c', vec),
                OptArray::U8(vec) => OptArray::fmt_opt_array(f, 'C', vec),
                OptArray::I16(vec) => OptArray::fmt_opt_array(f, 's', vec),
                OptArray::U16(vec) => OptArray::fmt_opt_array(f, 'S', vec),
                OptArray::I32(vec) => OptArray::fmt_opt_array(f, 'i', vec),
                OptArray::U32(vec) => OptArray::fmt_opt_array(f, 'I', vec),
                OptArray::F32(vec) => OptArray::fmt_opt_array(f, 'f', vec),
            },
        }
    }
}

impl OptArray {
    /// A helper function for displaying a [`SamOptValue::Array`].
    fn fmt_opt_array<T: Display>(f: &mut Formatter<'_>, arr_type: char, vals: &[T]) -> std::fmt::Result {
        write!(f, "B:{arr_type}")?;

        for v in vals {
            write!(f, ",{v}")?;
        }

        Ok(())
    }
}
