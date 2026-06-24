use crate::{data::sam::SamData, private::Sealed};
use std::{cmp::Ordering, collections::HashMap};

/// An extension trait adding coordinate sorting methods to collections of
/// [`SamData`].
///
/// Coordinate sorting uses the [`SamData::rname`] field as the major sort key,
/// with reference order defined by the order of `@SQ` lines in the supplied SAM
/// headers. The [`SamData::pos`] field is used as the minor sort key. *Zoe*
/// additionally sorts records with equal `RNAME` and `POS` by
/// [`SamData::qname`] so the output order is deterministic.
///
/// If an alignment has a non-`*` `RNAME` that is not listed in the supplied
/// `@SQ` headers, it is sorted after references that are present in the headers
/// but before alignments whose `RNAME` is `*`.
pub trait SamDataSort: Sealed {
    /// Sorts SAM records in place by coordinate.
    ///
    /// The `headers` argument may be any iterable of raw SAM header lines, such
    /// as the strings stored in the [`Header`] variant. Only `@SQ` lines and
    /// their `SN` fields are used to determine reference order.
    ///
    /// [`Header`]: crate::data::sam::SamRow::Header
    fn coordinate_sort<I, H>(&mut self, headers: I)
    where
        I: IntoIterator<Item = H>,
        H: AsRef<str>;
}

impl SamDataSort for [SamData] {
    #[inline]
    fn coordinate_sort<I, H>(&mut self, headers: I)
    where
        I: IntoIterator<Item = H>,
        H: AsRef<str>, {
        let reference_order = headers
            .into_iter()
            .filter_map(|header| sq_sequence_name(header.as_ref()).map(str::to_owned))
            .enumerate()
            .map(|(index, rname)| (rname, index))
            .collect();

        self.sort_by(|a, b| compare_sam_data(a, b, &reference_order));
    }
}

/// Extracts the `SN` field from an `@SQ` header line.
fn sq_sequence_name(header: &str) -> Option<&str> {
    let mut fields = header.split('\t');

    if fields.next()? != "@SQ" {
        return None;
    }

    fields.find_map(|field| field.strip_prefix("SN:"))
}

/// Compares two SAM records by coordinate sort order.
fn compare_sam_data(a: &SamData, b: &SamData, reference_order: &HashMap<String, usize>) -> Ordering {
    let a_rname = rname_sort_key(a.rname.as_str(), reference_order);
    let b_rname = rname_sort_key(b.rname.as_str(), reference_order);

    a_rname
        .cmp(&b_rname)
        .then_with(|| a.pos.cmp(&b.pos))
        .then_with(|| a.qname.cmp(&b.qname))
}

/// Returns the sortable bucket and order for an `RNAME`.
///
/// Known references sort by `@SQ` order, unknown named references follow, and
/// `*` sorts last.
fn rname_sort_key(rname: &str, reference_order: &HashMap<String, usize>) -> (usize, usize) {
    if rname == "*" {
        (2, 0)
    } else if let Some(&ref_idx) = reference_order.get(rname) {
        (0, ref_idx)
    } else {
        (1, 0)
    }
}
