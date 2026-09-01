//! BAM binning support for encoded alignment records.

use crate::data::records::bam::error::BamRecordError;

/// Reserved bin for unmapped reads that have no reference coordinate.
const UNPLACED_UNMAPPED_BIN: u16 = 4680;
/// BAI coordinate limit for indexable files
const BAI_COORDINATE_LIMIT: u32 = 1 << 29;

/// One level in the BAM/BAI bin hierarchy.
struct BinLevel {
    /// Shifting a coordinate by this amount gives the coordinate's bin index at
    /// that level.
    shift:     u32,
    /// The first global BAM bin number assigned to this level. It is the count
    /// of all bins in coarser levels.
    first_bin: u32,
}

/// The root bin covers the full classic BAI coordinate span.
const ROOT_BIN: u16 = 0;

/// BAM/BAI bin levels in a six-level, 8-ary bin tree.
///
/// Each coarser level is 3 bits wider than the previous one. The levels are
/// ordered finest-to-coarsest so `reg2bin` returns the smallest bin that fully
/// contains the alignment interval.
const BIN_LEVELS: [BinLevel; 5] = [
    BinLevel {
        shift:     14,
        first_bin: 4681,
    },
    BinLevel {
        shift:     17,
        first_bin: 585,
    },
    BinLevel {
        shift:     20,
        first_bin: 73,
    },
    BinLevel {
        shift:     23,
        first_bin: 9,
    },
    BinLevel {
        shift:     26,
        first_bin: 1,
    },
];

/// Computes the BAI-compatible bin for an alignment.
///
/// `pos0` is the 0-based start position and `ref_span` is the number of
/// reference bases consumed by the CIGAR. Records without a coordinate (`pos0`
/// as `-1`) use BAM's reserved unplaced-unmapped bin `4680`. Records with the
/// unmapped flag set but with a coordinate are binned as a one-base interval at
/// `pos0`.
///
/// Coordinate-bearing records must have a 0-based half-open interval within
/// `[0, 2^29)`.
pub(super) fn compute_bin(pos0: i32, ref_span: u32, flag: u16) -> Result<u16, BamRecordError> {
    match pos0 {
        ..=-2 => Err(BamRecordError::BinningOutOfRange),
        -1 => Ok(UNPLACED_UNMAPPED_BIN),
        pos0 @ 0..=i32::MAX => {
            let beg = pos0.cast_unsigned();
            let effective_ref_span = if (flag & 0x4) != 0 || ref_span == 0 { 1 } else { ref_span };

            let end = beg.checked_add(effective_ref_span).ok_or(BamRecordError::BinningOutOfRange)?;

            if end > BAI_COORDINATE_LIMIT {
                return Err(BamRecordError::BinningOutOfRange);
            }
            reg2bin(beg, end)
        }
    }
}

/// Maps a half-open 0-based reference interval to the BAM/BAI bin hierarchy.
///
/// `beg` is inclusive and `end` is exclusive. The BAM/BAI bin calculation tests
/// whether both covered endpoints fall in the same bin.
///
/// The returned bin is the finest hierarchy level whose bin fully contains the
/// interval. If no non-root level contains both endpoints, the interval belongs
/// to the root bin.
fn reg2bin(beg: u32, end: u32) -> Result<u16, BamRecordError> {
    let last = end.checked_sub(1).ok_or(BamRecordError::BinningOutOfRange)?;

    for level in BIN_LEVELS {
        if (beg >> level.shift) == (last >> level.shift) {
            let bin = level.first_bin + (beg >> level.shift);

            return u16::try_from(bin).map_err(|_| BamRecordError::BinningOutOfRange);
        }
    }

    Ok(ROOT_BIN)
}
