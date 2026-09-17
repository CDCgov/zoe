//! BAM binning support for encoded alignment records.

use crate::data::records::bam::error::BamRecordError;
use std::ops::Range;

/// Reserved bin for unmapped reads that have no reference coordinate.
const UNPLACED_UNMAPPED_BIN: u16 = 4680;
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

/// Computes the BAI-compatible bin for a checked reference interval.
///
/// Records without a coordinate (`None`) use BAM's reserved unplaced-unmapped
/// bin `4680`. Coordinate-bearing records must have an interval within `[0,
/// 2^29)`.
pub(super) fn compute_bin(interval: Option<&Range<u32>>) -> Result<u16, BamRecordError> {
    let Some(interval) = interval else {
        return Ok(UNPLACED_UNMAPPED_BIN);
    };

    let last = interval.end.checked_sub(1).ok_or(BamRecordError::BinningOutOfRange)?;

    for level in BIN_LEVELS {
        if (interval.start >> level.shift) == (last >> level.shift) {
            let bin = level.first_bin + (interval.start >> level.shift);

            return u16::try_from(bin).map_err(|_| BamRecordError::BinningOutOfRange);
        }
    }

    Ok(ROOT_BIN)
}
