use crate::alignment::phmm::state::PhmmBacktrackFlags;

mod global;

/// A traceback matrix storing only a diagonal band of the Viterbi DP table.
///
/// Rows correspond to pHMM layers and columns correspond to the number of
/// query residues consumed. Each row reserves `2 * band_width + 1` entries,
/// including unused entries where the band extends beyond the query bounds.
struct BandedViterbiTraceback {
    data:       Vec<PhmmBacktrackFlags>,
    band_width: usize,
}

impl BandedViterbiTraceback {
    #[inline]
    #[must_use]
    fn new(rows: usize, band_width: usize) -> Self {
        let band_full_width = 2 * band_width + 1;
        Self {
            data: vec![PhmmBacktrackFlags::new(); rows * band_full_width],
            band_width,
        }
    }

    /// Returns the flattened index for the DP coordinate `(j, i)`.
    #[inline]
    #[must_use]
    fn index(&self, i: usize, j: usize) -> usize {
        let band_full_width = 2 * self.band_width + 1;
        let num_cols_skipped = j.saturating_sub(self.band_width);
        debug_assert!(i >= num_cols_skipped);
        debug_assert!(i - num_cols_skipped < band_full_width);
        j * band_full_width + i - num_cols_skipped
    }

    #[inline]
    #[must_use]
    fn get(&self, i: usize, j: usize) -> PhmmBacktrackFlags {
        self.data[self.index(i, j)]
    }

    #[inline]
    #[must_use]
    fn get_mut(&mut self, i: usize, j: usize) -> &mut PhmmBacktrackFlags {
        let index = self.index(i, j);
        &mut self.data[index]
    }
}
