//! Structs and methods for building an in-memory BAI index (while BAM records
//! are streamed).

use crate::data::{
    bam::error::{BamEncodingError, BamError, NumberSizeTarget},
    err::ResultWithErrorContext,
};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    ops::Range,
    path::Path,
};

/// Number of coordinate bits represented by one linear-index window.
const LINEAR_INDEX_SHIFT: u32 = 14;

/// Accumulates BAI data and writes it to a companion output stream.
pub(super) struct BaiWriter {
    inner:           BufWriter<File>,
    index:           BaiIndex,
    last_coordinate: Option<(u8, i32, i32)>,
}

impl BaiWriter {
    /// Creates an index writer for `path`.
    pub(super) fn from_path(path: impl AsRef<Path>) -> Result<Self, BamError> {
        let file = File::create(&path).with_path_context("Cannot create .bai file to write", path)?;
        Ok(Self {
            inner:           BufWriter::new(file),
            index:           BaiIndex::default(),
            last_coordinate: None,
        })
    }

    /// Serializes the accumulated BAI index and flushes the output stream.
    pub(super) fn write_index(mut self) -> Result<(), BamError> {
        const BAI_MAGIC: &[u8; 4] = b"BAI\x01";

        self.inner.write_all(BAI_MAGIC)?;

        let n_ref = i32::try_from(self.index.references.len()).map_err(|_| BamEncodingError::SizeOverflow {
            field:  "Number of references",
            target: NumberSizeTarget::MaxExclusive(1 << 31),
        })?;

        self.inner.write_all(&n_ref.to_le_bytes())?;

        for (ref_id, reference) in self.index.references.iter().enumerate() {
            reference
                .write_to(&mut self.inner)
                .with_context(format!("Error writing BAI reference {ref_id}"))?;
        }

        self.inner.write_all(&self.index.n_no_coord.to_le_bytes())?;

        self.inner.flush()?;

        Ok(())
    }

    /// Checks if the incoming coordinate is in the correct sort order.
    pub(super) fn is_coordinate_ordered(&self, ref_id: i32, pos0: i32) -> bool {
        let coordinate = coordinate(ref_id, pos0);

        match self.last_coordinate {
            Some(last) => coordinate >= last,
            None => true,
        }
    }

    /// Updates the last seen record coordinate.
    pub(super) fn update_last_coord(&mut self, ref_id: i32, pos0: i32) {
        self.last_coordinate = Some(coordinate(ref_id, pos0));
    }

    /// Pre-allocates the index based on the number of references.
    ///
    /// This is called immediately after writing the header, before any records
    /// are written.
    pub(super) fn initialize_index(&mut self, n_ref: usize) {
        self.index.references = (0..n_ref).map(|_| ReferenceIndex::default()).collect();
        self.index.n_no_coord = 0;
    }

    /// Adds the current record meta-data to the index.
    pub(super) fn add_record_to_index(
        &mut self, ref_id: i32, bin: u16, ref_interval: Option<&Range<u32>>, chunk: (u64, u64),
    ) {
        let Ok(ref_id) = usize::try_from(ref_id) else {
            self.index.n_no_coord += 1;
            return;
        };

        let Some(reference) = self.index.references.get_mut(ref_id) else {
            return;
        };

        reference.bins.entry(bin).or_default().push(chunk);
        if let Some(ref_interval) = ref_interval {
            reference.update_linear_index(ref_interval, chunk);
        }
    }
}

/// In-memory BAI index accumulated from encoded BAM records.
#[derive(Default)]
struct BaiIndex {
    references: Vec<ReferenceIndex>,
    n_no_coord: u64,
}

/// Index entries for one reference sequence.
#[derive(Default)]
struct ReferenceIndex {
    bins:         HashMap<u16, Vec<(u64, u64)>>,
    linear_index: Vec<u64>,
}

impl ReferenceIndex {
    fn update_linear_index(&mut self, ref_interval: &Range<u32>, chunk: (u64, u64)) {
        let first_window = (ref_interval.start >> LINEAR_INDEX_SHIFT) as usize;
        let last_window = ((ref_interval.end - 1) >> LINEAR_INDEX_SHIFT) as usize;

        if self.linear_index.len() <= last_window {
            self.linear_index.resize(last_window + 1, 0);
        }
        for window in &mut self.linear_index[first_window..=last_window] {
            if *window == 0 || chunk.0 < *window {
                *window = chunk.0;
            }
        }
    }

    // TODO: samtools merges adjacent/overlapping chunks within a bin before
    // writing. We currently write one chunk per record, which is spec-valid
    // but larger than samtools output. Optimize once round-trip tests pass.
    fn write_to<W: Write>(&self, writer: &mut W) -> Result<(), BamError> {
        let n_bin = i32::try_from(self.bins.len()).map_err(|_| BamEncodingError::SizeOverflow {
            field:  "Number of bins",
            target: NumberSizeTarget::MaxExclusive(1usize << 31),
        })?;

        writer.write_all(&n_bin.to_le_bytes())?;

        let mut bins = self.bins.iter().collect::<Vec<_>>();
        bins.sort_unstable_by_key(|&(&bin, _)| bin);

        for (&bin, chunks) in bins {
            let bin = u32::from(bin);
            writer.write_all(&bin.to_le_bytes())?;

            let n_chunk = i32::try_from(chunks.len()).map_err(|_| BamEncodingError::SizeOverflow {
                field:  "Number of chunks",
                target: NumberSizeTarget::MaxExclusive(1usize << 31),
            })?;

            writer.write_all(&n_chunk.to_le_bytes())?;

            for &chunk in chunks {
                writer.write_all(&chunk.0.to_le_bytes())?;
                writer.write_all(&chunk.1.to_le_bytes())?;
            }
        }

        let lin_indx_len = i32::try_from(self.linear_index.len()).map_err(|_| BamEncodingError::SizeOverflow {
            field:  "Linear index length",
            target: NumberSizeTarget::MaxExclusive(1usize << 31),
        })?;
        writer.write_all(&lin_indx_len.to_le_bytes())?;

        for ioffset in &self.linear_index {
            writer.write_all(&ioffset.to_le_bytes())?;
        }

        Ok(())
    }
}

fn coordinate(ref_id: i32, pos0: i32) -> (u8, i32, i32) {
    if ref_id < 0 { (1, 0, 0) } else { (0, ref_id, pos0) }
}
