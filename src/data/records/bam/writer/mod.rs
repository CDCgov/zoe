//! Serialize [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) records
//! into a [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) stream.
//!
//! The streaming writer [`BamWriter`] accepts SAM header lines and
//! [`GetSamFields`] records and emits BGZF-wrapped BAM.

use crate::data::{
    bam::{
        encoder::PreparedBamRecord,
        error::{BamError, BamRecordError},
        header::Header,
        writer::{bai::BaiWriter, bgzf::BgzfWriter},
    },
    err::ResultWithErrorContext,
    sam::GetSamFields,
};
use std::{
    io::Write,
    {fs::File, path::Path},
};

mod bai;
mod bgzf;

pub use bgzf::BlockCompressor;
pub use bgzf::NoCompression;

/// Streaming writer for BAM records backed by a BGZF output stream.
///
/// Compression is controlled by the type parameter `C`, which implements
/// [`BlockCompressor`]. By default this type is [`NoCompression`], which always
/// writes stored-DEFLATE BGZF blocks. Downstream crates can provide their own
/// compression backend by implementing [`BlockCompressor`] and passing it to
/// [`with_compressor`].
///
/// Choose a constructor based on how much control is needed:
///
/// - [`from_path`]: create or truncate a BAM file at a filesystem path and use
///   the default stored-DEFLATE strategy.
/// - [`from_writer`]: wrap an existing [`Write`] sink and use the default
///   stored-DEFLATE strategy.
/// - [`with_compressor`]: replace the default compression strategy with a
///   caller-provided one.
/// - [`with_bai`]: enable a companion BAI index at a caller-provided path.
///
/// After construction, add all SAM header lines with [`write_header_line`],
/// then serialize [`GetSamFields`] records with [`write_record`]. The BAM
/// header is written lazily just before the first record, or during [`finish`]
/// if no records were written.
///
/// Calling [`finish`] is recommended because it is the only way to observe
/// errors from final header writing, pending BGZF block flushing, EOF marker
/// writing, and optional BAI index writing. Dropping the writer also attempts
/// finalization, but any error is ignored.
///
/// ## Restrictions Compared to SAM Format
///
/// The BAM writer is more restrictive than the SAM reader in several ways:
///
/// - **Header references**: `@SQ` lines must contain unique, NUL-free `SN`
///   values and positive signed-32-bit `LN` values. Other header lines are
///   preserved as text but are not parsed for record encoding.
/// - **QNAME**: Must match SAM/BAM's `[!-?A-~]{1,254}` rule.
/// - **RNAME**: Any value other than `*` must be declared in an `@SQ` header
///   line via [`write_header_line`].
/// - **SEQ/QUAL/CIGAR consistency**: If SEQ is not missing, SEQ length must
///   match the query-consuming length of a non-empty CIGAR. QUAL must be
///   missing when SEQ is missing; otherwise QUAL may be missing or must match
///   SEQ length.
/// - **CIGAR operations**: CIGAR strings must parse into valid operations with
///   non-zero increments, and each operation's length cannot exceed 268,435,455
///   (BAM's 28-bit CIGAR length limit). Records with more than 65,535 encoded
///   CIGAR operations use BAM's long-CIGAR placeholder and a synthesized
///   `CG:B:I` auxiliary field.
/// - **Auxiliary fields**: Optional fields must parse as SAM `A`, `i`, `f`,
///   `Z`, `H`, or `B` values. Strings (`Z` type) and hex data (`H` type) cannot
///   contain embedded NUL bytes. Integer values must fit BAM's signed or
///   unsigned 32-bit auxiliary encodings.
/// - **Flags**: Unsupported, reserved, paired-template, and mate-dependent flag
///   bits are cleared. The preserved bits are `0x4`, `0x10`, `0x100`, `0x200`,
///   `0x400`, and `0x800`.
/// - **Mate information**: The `RNEXT`, `PNEXT`, and `TLEN` fields from the
///   input [`GetSamFields`] are *currently* ignored; BAM records are always
///   written with `RNEXT = -1`, `PNEXT = -1`, and `TLEN = 0`. In the future
///   *Zoe* may implement these fields in SAM/BAM.
/// - **Sequence encoding**: `U` is encoded as `T`, and unrecognized sequence
///   symbols are encoded as `N`. The `U` to `T` mapping is an intentional *Zoe*
///   choice that deviates from the SAM/BAM spec.
/// - **Reference IDs**: Reference IDs must fit in BAM's signed 32-bit `refID`
///   fields.
/// - **Coordinates and indexing**: This writer uses classic BAI-compatible BAM
///   binning. For every mapped or placed-unmapped record, the reference
///   interval `[beg, end)` must lie within `[0, 2^29)`. See [Section
///   5](https://samtools.github.io/hts-specs/SAMv1.pdf).
///
/// [`from_path`]: BamWriter::from_path
/// [`from_writer`]: BamWriter::from_writer
/// [`with_compressor`]: BamWriter::with_compressor
/// [`with_bai`]: BamWriter::with_bai
/// [`write_header_line`]: BamWriter::write_header_line
/// [`write_record`]: BamWriter::write_record
/// [`finish`]: BamWriter::finish
pub struct BamWriter<W: Write = File, C: BlockCompressor = NoCompression> {
    /// BGZF output stream receiving BAM bytes.
    bgzf:        Option<BgzfWriter<W, C>>,
    /// Accumulated header lines and parsed reference dictionary.
    header:      Header,
    /// Whether the BAM header has already been serialized to `bgzf`.
    header_done: bool,
    /// Optional companion BAI index writer.
    bai:         Option<BaiWriter>,
}

impl BamWriter<File, NoCompression> {
    /// Creates a new BAM writer at `bam_path` using [`NoCompression`].
    ///
    /// The file is created or truncated immediately, while the BAM header is
    /// buffered until [`write_record`] or [`finish`] is called.
    ///
    /// ## Errors
    ///
    /// Returns an error if `bam_path` cannot be created for writing.
    /// Compression uses the default stored-DEFLATE strategy.
    ///
    /// [`write_record`]: BamWriter::write_record
    /// [`finish`]: BamWriter::finish
    pub fn from_path(bam_path: impl AsRef<Path>) -> Result<Self, BamError> {
        let file = File::create(&bam_path).with_path_context("Cannot create BAM file to write", bam_path)?;
        Ok(Self::from_writer(file))
    }
}

impl<W: Write> BamWriter<W, NoCompression> {
    /// Wraps `writer` in a BAM writer using [`NoCompression`].
    ///
    /// The BAM header is buffered until [`write_record`] or [`finish`] is
    /// called. Call [`with_compressor`] to replace the default stored-DEFLATE
    /// strategy before writing data.
    ///
    /// [`write_record`]: BamWriter::write_record
    /// [`finish`]: BamWriter::finish
    /// [`with_compressor`]: BamWriter::with_compressor
    pub fn from_writer(writer: W) -> Self {
        BamWriter {
            bgzf:        Some(BgzfWriter::new(writer)),
            header:      Header::default(),
            header_done: false,
            bai:         None,
        }
    }

    /// Replaces the default [`NoCompression`] strategy with `compressor`.
    ///
    /// ## Errors
    ///
    /// Returns an error if attempting to set the compression strategy after
    /// serialization has already begun.
    pub fn with_compressor<C: BlockCompressor>(mut self, compressor: C) -> Result<BamWriter<W, C>, BamError> {
        if self.header_done {
            return Err(BamError::CannotChangeConfiguration);
        }

        // `self` implements `Drop`, so fields must be extracted via
        // `take`/`mem::take` rather than moved out directly
        Ok(BamWriter {
            bgzf:        self.bgzf.take().map(|b| b.with_compressor(compressor)),
            header:      std::mem::take(&mut self.header),
            header_done: self.header_done,
            bai:         self.bai.take(),
        })
    }
}

impl<W: Write, C: BlockCompressor> BamWriter<W, C> {
    /// Enables BAI indexing and stores the result in `path` when finalized.
    ///
    /// ## Errors
    ///
    /// Returns an error if `path` cannot be created or if attempting to enable
    /// indexing after serialization has already begun.
    pub fn with_bai(mut self, path: impl AsRef<Path>) -> Result<Self, BamError> {
        if self.header_done {
            return Err(BamError::CannotChangeConfiguration);
        }

        let bai = BaiWriter::from_path(path)?;
        self.bai = Some(bai);
        Ok(self)
    }

    /// Adds one SAM header line to the pending BAM header.
    ///
    /// This method must be called before the BAM header is serialized. Header
    /// lines are preserved as text in the BAM header block; `@SQ` lines are
    /// also parsed for the reference dictionary used by [`write_record`]. The
    /// header is serialized after the first record is successfully encoded.
    ///
    /// ## Errors
    ///
    /// Returns an error if the BAM header has already been written, or if a
    /// `@SQ` line is malformed, duplicated, or exceeds BAM's limits. Specific
    /// limitations for `@SQ` lines include:
    ///
    /// - The `SN` (sequence name) field cannot contain an embedded NUL byte
    /// - The `SN` field must be unique across all `@SQ` lines
    /// - The `LN` (sequence length) field must be a positive signed 32-bit
    ///   integer
    /// - The total number of reference sequences cannot exceed 2,147,483,647
    ///
    /// [`write_record`]: BamWriter::write_record
    pub fn write_header_line(&mut self, header_line: &str) -> Result<(), BamError> {
        if self.header_done {
            Err(BamError::HeaderAlreadyWritten)
        } else {
            self.header.push_header_line(header_line)?;
            Ok(())
        }
    }

    /// Serializes one [`GetSamFields`] record as BAM.
    ///
    /// On the first call, this validates and encodes `record` before writing
    /// the accumulated BAM header. Any `RNAME` value other than `*` must have
    /// been introduced by an earlier `@SQ` line passed to
    /// [`write_header_line`].
    ///
    /// ## Errors
    ///
    /// Returns an error if the BAM header cannot be written, if `record` cannot
    /// be represented in BAM, or if writing to the output fails. Common reasons
    /// a record cannot be represented include:
    ///
    /// - `QNAME` does not match SAM/BAM's `[!-?A-~]{1,254}` rule
    /// - `RNAME` is an empty string
    /// - Non-empty `SEQ` length does not match the query-consuming length of a
    ///   non-empty CIGAR
    /// - `QUAL` is not missing when `SEQ` is missing
    /// - `QUAL` is neither missing nor the same length as `SEQ`
    /// - CIGAR is malformed, contains zero-length operations, or contains
    ///   operations with lengths exceeding BAM's 28-bit bounds
    /// - Auxiliary fields are malformed, contain out-of-range integer values,
    ///   contain embedded NUL bytes in Z or H values, or contain a reserved
    ///   `CG` tag.
    /// - Any field encoding would exceed BAM's representable limits. See
    ///   [Section 4.2](https://samtools.github.io/hts-specs/SAMv1.pdf) of the
    ///   BAM file format specs.
    ///
    /// [`write_header_line`]: BamWriter::write_header_line
    pub fn write_record(&mut self, record: impl GetSamFields) -> Result<(), BamError> {
        if self.bgzf.is_none() {
            return Err(BamError::WriterFinalized);
        }

        let prepared_record =
            PreparedBamRecord::new(&self.header, &record).map_err(|source| BamError::record(record.qname(), source))?;

        if let Some(bai) = self.bai.as_ref()
            && !bai.is_coordinate_ordered(prepared_record.ref_id, prepared_record.pos0)
        {
            return Err(BamError::record(record.qname(), BamRecordError::UnsortedRecord));
        }

        let encoded_record = prepared_record.encode(record.qname())?;

        let mut bgzf = self.bgzf.take().ok_or(BamError::WriterFinalized)?;
        self.write_header_if_pending(&mut bgzf)?;

        let chunk_beg = bgzf.virtual_offset()?;
        bgzf.write_all(&encoded_record).map_err(BamError::from)?;
        let chunk_end = bgzf.virtual_offset()?;

        if let Some(bai) = self.bai.as_mut() {
            bai.add_record_to_index(
                prepared_record.ref_id,
                prepared_record.bin,
                prepared_record.ref_interval.as_ref(),
                (chunk_beg, chunk_end),
            );
            bai.update_last_coord(prepared_record.ref_id, prepared_record.pos0);
        }

        self.bgzf = Some(bgzf);

        Ok(())
    }

    /// Attempts to finalize the BAM stream and append the mandatory BGZF EOF
    /// marker.
    ///
    /// If no records were written, this still serializes the accumulated BAM
    /// header before finalizing the file. Once finalization starts, [`Drop`]
    /// will not retry it, even if this method returns an error.
    ///
    /// ## Errors
    ///
    /// Returns an error if the pending BAM header or BGZF stream cannot be
    /// finalized, or if the optional BAI index cannot be written.
    pub fn finish(mut self) -> Result<(), BamError> {
        let Some(mut bgzf) = self.bgzf.take() else {
            return Ok(());
        };

        self.write_header_if_pending(&mut bgzf)?;
        bgzf.finish().map_err(BamError::from)?;

        if let Some(bai) = self.bai.take() {
            bai.write_index().with_context("Error writing BAI file")?;
        }

        Ok(())
    }

    /// Writes the accumulated header if pending. `header_done` is set to true
    /// on an attempted write.
    fn write_header_if_pending(&mut self, bgzf: &mut BgzfWriter<W, C>) -> Result<(), BamError> {
        if !self.header_done {
            self.header_done = true;
            self.header.write_to(bgzf)?;

            if let Some(bai) = self.bai.as_mut() {
                bai.initialize_index(self.header.ref_count());
            }
        }

        Ok(())
    }
}

impl<W: Write, C: BlockCompressor> Drop for BamWriter<W, C> {
    /// Attempts to write any pending header, BGZF EOF marker, and BAI index,
    /// ignoring all errors.
    fn drop(&mut self) {
        let Some(mut bgzf) = self.bgzf.take() else {
            return;
        };

        if self.write_header_if_pending(&mut bgzf).is_err() {
            return;
        }

        if bgzf.finish().is_err() {
            return;
        }

        if let Some(bai) = self.bai.take() {
            let _ = bai.write_index();
        }
    }
}
