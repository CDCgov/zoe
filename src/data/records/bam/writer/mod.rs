//! Serialize [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) records
//! into a [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) stream.
//!
//! The streaming writer [`BamWriter`] accepts SAM header lines and [`SamData`]
//! records and emits BGZF-wrapped BAM.

use crate::data::{
    bam::{encoder::PreparedBamRecord, error::BamError, header::Header, writer::bgzf::BgzfWriter},
    err::ResultWithErrorContext,
    sam::SamData,
};
use std::{
    io::Write,
    {fs::File, path::Path},
};

mod bgzf;

pub use bgzf::BlockCompressor;
pub use bgzf::NoCompression;

/// Streaming writer for BAM records backed by a BGZF output stream.
///
/// Compression is controlled by the type parameter `C`, which implements
/// [`BlockCompressor`]. By default this type is [`NoCompression`], which always
/// writes stored-DEFLATE BGZF blocks. Downstream crates can provide their own
/// compression backend by implementing [`BlockCompressor`] and constructing a
/// writer with [`from_path_with_compressor`] or
/// [`from_writer_with_compressor`].
///
/// Choose a constructor based on how much control is needed:
///
/// - [`from_path`]: create or truncate a BAM file at a filesystem path and use
///   the default stored-DEFLATE strategy.
/// - [`from_path_with_compressor`]: create or truncate a BAM file at a path and
///   use a caller-provided compression strategy.
/// - [`from_writer_with_compressor`]: wrap an existing [`Write`] sink and use a
///   caller-provided compression strategy.
///
/// After construction, add all SAM header lines with [`write_header_line`],
/// then serialize [`SamData`] records with [`write_record`]. The BAM header is
/// written lazily just before the first record, or during [`finish`] if no
/// records were written.
///
/// Calling [`finish`] is recommended because it is the only way to observe
/// errors from final header writing, pending BGZF block flushing, and EOF
/// marker writing. Dropping the writer also attempts finalization, but any
/// error is ignored.
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
///   input [`SamData`] are *currently* ignored; BAM records are always written
///   with `RNEXT = -1`, `PNEXT = -1`, and `TLEN = 0`. In the future *Zoe* may
///   implement these fields in SAM/BAM.
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
/// [`from_path_with_compressor`]: BamWriter::from_path_with_compressor
/// [`from_writer_with_compressor`]: BamWriter::from_writer_with_compressor
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
        Self::from_path_with_compressor(bam_path, NoCompression)
    }
}

impl<C: BlockCompressor> BamWriter<File, C> {
    /// Creates a new BAM writer at `bam_path` using `compressor`.
    ///
    /// This constructor is the path-based variant for callers who want to
    /// control compression while still letting [`BamWriter`] open the output
    /// file.
    ///
    /// The file is created or truncated immediately, while the BAM header is
    /// buffered until [`write_record`] or [`finish`] is called.
    ///
    /// ## Errors
    ///
    /// Returns an error if `bam_path` cannot be created for writing.
    ///
    /// [`write_record`]: BamWriter::write_record
    /// [`finish`]: BamWriter::finish
    pub fn from_path_with_compressor(bam_path: impl AsRef<Path>, compressor: C) -> Result<Self, BamError> {
        let file = File::create(&bam_path).with_path_context("Cannot create BAM file to write", bam_path)?;
        Ok(Self::from_writer_with_compressor(file, compressor))
    }
}

impl<W: Write, C: BlockCompressor> BamWriter<W, C> {
    /// Creates a new BAM writer over an existing writer and compressor.
    ///
    /// The wrapped writer receives BGZF blocks, and compression decisions are
    /// delegated to `compressor`. It is not recommended to pass a buffered
    /// writer, since [`BamWriter`] performs its own buffering.
    pub fn from_writer_with_compressor(writer: W, compressor: C) -> Self {
        BamWriter {
            bgzf:        Some(BgzfWriter::with_compressor(writer, compressor)),
            header:      Header::default(),
            header_done: false,
        }
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

    /// Serializes one [`SamData`] record as BAM.
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
    pub fn write_record(&mut self, record: &SamData) -> Result<(), BamError> {
        if self.bgzf.is_none() {
            return Err(BamError::WriterFinalized);
        }

        let prepared_record = PreparedBamRecord::new(&self.header, record)
            .map_err(|source| BamError::record(record.qname.as_str(), source))?;
        let encoded_record = prepared_record.encode(&record.qname)?;

        let mut bgzf = self.bgzf.take().ok_or(BamError::WriterFinalized)?;
        self.write_header_if_pending(&mut bgzf)?;
        bgzf.write_all(&encoded_record).map_err(BamError::from)?;
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
    /// Returns an error if the pending BAM header cannot be written or if the
    /// BGZF stream cannot be finalized.
    pub fn finish(mut self) -> Result<(), BamError> {
        let Some(mut bgzf) = self.bgzf.take() else {
            return Ok(());
        };

        self.write_header_if_pending(&mut bgzf)?;
        bgzf.finish().map_err(BamError::from)
    }

    /// Writes the accumulated header if pending. `header_done` is set to true
    /// on an attempted write.
    fn write_header_if_pending(&mut self, bgzf: &mut BgzfWriter<W, C>) -> Result<(), BamError> {
        if !self.header_done {
            self.header_done = true;
            self.header.write_to(bgzf)?;
        }

        Ok(())
    }
}

impl<W: Write, C: BlockCompressor> Drop for BamWriter<W, C> {
    /// Attempts to write any pending header and BGZF EOF marker, ignoring all
    /// errors.
    fn drop(&mut self) {
        let Some(mut bgzf) = self.bgzf.take() else {
            return;
        };

        if self.write_header_if_pending(&mut bgzf).is_ok() {
            let _ = bgzf.finish();
        }
    }
}
