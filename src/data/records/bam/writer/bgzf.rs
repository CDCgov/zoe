use crate::data::err::ResultWithErrorContext;
use std::io::Write;

/// Maximum size of one complete BGZF block, in bytes.
const MAX_BGZF_BLOCK_SIZE: usize = 65_536;
/// Bytes in the fixed BGZF/gzip header, including the `BC` subfield and `BSIZE`
/// value.
const BGZF_HEADER_SIZE: usize = 18;
/// Bytes in the gzip footer: the CRC-32 and uncompressed input size.
const GZIP_FOOTER_SIZE: usize = 8;
/// Maximum raw-DEFLATE (compressed) stream length that fits between a BGZF
/// header and gzip footer.
const MAX_DEFLATE_SIZE: usize = MAX_BGZF_BLOCK_SIZE - BGZF_HEADER_SIZE - GZIP_FOOTER_SIZE;
/// Bytes added by a final
/// [stored-DEFLATE](https://en.wikipedia.org/wiki/Deflate#:~:text=(sometimes%20called-,stored,-).%20Any%20bits%20up)
/// (uncompressed data) block: its header, `LEN`, and `NLEN`.
const STORED_DEFLATE_OVERHEAD: usize = 5;
/// Maximum uncompressed payload buffered for one BGZF block.
///
/// A
/// [stored-DEFLATE](https://en.wikipedia.org/wiki/Deflate#:~:text=(sometimes%20called-,stored,-).%20Any%20bits%20up)
/// (uncompressed data) block of this size exactly fits the BGZF block limit.
/// Compressed blocks use this same input size and fall back to stored DEFLATE
/// when their compressed raw-DEFLATE stream does not fit.
const MAX_BGZF_PAYLOAD: usize = MAX_DEFLATE_SIZE - STORED_DEFLATE_OVERHEAD;

/// Compression hook for producing a raw-DEFLATE (compressed) stream from one
/// BGZF payload.
///
/// Implementors may return `Ok(None)` to indicate that the payload should be
/// written as
/// [stored-DEFLATE](https://en.wikipedia.org/wiki/Deflate#:~:text=(sometimes%20called-,stored,-).%20Any%20bits%20up)
/// (uncompressed data) instead.
pub trait BlockCompressor {
    /// Attempts to compress `input` into `output` as a complete raw-DEFLATE
    /// (compressed) stream.
    ///
    /// On success, returns `Some(len)` where `len` is the number of bytes in
    /// `output` that make up the encoded stream. Returning `None` requests the
    /// [stored-DEFLATE](https://en.wikipedia.org/wiki/Deflate#:~:text=(sometimes%20called-,stored,-).%20Any%20bits%20up)
    /// (uncompressed data) fallback path.
    ///
    /// # Errors
    ///
    /// Returns an error if compression fails for the current payload.
    fn compress(&mut self, input: &[u8], output: &mut Vec<u8>) -> std::io::Result<Option<usize>>;
}

/// Default compressor that never emits compressed bytes.
pub struct NoCompression;

impl BlockCompressor for NoCompression {
    #[inline]
    fn compress(&mut self, _input: &[u8], _output: &mut Vec<u8>) -> std::io::Result<Option<usize>> {
        Ok(None)
    }
}

/// Buffered writer that repackages a byte stream into BGZF blocks.
///
/// BAM bytes are written through this adapter so record serialization can
/// stream directly to the output while the writer handles DEFLATE block
/// boundaries and the mandatory BGZF EOF marker.
///
/// Compression is delegated to a pluggable [`BlockCompressor`]. The default
/// type parameter uses [`NoCompression`], which always emits
/// [stored-DEFLATE](https://en.wikipedia.org/wiki/Deflate#:~:text=(sometimes%20called-,stored,-).%20Any%20bits%20up)
/// (uncompressed data) blocks.
pub(super) struct BgzfWriter<W: Write, C: BlockCompressor = NoCompression> {
    /// Wrapped writer that receives complete BGZF blocks.
    inner:            W,
    /// Pending uncompressed payload bytes for the next BGZF block.
    payload:          Vec<u8>,
    /// Reusable scratch buffer for one encoded BGZF block.
    block:            Vec<u8>,
    /// Pluggable raw-DEFLATE (compressed) encoder.
    compressor:       C,
    /// Reusable raw-DEFLATE (compressed) output storage for compressed blocks.
    compressed_block: Vec<u8>,
}

impl<W: Write, C: BlockCompressor> BgzfWriter<W, C> {
    /// Creates an empty BGZF writer over `inner`.
    pub(super) fn with_compressor(inner: W, compressor: C) -> Self {
        Self {
            inner,
            payload: Vec::with_capacity(MAX_BGZF_PAYLOAD),
            block: Vec::with_capacity(MAX_BGZF_BLOCK_SIZE),
            compressor,
            compressed_block: Vec::with_capacity(MAX_DEFLATE_SIZE),
        }
    }

    /// Attempts to flush any buffered payload, write the BGZF EOF marker, and
    /// flush the wrapped writer.
    ///
    /// ## Errors
    ///
    /// Returns an error if writing a pending BGZF block, writing the EOF
    /// marker, or flushing the wrapped writer fails.
    pub(super) fn finish(mut self) -> std::io::Result<()> {
        /// BGZF EOF marker block.
        const BGZF_EOF_MARKER: [u8; 28] = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00,
            0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        self.flush_block()?;
        self.inner
            .write_all(&BGZF_EOF_MARKER)
            .with_context("Error writing BGZF EOF marker.")?;
        self.inner.flush()?;

        Ok(())
    }

    /// Encodes the buffered payload as one BGZF block and clears the buffer.
    fn flush_block(&mut self) -> std::io::Result<()> {
        if self.payload.is_empty() {
            return Ok(());
        }
        self.write_bgzf_block()?;
        self.payload.clear();
        Ok(())
    }

    /// Writes one payload as a BGZF block, using compressed raw DEFLATE when
    /// worthwhile, and
    /// [stored-DEFLATE](https://en.wikipedia.org/wiki/Deflate#:~:text=(sometimes%20called-,stored,-).%20Any%20bits%20up)
    /// (uncompressed data) otherwise.
    fn write_bgzf_block(&mut self) -> std::io::Result<()> {
        self.compressed_block.clear();
        let deflate_len = self
            .compressor
            .compress(&self.payload, &mut self.compressed_block)
            .with_context("Failed to compress a BGZF block.")?;

        // DEFLATE is allowed to expand incompressible payloads. Prefer the
        // existing stored representation unless compression makes the block
        // smaller. If compression is unavailable or unsuitable, fall back to
        // [stored-DEFLATE](https://en.wikipedia.org/wiki/Deflate#:~:text=(sometimes%20called-,stored,-).%20Any%20bits%20up)
        // (uncompressed data).
        if let Some(deflate_len) = deflate_len
            && deflate_len < self.payload.len() + STORED_DEFLATE_OVERHEAD
        {
            return self.write_compressed_bgzf_block(deflate_len);
        }

        self.write_stored_bgzf_block()
    }

    /// Writes one `payload` as a single
    /// [stored-DEFLATE](https://en.wikipedia.org/wiki/Deflate#:~:text=(sometimes%20called-,stored,-).%20Any%20bits%20up)
    /// (uncompressed data) BGZF block.
    ///
    /// The `payload` is expected to be no larger than [`MAX_BGZF_PAYLOAD`].
    fn write_stored_bgzf_block(&mut self) -> std::io::Result<()> {
        /// DEFLATE header byte for a final (the only) stored block: bit `0` is
        /// `BFINAL = 1`,bits `1..=2` are `BTYPE = 00` and bits `3..=7` are zero
        /// padding. Least-significant bit ordering.
        const STORED_DEFLATE_HEADER: u8 = 0b0000_0001;

        let payload = &self.payload;
        let len = u16::try_from(payload.len())
            .map_err(std::io::Error::other)
            .with_context("BGZF payload length cannot be represented with u16.")?;
        let nlen = !len;
        // (BGZF header len `12` + gzip header len `6` + deflate header len `5`
        // + payload len + gzip footer len `8`) - 1
        let bsize = len
            .checked_add(30)
            .ok_or_else(|| std::io::Error::other("BGZF block size exceeds 64 KiB"))?;
        let crc = crc32_table(payload);

        self.block.clear();

        // BGZF + gzip header
        Self::append_bgzf_header(&mut self.block, bsize);

        // stored DEFLATE block header
        self.block.extend_from_slice(&[STORED_DEFLATE_HEADER]);
        self.block.extend_from_slice(&len.to_le_bytes());
        self.block.extend_from_slice(&nlen.to_le_bytes());

        // payload
        self.block.extend_from_slice(payload);

        // gzip footer
        Self::append_bgzf_footer(&mut self.block, crc, len);

        self.inner.write_all(&self.block)
    }

    /// Writes one payload as a BGZF block containing an already-finished raw
    /// DEFLATE stream.
    fn write_compressed_bgzf_block(&mut self, deflate_len: usize) -> std::io::Result<()> {
        let deflate = self
            .compressed_block
            .get(..deflate_len)
            .ok_or_else(|| std::io::Error::other("Compressed DEFLATE length exceeds output buffer."))?;

        let payload = &self.payload;
        let len = u16::try_from(payload.len())
            .map_err(std::io::Error::other)
            .with_context("BGZF payload length cannot be represented with u16.")?;
        let deflate_len = u16::try_from(deflate.len())
            .map_err(std::io::Error::other)
            .with_context("Raw DEFLATE length cannot be represented with u16.")?;
        // (BGZF header len `12` + gzip header len `6` + raw DEFLATE len
        // + gzip footer len `8`) - 1
        let bsize = deflate_len
            .checked_add(25)
            .ok_or_else(|| std::io::Error::other("BGZF block size exceeds 64 KiB"))?;
        let crc = crc32_table(payload);

        self.block.clear();

        // BGZF + gzip header
        Self::append_bgzf_header(&mut self.block, bsize);

        // raw DEFLATE data
        self.block.extend_from_slice(deflate);

        // gzip footer
        Self::append_bgzf_footer(&mut self.block, crc, len);

        self.inner.write_all(&self.block)
    }

    /// Appends the fixed BGZF/gzip prefix, including `BSIZE`.
    fn append_bgzf_header(block: &mut Vec<u8>, bsize: u16) {
        /// gzip `IDentifier1`.
        const GZIP_ID1: u8 = 31;
        /// gzip `IDentifier2`.
        const GZIP_ID2: u8 = 139;
        /// gzip `Compression Method` for DEFLATE.
        const GZIP_CM_DEFLATE: u8 = 8;
        /// gzip `FLaGs`: the BGZF extra field is present.
        const GZIP_FLG_FEXTRA: u8 = 4;
        /// gzip `eXtra FLags`.
        const GZIP_XFL: u8 = 0;
        /// gzip `eXtra LENgth`: the BGZF extra field is six bytes.
        const BGZF_XLEN: [u8; 2] = 6_u16.to_le_bytes();
        /// BGZF extra subfield identifier.
        const BGZF_SI: [u8; 2] = *b"BC";
        /// BGZF extra subfield length: the `BSIZE` value is a `u16`.
        const BGZF_SLEN: [u8; 2] = 2_u16.to_le_bytes();

        // gzip `Modification TIME`; `mtime = 0` means no time stamp is
        // available.
        let mtime = 0_u32.to_le_bytes();
        // gzip `Operating System`; `os = 255` means unknown.
        let os = 255_u8;

        // BGZF Header
        block.extend_from_slice(&[GZIP_ID1, GZIP_ID2, GZIP_CM_DEFLATE, GZIP_FLG_FEXTRA]);
        block.extend_from_slice(&mtime);
        block.extend_from_slice(&[GZIP_XFL, os]);
        block.extend_from_slice(&BGZF_XLEN);

        // gzip header
        block.extend_from_slice(&BGZF_SI);
        block.extend_from_slice(&BGZF_SLEN);
        block.extend_from_slice(&bsize.to_le_bytes());
    }

    /// Appends the gzip footer with CRC-32 and uncompressed payload length.
    fn append_bgzf_footer(block: &mut Vec<u8>, crc: u32, uncompressed_len: u16) {
        block.extend_from_slice(&crc.to_le_bytes());
        block.extend_from_slice(&u32::from(uncompressed_len).to_le_bytes());
    }
}

impl<W: Write, C: BlockCompressor> Write for BgzfWriter<W, C> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.write_all(buf)?;
        Ok(buf.len())
    }

    /// Buffers all bytes from `buf`, flushing full BGZF blocks as needed.
    fn write_all(&mut self, mut buf: &[u8]) -> std::io::Result<()> {
        while !buf.is_empty() {
            let take = (MAX_BGZF_PAYLOAD - self.payload.len()).min(buf.len());
            self.payload.extend_from_slice(&buf[..take]);
            buf = &buf[take..];
            if self.payload.len() == MAX_BGZF_PAYLOAD {
                self.flush_block()?;
            }
        }
        Ok(())
    }

    /// Flushes any pending BGZF block and then flushes the wrapped writer.
    fn flush(&mut self) -> std::io::Result<()> {
        self.flush_block()?;
        self.inner.flush()
    }
}

/// Returns the `CRC-32` checksum written in each BGZF block footer based on a
/// lookup table for each byte.
fn crc32_table(bytes: &[u8]) -> u32 {
    /// Initial CRC-32 polynomial.
    const CRC_INIT: u32 = 0xFFFF_FFFF;
    /// CRC-32 lookup table for each byte.
    const CRC32_TABLE: [u32; 256] = make_crc_table();

    let mut crc = CRC_INIT;

    for &byte in bytes {
        let index = ((crc ^ u32::from(byte)) & 0xff) as usize;
        crc = (crc >> 8) ^ CRC32_TABLE[index];
    }

    !crc
}

/// Computes the lookup table for the `CRC-32` checksum for each byte.
const fn make_crc_table() -> [u32; 256] {
    /// IEEE CRC-32 generator polynomial.
    const POLY: u32 = 0xEDB8_8320;

    let mut table = [0u32; 256];

    let mut i = 0;
    while i < 256 {
        let mut crc = i;

        let mut bit = 0;
        while bit < 8 {
            if (crc & 1) != 0 {
                crc = (crc >> 1) ^ POLY;
            } else {
                crc >>= 1;
            }
            bit += 1;
        }

        table[i as usize] = crc;
        i += 1;
    }

    table
}
