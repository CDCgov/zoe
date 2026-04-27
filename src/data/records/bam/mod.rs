//! Write [BAM](https://en.wikipedia.org/wiki/BAM_(file_format)) files from
//! [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) data.
//!
//! This module provides [`BamWriter`] for streaming [`SamData`] records to a
//! BAM output stream.
//!
//! Output is BGZF-wrapped BAM with the EOF marker. Compression is controlled by
//! the writer: the default [`BamWriter`] uses stored-DEFLATE, and downstream
//! crates can provide custom raw-DEFLATE backends through [`BlockCompressor`].
//!
//! Alignment records are validated more strictly than the SAM reader alone can
//! validate, because several SAM values have narrower or more structured
//! encodings in BAM.
//!
//! [`SamData`]: crate::data::sam::SamData
//! [`BamWriter`]: crate::data::bam::writer::BamWriter
//! [`BlockCompressor`]: crate::data::records::bam::writer::BlockCompressor

mod encoder;
pub mod error;
mod header;
pub mod writer;
