use zoe::data::{
    bam::{
        error::BamError,
        writer::{BamWriter, BlockCompressor},
    },
    sam::{SAMReader, SamRow},
};

fn main() -> Result<(), BamError> {
    let sam_path = "examples/example.sam";
    let bam_file = "examples/example.bam";

    let sam_reader = SAMReader::from_path(sam_path)?;
    let mut bam_writer = BamWriter::from_path_with_compressor(bam_file, CustomCompressor)?;

    for line in sam_reader {
        let line = line.map_err(BamError::from)?;

        match line {
            SamRow::Header(header_line) => bam_writer.write_header_line(&header_line)?,
            SamRow::Data(record) => {
                bam_writer.write_record(&record)?;
            }
        }
    }

    bam_writer.finish()
}

struct CustomCompressor;

impl BlockCompressor for CustomCompressor {
    /// Define your custom function to compress `input` into `output` as a
    /// complete raw-DEFLATE stream.
    ///
    /// On success, return `Some(len)` where `len` is the number of bytes in
    /// `output` that make up the encoded stream. Returning `None` requests the
    /// stored-DEFLATE fallback path.
    ///
    /// This example mimics the existing [`NoCompression`] option.
    fn compress(&mut self, _input: &[u8], _output: &mut Vec<u8>) -> std::io::Result<Option<usize>> {
        // Placeholder currently defaults to stored-DEFLATE
        Ok(None)
    }
}
