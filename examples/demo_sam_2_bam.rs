use zoe::data::{
    bam::{
        error::BamError,
        writer::{BamWriter, BlockCompressor},
    },
    sam::{SAMReader, SamDataSort, SamRow},
};

fn main() -> Result<(), BamError> {
    let sam_path = "examples/example.sam";
    let bam_path = "examples/example.bam";
    let bai_path = "examples/example.bai";

    let sam_reader = SAMReader::from_path(sam_path)?;

    // Collect and coordinate sort the SAM records.
    let mut headers = Vec::new();
    let mut records = Vec::new();
    for line in sam_reader {
        match line.map_err(BamError::from)? {
            SamRow::Header(header) => headers.push(header),
            SamRow::Data(record) => records.push(record),
        }
    }
    let coordinate_header = records.coordinate_sort(&headers);

    // Create the `BamWriter` with a custom compression backend and with bai
    // enabled.
    let mut bam_writer = BamWriter::from_path(bam_path)?
        .with_compressor(CustomCompressor)?
        .with_bai(bai_path)?;

    // The first header line is `"@HD\tVN:1.6\tSO:coordinate"` which indicates
    // the records have been coordinate sorted.
    bam_writer.write_header_line(coordinate_header)?;
    for header in headers.iter().filter(|header| !header.starts_with("@HD\t")) {
        bam_writer.write_header_line(header)?;
    }

    for record in records {
        bam_writer.write_record(&record)?;
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
