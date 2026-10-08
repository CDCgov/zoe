use std::env;
use zoe::{data::err::ErrorWithContext, prelude::*};

fn main() -> Result<(), ErrorWithContext> {
    let args: Vec<String> = env::args().collect();

    let path = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fastq".to_owned()
    };

    let fastq_reader = FastQReader::from_path(path)?;

    for record in fastq_reader {
        let mut record = record?;
        record.recode_dna_reads();
        print!("{record}");
    }

    Ok(())
}
