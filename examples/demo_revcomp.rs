use std::env;
use zoe::{data::err::ErrorWithContext, prelude::*};

fn main() -> Result<(), ErrorWithContext> {
    let args: Vec<String> = env::args().collect();

    let path = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let fasta_reader = FastaReader::from_path(path)?;

    for record in fasta_reader {
        let mut record = record?;
        record.reverse_complement();
        print!("{record}");
    }

    Ok(())
}
