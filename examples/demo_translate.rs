use std::env;
use zoe::{
    data::err::{ErrorWithContext, ResultWithErrorContext},
    prelude::*,
};

fn main() -> Result<(), ErrorWithContext> {
    let args: Vec<String> = env::args().collect();

    let path = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let reader = FastaReader::from_path(&path).with_path_context("Failed to open FASTA file", path)?;

    for record in reader {
        let record = record?;
        let translated = record.translate();
        print!("{translated}")
    }

    Ok(())
}
