use std::env;
use zoe::data::fasta::*;

fn main() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let fasta_reader = FastaReader::from_filename(&filename)?;

    fasta_reader
        .map(|r| r.translate())
        .for_each(|fasta_record| print!("{fasta_record}"));

    Ok(())
}
