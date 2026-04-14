use std::env;
use zoe::prelude::*;

fn main() {
    let args: Vec<String> = env::args().collect();

    let path = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    FastaReader::from_path(path)
        .unwrap_or_fail()
        .flatten()
        .map(|r| r.translate())
        .for_each(|fasta_record| print!("{fasta_record}"));
}
