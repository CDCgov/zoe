use std::env;
use zoe::prelude::*;

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    FastaReader::from_filename(filename)
        .unwrap_or_die("Translate file error!")
        .flatten()
        .map(|r| r.translate())
        .for_each(|fasta_record| print!("{fasta_record}"));
}
