use std::env;
use zoe::prelude::*;

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let fasta_reader = FastaReader::from_filename(&filename)
        .unwrap_or_die("RevComp file error!")
        .filter_map(|f| f.ok());

    for mut record in fasta_reader {
        record.reverse_complement();
        print!("{record}");
    }
}
