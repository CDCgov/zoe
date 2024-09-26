use std::env;
use zoe::prelude::*;

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let fasta_reader = FastaReader::from_filename(filename).unwrap_or_die("RevComp file error!");

    for record in fasta_reader {
        let mut record = record.unwrap_or_die("RevComp FastaReader error.");
        record.reverse_complement();
        print!("{record}");
    }
}
