use std::env;
use zoe::prelude::*;

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fastq".to_owned()
    };

    let fastq_reader = FastQReader::from_filename(filename)
        .unwrap_or_die("FastQ file error!")
        .flatten();

    for mut record in fastq_reader {
        record.recode_iupac_to_actgn();
        print!("{record}");
    }
}
