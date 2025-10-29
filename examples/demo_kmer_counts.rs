use std::env;
use zoe::{kmer::encoders::three_bit::ThreeBitKmerCounter, prelude::*};

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let fasta_reader = FastaReader::from_filename(filename).unwrap_or_die("k-mer counts file error!");
    let mut kmer_counts = ThreeBitKmerCounter::<2>::new(2).unwrap_or_die("Cannot build k-mer counter.");

    for record in fasta_reader {
        let nt = record.unwrap_or_die("FastaReader error.").sequence;
        kmer_counts.tally_from_sequence(&nt);
    }

    for (kmer, count) in kmer_counts {
        println!("{kmer}\t{count}");
    }
}
