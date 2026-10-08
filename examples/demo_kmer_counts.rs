use std::env;
use zoe::{
    data::err::{ErrorWithContext, ResultWithErrorContext},
    kmer::encoders::three_bit::ThreeBitKmerCounter,
    prelude::*,
};

fn main() -> Result<(), ErrorWithContext> {
    let args: Vec<String> = env::args().collect();

    let path = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let fasta_reader = FastaReader::from_path(path)?;
    let mut kmer_counts = ThreeBitKmerCounter::<2>::new(2).with_context("Cannot build k-mer counter")?;

    for record in fasta_reader {
        let nt = record?.sequence;
        kmer_counts.tally_from_sequence(&nt);
    }

    for (kmer, count) in kmer_counts {
        println!("{kmer}\t{count}");
    }

    Ok(())
}
