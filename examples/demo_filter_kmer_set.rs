use std::env;
use zoe::{
    data::err::{ErrorWithContext, ResultWithErrorContext},
    kmer::encoders::three_bit::ThreeBitKmerSet,
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
    let mut kmer_set = ThreeBitKmerSet::<2>::new(2).with_context("Failed to build k-mer set")?;

    for record in fasta_reader {
        let nt = record?.sequence;
        kmer_set.insert_from_sequence(&nt);
    }

    let bad_kmer = kmer_set.encoder().encode_kmer(b"NN");
    let mut filtered_kmer_set = ThreeBitKmerSet::<2>::new(2).with_context("Cannot build k-mer counter")?;
    for kmer in kmer_set.iter_encoded() {
        if kmer != bad_kmer {
            filtered_kmer_set.insert_kmer(&kmer);
        }
    }

    // Check to see if it worked
    for kmer in filtered_kmer_set {
        println!("{kmer}")
    }

    Ok(())
}
