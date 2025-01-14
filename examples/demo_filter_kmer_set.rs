use std::env;

use zoe::{kmer::ThreeBitKmerSet, prelude::*};

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let fasta_reader = FastaReader::from_filename(filename).unwrap_or_die("k-mer counts file error!");
    let mut kmer_set = ThreeBitKmerSet::<2>::new(2).unwrap_or_die("Cannot build k-mer counter.");

    for record in fasta_reader {
        let nt = record.unwrap_or_die("FastaReader error.").sequence;
        kmer_set.insert_from_sequence(&nt);
    }

    let bad_kmer = kmer_set.encoder().encode_kmer(b"NN");
    let mut filtered_kmer_set = ThreeBitKmerSet::<2>::new(2).unwrap_or_die("Cannot build k-mer counter.");
    for kmer in kmer_set.iter_encoded() {
        if *kmer != bad_kmer {
            filtered_kmer_set.insert_encoded_kmer(*kmer);
        }
    }

    // Check to see if it worked
    for kmer in filtered_kmer_set {
        println!("{kmer}")
    }
}
