use std::env;
use zoe::{kmer::ThreeBitKmerSet, prelude::*};

// A toy example for performing several subsequent trimming operations cheaply
// using views:
// 1. Remove `G`s from the beginning or end if at least 5 are present in a row
// 2. Hard trim 2 bases from the beginning and end
// 3. Search for a 17-mers of a primer within the first 30 bases, and trim
//    everything before the last such occurrence

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fastq".to_owned()
    };

    let primer = b"TGATAGTTTTAGAGTTAGGTAG";
    let mut kmer_set = ThreeBitKmerSet::<17>::new(17).unwrap();
    kmer_set.insert_from_sequence_one_mismatch(primer);

    let iterator = FastQReader::from_filename(filename)
        .unwrap_or_die("Translate file error!")
        .flatten();

    for fastq_record in iterator {
        let mut fastq_view = fastq_record.as_view();

        let num_g_start = fastq_view.sequence.iter().take_while(|&&x| x == b'G').count();
        if num_g_start >= 5 {
            fastq_view.restrict(num_g_start..);
        }

        let num_g_end = fastq_view.sequence.iter().rev().take_while(|&&x| x == b'G').count();
        if num_g_end >= 5 {
            fastq_view.restrict(..fastq_view.len() - num_g_end);
        }

        fastq_view.restrict(2..fastq_view.len() - 2);

        let search_region_len = fastq_view.len().min(30);
        if let Some(kmer_pos) = kmer_set.find_kmers_rev(fastq_view.sequence.slice(..search_region_len)) {
            fastq_view.restrict(kmer_pos.end..);
        }

        println!("{fastq_view}")
    }
}
