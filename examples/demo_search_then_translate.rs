use std::env;
use zoe::{kmer::ThreeBitKmerSet, prelude::*};

// Search a sequence for the first start codon and the subsequent stop codon.
// Translate the coding region between them, and also calculate the average
// quality.

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example_search_then_translate.fastq".to_owned()
    };

    let start_codon: Nucleotides = b"ATG".into();
    let stop_codons = {
        let mut out = ThreeBitKmerSet::<3>::new(3).unwrap();
        out.insert_kmer(b"TAG");
        out.insert_kmer(b"TGA");
        out.insert_kmer(b"TAA");
        out
    };

    let iterator = FastQReader::from_filename(filename)
        .unwrap_or_die("Translate file error!")
        .flatten();

    for fastq_record in iterator {
        let Some(start_codon_position) = fastq_record.sequence.find_substring(&start_codon) else {
            continue;
        };
        let fastq_record = fastq_record.slice(start_codon_position.end..);

        let Some(stop_codon_position) = stop_codons.find_kmers(&fastq_record.sequence) else {
            continue;
        };
        let fastq_record = fastq_record.slice(..stop_codon_position.start);

        let header = fastq_record.header;
        let protein = fastq_record.translate();
        let Some(average_quality) = fastq_record.quality.arithmetic_mean().map(|x| x.as_f32()) else {
            continue;
        };

        println!("{header}: {protein} with average quality {average_quality}");
    }
}
