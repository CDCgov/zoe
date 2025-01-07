#![feature(array_chunks)]

use std::{collections::HashSet, env};
use zoe::{composition::GcContent, prelude::*};

// Given a FastQ file containing sequences with a single coding region,
// calculate the percent GC content of the coding regions and the noncoding
// regions, exclude start and stop codons.

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example_gc_coding_noncoding.fastq".to_owned()
    };

    let start_codon: Nucleotides = b"ATG".into();
    let stop_codons = {
        let mut out = HashSet::new();
        out.insert(b"TAG");
        out.insert(b"TGA");
        out.insert(b"TAA");
        out
    };

    let mut gc_coding = 0;
    let mut total_coding = 0;
    let mut gc_noncoding = 0;
    let mut total_noncoding = 0;

    let iterator = FastQReader::from_filename(filename)
        .unwrap_or_die("Translate file error!")
        .flatten();

    for fastq_record in iterator {
        let sequence_view = fastq_record.sequence;

        let Some(start_codon_position) = sequence_view.find_substring(&start_codon) else {
            continue;
        };

        let Some(num_codons_before_stop) = sequence_view
            .slice(start_codon_position.end..)
            .as_bytes()
            .array_chunks::<3>()
            .position(|x| stop_codons.contains(x))
        else {
            continue;
        };
        let stop_codon_start = start_codon_position.end + num_codons_before_stop * 3;
        let stop_codon_position = stop_codon_start..stop_codon_start + 3;

        let before_coding = sequence_view.slice(..start_codon_position.start);
        let coding = sequence_view.slice(start_codon_position.end..stop_codon_position.start);
        let after_coding = sequence_view.slice(stop_codon_position.end..);

        gc_coding += coding.gc_content();
        total_coding += coding.len();

        gc_noncoding += before_coding.gc_content() + after_coding.gc_content();
        total_noncoding += before_coding.len() + after_coding.len();
    }

    let percent_gc_coding = (gc_coding as f32) / (total_coding as f32) * 100.0;
    let percent_gc_noncoding = (gc_noncoding as f32) / (total_noncoding as f32) * 100.0;

    println!("Percent GC Coding: {percent_gc_coding}");
    println!("Percent GC Non-Coding: {percent_gc_noncoding}");
}
