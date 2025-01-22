use std::{env, ops::Range};
use zoe::prelude::*;

// Given a FASTQ file containing sequences, calculate the percent GC content of
// the coding regions (anywhere between a start and stop codon, or after a start
// codon to the end of the sequence) and the noncoding regions, exclude the
// start and stop codons themselves. Assumes T is used rather than U.

const STOP_CODONS: [&[u8]; 3] = [b"TAA", b"TAG", b"TGA"];

fn find_stop_codon(sequence: &NucleotidesView, starting: usize) -> Option<Range<usize>> {
    sequence[starting..]
        .chunks_exact(3)
        .position(|x| STOP_CODONS.contains(&x))
        .map(|num_codons_before_stop| {
            let stop_codon_start = starting + num_codons_before_stop * 3;
            stop_codon_start..stop_codon_start + 3
        })
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fastq".to_owned()
    };

    let mut gc_coding = 0;
    let mut total_coding = 0;
    let mut gc_noncoding = 0;
    let mut total_noncoding = 0;

    let iterator = FastQReader::from_filename(filename)
        .unwrap_or_die("Translate file error!")
        .flatten();

    for mut fastq_record in iterator {
        // Create a view, whose bounds we will change to reflect the portion of
        // the sequence we are considering (this avoids needing to have an
        // indexing variable)
        fastq_record.sequence.recode_iupac_to_actgn_uc();
        let mut sequence = fastq_record.sequence.as_view();

        // Each execution of the loop represents one coding region
        while let Some(start_codon_position) = sequence.find_substring(b"ATG") {
            // Tally sequence before start codon
            let before_coding = sequence.slice(..start_codon_position.start);
            gc_noncoding += before_coding.gc_content();
            total_noncoding += before_coding.len();

            // Find the stop codon position, respecting the reading frame, then
            // tally the coding region
            if let Some(stop_codon_position) = find_stop_codon(&sequence, start_codon_position.end) {
                let coding = sequence.slice(start_codon_position.end..stop_codon_position.start);
                gc_coding += coding.gc_content();
                total_coding += coding.len();
                // The view is restricted to only hold the unprocessed portion
                sequence.restrict(stop_codon_position.end..);
            } else {
                let coding = sequence.slice(start_codon_position.end..);
                gc_coding += coding.gc_content();
                total_coding += coding.len();
                // The view is restricted to be empty
                sequence.restrict(sequence.len()..);
            }
        }

        // Process any remaining bases after the last coding region
        gc_noncoding += sequence.gc_content();
        total_noncoding += sequence.len();
    }

    let percent_gc_coding = (gc_coding as f32) / (total_coding as f32) * 100.0;
    let percent_gc_noncoding = (gc_noncoding as f32) / (total_noncoding as f32) * 100.0;

    println!("Percent GC Coding: {percent_gc_coding}");
    println!("Percent GC Non-Coding: {percent_gc_noncoding}");
}
