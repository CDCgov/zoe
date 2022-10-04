use irma::data::fasta::*;
use std::env;

fn main() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();

    let filename = if args.len() == 2 {
        args[1].clone()
    } else {
        "examples/example.fasta".to_owned()
    };

    let fasta_reader = FastaReader::from_filename(&filename)?;

    for mut record in fasta_reader {
        record.reverse_complement();
        print!("{record}");
    }  
    Ok(()) 
}
