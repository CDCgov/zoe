
use irma::data::types::cigar::*;

fn main() {
    let cigar: Cigar = "3S10M2I2D3M4H4P".into();
    let length  = cigar.match_length();
    println!("{length} for {cigar}");

    println!("\n impl Into Iterator method:");
    for (inc,op) in cigar.into_iter_tuple() {
        println!("({},{})", inc, op as char);
    }

    println!("\n Ciglets:");
    for cig in cigar.into_iter() {
        println!("{cig}");
    }
}
