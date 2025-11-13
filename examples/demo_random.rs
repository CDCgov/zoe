// Requires the `random` feature.

use std::env;
use zoe::prelude::*;

fn main() {
    let args: Vec<String> = env::args().collect();

    let (samples, length) = if args.len() == 3 {
        (args[1].parse().unwrap(), args[2].parse().unwrap())
    } else {
        println!("Usage:\n\t{} <samples> <length>\n\nSetting to: 4 4", &args[0]);
        (4, 4)
    };

    const ALPHA: &[u8] = b"AGTCN";
    for seed in 0..samples {
        let s = rand_sequence(ALPHA, length, seed);
        println!("{}", String::from_utf8_lossy(&s));
    }
}
