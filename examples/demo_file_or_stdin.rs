use std::{
    env,
    fs::File,
    io::{Read, stdin},
};
use zoe::{
    define_whichever,
    prelude::{FastQReader, Len},
};

// An example of reading FASTQ data from either a file or stdin, and tallying
// the total number of bases. Zoe provides a macro `define_whichever` to make
// this easier.

define_whichever! {
    // This is similar to the excellent `Either` crate, but here you can define
    // as many variants as you like.
    #[doc = "An enum allowing either reading from a file or from stdin"]
    pub(crate) enum ReadFileStdin {
        File(std::fs::File),
        Stdin(std::io::Stdin),
    }

    // Methods are implemented automatically by define_whichever
    impl Read for ReadFileStdin {}
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();

    let reader = if args.len() == 2 {
        FastQReader::new(ReadFileStdin::File(File::open(args[1].clone())?))
    } else {
        FastQReader::new(ReadFileStdin::Stdin(stdin()))
    };

    let mut total_bases = 0;
    for record in reader {
        total_bases += record?.sequence.len();
    }

    println!("Total bases: {total_bases}");

    Ok(())
}
