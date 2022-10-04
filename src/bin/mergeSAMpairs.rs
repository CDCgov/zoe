#![allow(non_snake_case,
         unused_variables,
         unused_mut,
         dead_code,
         unused_assignments)]

const PROGRAM: &'static str = "mergeSAMpairs";

use irma::data::{fasta::*, sam::*};

use lazy_static::lazy_static;
use regex::Regex;
//use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;


fn finish(message: &str) -> ! {
    eprintln!("\n{message}\n");
    std::process::exit(0);
}

fn die(message: &str) -> ! {
    eprintln!("\n{PROGRAM} ERROR! {message}\n");
    std::process::exit(1);
}

// Courtesy: https://doc.rust-lang.org/rust-by-example/std_misc/file/read_lines.html
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn get_first_fasta_record(filename: &str) -> Option<FastaSeq> {
    let file = match File::open(filename) {
        Err(why) => die(&format! {"Couldn't open {filename}: {why}"}),
        Ok(file) => file,
    };

    let reader = FastaReader::new(file);

    for record in reader.into_iter() {
        return Some(record);
    }    
    /* *
    // Can haz buffered reader
    let mut b_reader = BufReader::new(file);
    let mut buf = vec![];

    // For writing
    //let stdout = std::io::stdout();
    //let mut w_buf = std::io::BufWriter::with_capacity(32 * 1024, stdout.lock());

    // Get FASTA data
    loop {
        // Look for header row
        let bytes = match b_reader.read_until(b'>', &mut buf) {
            Err(why) => die(&format! {"\nCouldn't read from file: {why}"}),
            Ok(bytes) => bytes,
        };

        if bytes == 0 {
            break;
        } else if buf.len() > 1 {
            // Inclusive split
            if buf.ends_with(b">") {
                buf.pop();
            }

            let mut lines = buf.split(|x| *x == b'\n' || *x == b'\r');
            let name = match lines.next() {
                Some(t) => t.to_vec(),
                None => b"UNKNOWN".to_vec(),
            };

            let sequence: Vec<u8> = lines.flatten().copied().collect();

            if sequence.len() > 0 {
                return Some(FastaSeq { name, sequence });
            }
        }
        buf.clear();
    }
    */

    return None;
}

fn main() {
    let args: Vec<String> = env::args().collect();


    if args.len() != 4 {
        finish(&format! {"Usage:\n\t{PROGRAM} <ref> <sam> <prefix> [-S|--use-serde] [-B|--bowtie-format]]\n"});
    }

    let reference_file = args[1].clone();
    let reference = get_first_fasta_record(&reference_file)
        .expect(&format! {"\n{PROGRAM} ERROR! See: {reference_file}"});
    println!("{reference}");

    let sam_file = args[2].clone();
    let mut sam_data: Vec<SamRowAligned> = Vec::new();
    if let Ok(lines) = read_lines(&sam_file) {
        for line in lines.filter_map(Result::ok) {
            if !line.starts_with('@') {
                let r: Vec<&str> = line.split('\t').collect();
                if r.len() >= 11 && r[2].as_bytes() == reference.name {
                    let row = SamRow { qname: r[0].to_owned(),
                                       flag: r[1].to_owned(),
                                       reference_name: r[2].to_owned(),
                                       position: r[3].parse().unwrap(),
                                       mapq: r[4].to_owned(),
                                       cigar: r[5].into(),
                                       mrnm: r[6].to_owned(),
                                       mpos: r[7].to_owned(),
                                       index_size: r[8].to_owned(),
                                       seq: r[9].as_bytes().to_owned(),
                                       qual: r[10].as_bytes().to_owned() };

                    lazy_static! {
                        static ref MOL_NAME_ID: Regex = Regex::new(r"?-u:.+?)[_ ]([12]):.+")
                            .expect("REGEX molecular ID didn't compile.");
                    }

                    let r_start: usize = row.position - 1;
                    let r_end: usize = r_start + row.cigar.match_length() - 1;

                    // bowtie is order based
                    let (q_name, q_side) =
                        if let Some(cap) = MOL_NAME_ID.captures_iter(&row.reference_name).next() {
                            (cap[1].to_owned(), cap[2].to_owned())
                        } else {
                            (row.reference_name.clone(), "1".to_owned())
                        };

                    sam_data.push(row.into());
                }
            } else {
                // REDIRECT to OSAM, arg[2]
                // open($OSAM,'>',$ARGV[2].'.sam') or die("Cannot open $ARGV[2].sam\n");
                println!("{line}");
            }
        }
    } else {
        die(&format! {"SAM file could not be opened: {sam_file}"});
    }
}

#[inline(always)]
fn is_base(b: u8) -> bool {
    b == b'A' || b == b'G' || b == b'T' || b == b'C'
}

/*my %pairs = ();
my %insByIndex = ();
foreach my $K ( 0 .. $#sam ) {

    my ($qMolID, $qSide);

    if ( $bowtieFormat ) {
        $qMolID = $qname;
        if ( defined($pairs{$qMolID}) ) {
            $qSide = 2;
        } else {
            $qSide = 1;
        }
    } elsif ( $qname =~ $REgetMolID ) {
        $qMolID = $1;
        $qSide = $2;
    }

#	if ( $REF_NAME eq $rn ) {
*/

/*
my %observations = ();
my $OSAM;

my ($dmv, $obs, $fmv, $tmv) = (0,0,0,0);
my ($insObs, $insErr) = (0,0);
foreach my $mID ( keys(%pairs) ) {
    my @mPairs = keys(%{$pairs{$mID}});
    if ( scalar(@mPairs) == 2 ) {
        my @a1 = @{$pairs{$mID}{'1'}};
        my @a2 = @{$pairs{$mID}{'2'}};
        my ($s1,$e1) = ($a1[3],$a1[4]);
        my ($s2,$e2) = ($a2[3],$a2[4]);

        # algorithm doesn't require full alignment, just contained region
        my $start = min($s1,$s2);
        my $end = max($e1,$e2);

        my $mSeq = '';
        my $cigars = '';
        my $qSeq = '';

        my $K1 = $a1[2];
        my $K2 = $a2[2];
        my @bases1 = split('',$a1[0]);
        my @bases2 = split('',$a2[0]);
        my @quals1 = unpack("c* i*",$a1[1]);
        my @quals2 = unpack("c* i*",$a2[1]);

        foreach my $i ( $start .. $end ) {
            # offset here
            # i (ref space) to match states of query q
            # if i < rq_start || i > rq_end return '.'
            # else return ri - OFFSET (qi)
            # OFFSET = rq_start
            # when rq_start = ri then qi = 0
            # when ri = rq_start + 1 then qi = 1
            # qi = (ri - qi_start)
            my $x = $bases1[$i];
            my $y = $bases2[$i];
            my $qx = $quals1[$i];
            my $qy = $quals2[$i];
            my $r = $REF_SEQ[$i];

            # overlap
            if ( $x ne '.' && $y ne '.' ) {
                $obs++;
                # no correction needed
                if ( $x eq $y ) {
                    if ( $x eq '-' ) {
                        $tmv++;
                        $cigars .= 'D';
                    } else {
                        if ( $x ne $r ) { $tmv++; }
                        $cigars .= 'M';
                        $mSeq .= $x;
                        $qSeq .= chr(max($qx,$qy));
                    }
                # reference corrected to r, x
                } elsif( $x eq $r ) {
                    $fmv++;
                    if ( $y eq '-' ) { $dmv++; }
                    $mSeq .= $x;
                    $qSeq .= chr($qx);
                    $cigars .= 'M';
                # reference corrected to r, y
                } elsif( $y eq $r ) {
                    $fmv++;
                    if ( $x eq '-' ) { $dmv++; }
                    $mSeq .= $y;
                    $qSeq .= chr($qy);
                    $cigars .= 'M';
                } else {
                    # reference cannot correct
                    $fmv++;
                    # x is canonical, r, x is not
                    if ( $x =~ $REisBase && $y !~ $REisBase ) {
                        $cigars .= 'M';
                        $mSeq .= $x;
                        $qSeq .= chr($qx);
                        if ( $y eq '-' ) { $dmv++; }
                    # y is canonical, r, x is not
                    } elsif ( $x !~ $REisBase && $y =~ $REisBase ) {
                        $cigars .= 'M';
                        $mSeq .= $y;
                        $qSeq .= chr($qy);
                        # x is also a deletion
                        if ( $x eq '-' ) { $dmv++; }
                    # x more likely
                    } elsif ( $qx > ($qy+4) ) {
                        $cigars .= 'M';
                        $mSeq .= $x;
                        $qSeq .= chr($qx);
                        if ( $y eq '-' ) { $dmv++; }
                    # y more likely
                    } elsif ( $qy > ($qx+4) ) {
                        $cigars .= 'M';
                        $mSeq .= $y;
                        $qSeq .= chr($qy);
                        if ( $x eq '-' ) { $dmv++; }
                    # disagreement, both bases, no information
                    } else {
                        $cigars .= 'M';
                        $mSeq .= 'N';
                        $qSeq .= chr(int(avg($qx,$qy)));
                    }
                }
            # y exists, use that
            } elsif ( $x eq '.' && $y ne '.' ) {
                if ( $y eq '-' ) {
                    $cigars .= 'D';
                } else {
                    $cigars .= 'M';
                    $mSeq .= $y;
                    $qSeq .= chr($qy);
                }
            # x exists, use that
            } elsif( $x ne '.' && $y eq '.' ) {
                if ( $x eq '-' ) {
                    $cigars .= 'D';
                } else {
                    $cigars .= 'M';
                    $mSeq .= $x;
                    $qSeq .= chr($qx);
                }
            # neither exist, so pad
            } else {
                $cigars .= 'N';
            }

            if ( defined($insByIndex{$K1}{$i}) && defined($insByIndex{$K2}{$i}) ) {
                my $ins1 = lc($insByIndex{$K1}{$i}[0]);
                my $ins2 = lc($insByIndex{$K2}{$i}[0]);
                $insObs++;
                if ( $ins1 eq $ins2 ) {
                    $mSeq .= $ins1;
                    my @qIns1 = split('',$insByIndex{$K1}{$i}[1]);
                    my @qIns2 = split('',$insByIndex{$K2}{$i}[1]);
                    my $qSeqNew = '';

                    # new quality is max quality
                    foreach my $qIndex ( 0 .. (length($ins1)-1) ) {
                        $qSeqNew .= chr(max(ord($qIns1[$qIndex]),ord($qIns2[$qIndex])));
                    }
                    $qSeq .= $qSeqNew;
                    $cigars .= 'I' x length($ins1);
                } elsif ( $ins2 =~ /$ins1/ ) {
                    # 1 in 2
                    $mSeq .= $ins1;
                    $qSeq .= $insByIndex{$K1}{$i}[1];
                    $cigars .= 'I' x length($ins1);
                    $insErr++;
                } elsif ( $ins1 =~ /$ins2/ ) {
                    # 2 in 1
                    $mSeq .= $ins2;
                    $qSeq .= $insByIndex{$K2}{$i}[1];
                    $cigars .= 'I' x length($ins2);
                    $insErr++;
                } else {
                    # total disagreement, so don't use it at all
                    $insErr++;
                }
            } elsif ( defined($insByIndex{$K1}{$i}) ) {
                my $w = '';
                if ( $i != $end ) {
                    $w = $bases2[$i+1]
                } else {
                    $w = '.';
                }

                # TO-DO: can ssw permit hanging insertions?
                if ( $y ne '.' && $w ne '.' ) {
                    $insObs++; $insErr++;
                } else {
                    my $ins1 = lc($insByIndex{$K1}{$i}[0]);

                    $mSeq .= $ins1;
                    $qSeq .= $insByIndex{$K1}{$i}[1];
                    $cigars .= 'I' x length($ins1);
                }
            } elsif ( defined($insByIndex{$K2}{$i}) ) {
                my $v = '';
                if ( $i != $end ) {
                    $v = $bases1[$i+1]
                } else {
                    $v = '.';
                }

                # if other seq i and i+1 are not null, then we have a mismatch
                if ( $x ne '.' && $v ne '.' ) {
                    $insObs++; $insErr++;
                # otherwise we can use the insertion from K2
                } else {
                    my $ins2 = lc($insByIndex{$K2}{$i}[0]);

                    $mSeq .= $ins2;
                    $qSeq .= $insByIndex{$K2}{$i}[1];
                    $cigars .= 'I' x length($ins2);
                }
            }
        }

        my $qname = $a1[5];
        my $mapq = int(avg($a1[6],$a2[6]));

        if ( ! $bowtieFormat ) {
            $qname =~ s/(.+?[_ ])[12](:.+)/${1}3${2}/;
        }
        print $OSAM $qname,"\t",'0',"\t",$REF_NAME,"\t",($start+1),"\t",$mapq;
        print $OSAM "\t",condenseCigar($cigars),"\t*\t0\t0\t",$mSeq,"\t",$qSeq,"\n";
    } else {
        my $K = $pairs{$mID}{$mPairs[0]}[2];
        print $OSAM $sam[$K],"\n";
    }
}

$observations{$REF_NAME}{'obs'} = $obs;
$observations{$REF_NAME}{'fmv'} = $fmv;
$observations{$REF_NAME}{'tmv'} = $tmv;
$observations{$REF_NAME}{'dmv'} = $dmv;
$observations{$REF_NAME}{'insObs'} = $insObs;
$observations{$REF_NAME}{'insErr'} = $insErr;

close($OSAM);
if ( $useStorable ) {
    store(\%observations,$ARGV[2].'.sto');
}
*/
