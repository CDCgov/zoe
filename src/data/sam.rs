use crate::data::types::cigar::Cigar;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct SamRow {
    pub qname:          String,
    pub flag:           String,
    pub reference_name: String,
    pub position:       usize,
    pub mapq:           String,
    pub cigar:          Cigar,
    pub mrnm:           String,
    pub mpos:           String,
    pub index_size:     String,
    pub seq:            Vec<u8>,
    pub qual:           Vec<u8>,
}

#[derive(Debug)]
pub struct SamReadSlice {
    pub bases:     Vec<u8>,
    pub qualities: Vec<u8>,
}

#[derive(Debug)]
pub struct SamRowAligned {
    pub qname:          String,
    pub flag:           String,
    pub reference_name: String,
    pub position:       usize,
    pub mapq:           String,
    pub cigar:          Cigar,
    pub mrnm:           String,
    pub mpos:           String,
    pub index_size:     String,
    pub seq:            Vec<u8>,
    pub qual:           Vec<u8>,
    pub aligned:        Vec<u8>,
    pub qaligned:       Vec<u8>,
    pub insertions:     HashMap<usize, SamReadSlice>,
}

impl SamRowAligned {
    fn new(row: SamRow, aligned: Vec<u8>, qaligned: Vec<u8>, insertions: HashMap<usize, SamReadSlice>) -> SamRowAligned {
        SamRowAligned {
            qname: row.qname,
            flag: row.flag,
            reference_name: row.reference_name,
            position: row.position,
            mapq: row.mapq,
            cigar: row.cigar,
            mpos: row.mpos,
            mrnm: row.mrnm,
            index_size: row.index_size,
            seq: row.seq,
            qual: row.qual,
            aligned,
            qaligned,
            insertions,
        }
    }
}

impl From<SamRow> for SamRowAligned {
    #![allow(non_snake_case)]
    fn from(row: SamRow) -> Self {
        let mut rpos = row.position - 1;
        let mut qpos = 0;

        let mut aln: Vec<u8> = Vec::new();
        let mut qAln: Vec<u8> = Vec::new();
        let mut insertions: HashMap<usize, SamReadSlice> = HashMap::new();

        for (inc, op) in row.cigar.into_iter_tuple() {
            match op {
                b'M' => {
                    for _ in 0..inc {
                        qAln.push(row.qual[qpos]);
                        aln.push(row.seq[qpos]);
                        qpos += 1;
                        rpos += 1;
                    }
                }
                b'D' => {
                    for _ in 0..inc {
                        qAln.push(b' ');
                        aln.push(b'-');
                    }
                    rpos += inc;
                }
                b'I' => {
                    qpos += inc;
                    let ins = SamReadSlice {
                        bases:     row.seq[qpos..(qpos + inc)].to_owned(),
                        qualities: row.qual[qpos..(qpos + inc)].to_owned(),
                    };
                    insertions.insert(rpos - 1, ins);
                }

                // Rather than hashing into it, let's add it as a part of the row struct
                //$insByIndex{$K}{$rpos-1} =
                // Rather than an array we want a struct or tuple
                // Something like samInsert.sequence .qualities
                // [substr($seq,$qpos,$inc),substr($qual,$qpos,$inc)];
                b'S' => qpos += inc,
                b'N' => {
                    for _ in 0..inc {
                        qAln.push(b' ');
                        aln.push(b'N');
                    }
                    rpos += inc;
                }
                b'H' => continue,
                _ => panic!("Extended CIGAR {op} not yet supported.\n"),
            }
        }

        SamRowAligned::new(row, aln, qAln, insertions)
        //$pairs{$qMolID}{$qSide} = [$aln,$qAln,$K,($pos-1),($rpos-1),$qname,$mapq];
    }
}
