
// Nucleotide reverse complement
const fn init_rev_comp() -> [u8; 256] {
    let mut comp = [0u8; 256];
    let mut i : usize = 0;
    while i < 256 {
        comp[i] = i as u8;    
        i += 1;
    }

    let forward = b"gcatrykmbvdhuGCATRYKMBVDHU";
    let reverse = b"cgtayrmkvbhdaCGTAYRMKVBHDA";
    let mut b : usize = 0;
    while b < 26 {
        comp[ forward[b] as usize ] = reverse[b];
        b += 1;
    }

    return comp;
}

pub(crate) const REV_COMP: [u8; 256] = init_rev_comp();
