use lazy_static::lazy_static;
use std::collections::HashMap;

// Nucleotide reverse complement
#[allow(clippy::cast_possible_truncation)]
pub(crate) const REV_COMP: [u8; 256] = {
    let mut comp = [0u8; 256];
    let mut i: usize = 0;
    while i < 256 {
        // usize cannot be greater than bound
        comp[i] = i as u8;
        i += 1;
    }

    let forward = b"gcatrykmbvdhuGCATRYKMBVDHU";
    let reverse = b"cgtayrmkvbhdaCGTAYRMKVBHDA";
    let mut b: usize = 0;
    while b < 26 {
        comp[forward[b] as usize] = reverse[b];
        b += 1;
    }

    comp
};

macro_rules! fill_map {
    ($( $key: expr => $val: expr ),*) => {{
        let mut map: HashMap<Vec<u8>, u8, ahash::RandomState> = HashMap::default();
        $( map.insert($key.to_vec(), $val); )*
        map
   }}
}

lazy_static! {
    pub(crate) static ref GENETIC_CODE: HashMap<Vec<u8>, u8, ahash::RandomState> = fill_map!(
        b"TAA"=>b'*', b"TAG"=>b'*', b"TAR"=>b'*', b"TGA"=>b'*', b"TRA"=>b'*', b"GCA"=>b'A', b"GCB"=>b'A', b"GCC"=>b'A', b"GCD"=>b'A', b"GCG"=>b'A', b"GCH"=>b'A',
        b"GCK"=>b'A', b"GCM"=>b'A', b"GCN"=>b'A', b"GCR"=>b'A', b"GCS"=>b'A', b"GCT"=>b'A', b"GCV"=>b'A', b"GCW"=>b'A', b"GCY"=>b'A', b"TGC"=>b'C', b"TGT"=>b'C',
        b"TGY"=>b'C', b"GAC"=>b'D', b"GAT"=>b'D', b"GAY"=>b'D', b"GAA"=>b'E', b"GAG"=>b'E', b"GAR"=>b'E', b"TTC"=>b'F', b"TTT"=>b'F', b"TTY"=>b'F', b"GGA"=>b'G',
        b"GGB"=>b'G', b"GGC"=>b'G', b"GGD"=>b'G', b"GGG"=>b'G', b"GGH"=>b'G', b"GGK"=>b'G', b"GGM"=>b'G', b"GGN"=>b'G', b"GGR"=>b'G', b"GGS"=>b'G', b"GGT"=>b'G',
        b"GGV"=>b'G', b"GGW"=>b'G', b"GGY"=>b'G', b"CAC"=>b'H', b"CAT"=>b'H', b"CAY"=>b'H', b"ATA"=>b'I', b"ATC"=>b'I', b"ATH"=>b'I', b"ATM"=>b'I', b"ATT"=>b'I',
        b"ATW"=>b'I', b"ATY"=>b'I', b"AAA"=>b'K', b"AAG"=>b'K', b"AAR"=>b'K', b"CTA"=>b'L', b"CTB"=>b'L', b"CTC"=>b'L', b"CTD"=>b'L', b"CTG"=>b'L', b"CTH"=>b'L',
        b"CTK"=>b'L', b"CTM"=>b'L', b"CTN"=>b'L', b"CTR"=>b'L', b"CTS"=>b'L', b"CTT"=>b'L', b"CTV"=>b'L', b"CTW"=>b'L', b"CTY"=>b'L', b"TTA"=>b'L', b"TTG"=>b'L',
        b"TTR"=>b'L', b"YTA"=>b'L', b"YTG"=>b'L', b"YTR"=>b'L', b"ATG"=>b'M', b"AAC"=>b'N', b"AAT"=>b'N', b"AAY"=>b'N', b"CCA"=>b'P', b"CCB"=>b'P', b"CCC"=>b'P',
        b"CCD"=>b'P', b"CCG"=>b'P', b"CCH"=>b'P', b"CCK"=>b'P', b"CCM"=>b'P', b"CCN"=>b'P', b"CCR"=>b'P', b"CCS"=>b'P', b"CCT"=>b'P', b"CCV"=>b'P', b"CCW"=>b'P',
        b"CCY"=>b'P', b"CAA"=>b'Q', b"CAG"=>b'Q', b"CAR"=>b'Q', b"AGA"=>b'R', b"AGG"=>b'R', b"AGR"=>b'R', b"CGA"=>b'R', b"CGB"=>b'R', b"CGC"=>b'R', b"CGD"=>b'R',
        b"CGG"=>b'R', b"CGH"=>b'R', b"CGK"=>b'R', b"CGM"=>b'R', b"CGN"=>b'R', b"CGR"=>b'R', b"CGS"=>b'R', b"CGT"=>b'R', b"CGV"=>b'R', b"CGW"=>b'R', b"CGY"=>b'R',
        b"MGA"=>b'R', b"MGG"=>b'R', b"MGR"=>b'R', b"AGC"=>b'S', b"AGT"=>b'S', b"AGY"=>b'S', b"TCA"=>b'S', b"TCB"=>b'S', b"TCC"=>b'S', b"TCD"=>b'S', b"TCG"=>b'S',
        b"TCH"=>b'S', b"TCK"=>b'S', b"TCM"=>b'S', b"TCN"=>b'S', b"TCR"=>b'S', b"TCS"=>b'S', b"TCT"=>b'S', b"TCV"=>b'S', b"TCW"=>b'S', b"TCY"=>b'S', b"ACA"=>b'T',
        b"ACB"=>b'T', b"ACC"=>b'T', b"ACD"=>b'T', b"ACG"=>b'T', b"ACH"=>b'T', b"ACK"=>b'T', b"ACM"=>b'T', b"ACN"=>b'T', b"ACR"=>b'T', b"ACS"=>b'T', b"ACT"=>b'T',
        b"ACV"=>b'T', b"ACW"=>b'T', b"ACY"=>b'T', b"GTA"=>b'V', b"GTB"=>b'V', b"GTC"=>b'V', b"GTD"=>b'V', b"GTG"=>b'V', b"GTH"=>b'V', b"GTK"=>b'V', b"GTM"=>b'V',
        b"GTN"=>b'V', b"GTR"=>b'V', b"GTS"=>b'V', b"GTT"=>b'V', b"GTV"=>b'V', b"GTW"=>b'V', b"GTY"=>b'V', b"TGG"=>b'W', b"TAC"=>b'Y', b"TAT"=>b'Y', b"TAY"=>b'Y',
        b"..."=>b'.', b"---"=>b'-', b"NNN"=>b'X'
    );
}

pub(crate) static PHYSIOCHEMICAL_FACTORS: [[Option<f32>; 256]; 256] = {
    const AA: [u8; 43] = [
        b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R', b'S', b'T', b'V', b'W',
        b'Y', b'X', b'a', b'c', b'd', b'e', b'f', b'g', b'h', b'i', b'k', b'l', b'm', b'n', b'p', b'q', b'r', b's', b't',
        b'v', b'w', b'y', b'x', b'-',
    ];

    const PCF: [[f64; 5]; 43] = [
        [-0.59, -1.3, -0.73, 1.57, -0.15],
        [-1.34, 0.47, -0.86, -1.02, -0.26],
        [1.05, 0.3, -3.66, -0.26, -3.24],
        [1.36, -1.45, 1.48, 0.11, -0.84],
        [-1.01, -0.59, 1.89, -0.4, 0.41],
        [-0.38, 1.65, 1.33, 1.05, 2.06],
        [0.34, -0.42, -1.67, -1.47, -0.08],
        [-1.24, -0.55, 2.13, 0.39, 0.82],
        [1.83, -0.56, 0.53, -0.28, 1.65],
        [-1.02, -0.99, -1.51, 1.27, -0.91],
        [-0.66, -1.52, 2.22, -1.01, 1.21],
        [0.95, 0.83, 1.3, -0.17, 0.93],
        [0.19, 2.08, -1.63, 0.42, -1.39],
        [0.93, -0.18, -3.01, -0.5, -1.85],
        [1.54, -0.06, 1.5, 0.44, 2.9],
        [-0.23, 1.4, -4.76, 0.67, -2.65],
        [-0.03, 0.33, 2.21, 0.91, 1.31],
        [-1.34, -0.28, -0.54, 1.24, -1.26],
        [-0.6, 0.01, 0.67, -2.13, -0.18],
        [0.26, 0.83, 3.1, -0.84, 1.51],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.59, -1.3, -0.73, 1.57, -0.15],
        [-1.34, 0.47, -0.86, -1.02, -0.26],
        [1.05, 0.3, -3.66, -0.26, -3.24],
        [1.36, -1.45, 1.48, 0.11, -0.84],
        [-1.01, -0.59, 1.89, -0.4, 0.41],
        [-0.38, 1.65, 1.33, 1.05, 2.06],
        [0.34, -0.42, -1.67, -1.47, -0.08],
        [-1.24, -0.55, 2.13, 0.39, 0.82],
        [1.83, -0.56, 0.53, -0.28, 1.65],
        [-1.02, -0.99, -1.51, 1.27, -0.91],
        [-0.66, -1.52, 2.22, -1.01, 1.21],
        [0.95, 0.83, 1.3, -0.17, 0.93],
        [0.19, 2.08, -1.63, 0.42, -1.39],
        [0.93, -0.18, -3.01, -0.5, -1.85],
        [1.54, -0.06, 1.5, 0.44, 2.9],
        [-0.23, 1.4, -4.76, 0.67, -2.65],
        [-0.03, 0.33, 2.21, 0.91, 1.31],
        [-1.34, -0.28, -0.54, 1.24, -1.26],
        [-0.6, 0.01, 0.67, -2.13, -0.18],
        [0.26, 0.83, 3.1, -0.84, 1.51],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let mut pcd = [[None; 256]; 256];

    let mut aa1: usize = 0;
    while aa1 < 43 {
        let mut aa2: usize = 0;
        while aa2 < 43 {
            pcd[AA[aa1] as usize][AA[aa2] as usize] = if AA[aa1].to_ascii_uppercase() == AA[aa2].to_ascii_uppercase() {
                Some(0.0)
            } else {
                let mut d: f64 = 0.0;
                let mut k: usize = 0;
                while k < 5 {
                    d += (PCF[aa1][k] - PCF[aa2][k]) * (PCF[aa1][k] - PCF[aa2][k]);
                    k += 1;
                }
                #[allow(clippy::cast_possible_truncation)]
                Some(crate::math::sqrt_baby(d) as f32)
            };

            aa2 += 1;
        }
        aa1 += 1;
    }

    pcd
};
