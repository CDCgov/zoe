//! Functionality for analyzing ambiguous IUPAC bases.

/// A look-up table for disambiguating IUPAC nucleotides.
///
/// See the [IUPAC definitions](https://www.bioinformatics.org/sms/iupac.html)
/// for more details.
pub struct DnaDisambiguation;

impl DnaDisambiguation {
    /// The look-up table used for disambiguation.
    const TABLE: &[u8; 256] = &IUPAC_DISAMBIGUATION;

    /// Returns whether the codon potentially codes for a stop codon under the
    /// standard genetic code (`TAA`, `TAG`, or `TGA`).
    ///
    /// Note that `NNN` returns `true`.
    #[must_use]
    pub fn maybe_std_stop_codon(codon: [u8; 3]) -> bool {
        // Check that first base is possibly `T`, second base is possibly `A` or
        // `G`, third base is possibly `A` or `G`, and either the second or
        // third base is possibly `A`.
        Self::maybe_t(codon[0]) && {
            let encoded1 = Self::TABLE[codon[1] as usize];
            let encoded2 = Self::TABLE[codon[2] as usize];
            (encoded1 & (MAYBE_A | MAYBE_G) > 0)
                && (encoded2 & (MAYBE_A | MAYBE_G) > 0)
                && ((encoded1 | encoded2) & MAYBE_A > 0)
        }
    }

    /// Returns whether the `base` potentially represents `A`.
    #[must_use]
    pub fn maybe_a(base: u8) -> bool {
        Self::TABLE[base as usize] & MAYBE_A > 0
    }

    /// Returns whether the `base` potentially represents `C`.
    #[must_use]
    pub fn maybe_c(base: u8) -> bool {
        Self::TABLE[base as usize] & MAYBE_C > 0
    }

    /// Returns whether the `base` potentially represents `G`.
    #[must_use]
    pub fn maybe_g(base: u8) -> bool {
        Self::TABLE[base as usize] & MAYBE_G > 0
    }

    /// Returns whether the `base` potentially represents `T`.
    #[must_use]
    pub fn maybe_t(base: u8) -> bool {
        Self::TABLE[base as usize] & MAYBE_T > 0
    }

    /// Returns whether the `base` potentially represents `A` or `C`.
    #[must_use]
    pub fn maybe_ac(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_A | MAYBE_C) > 0
    }

    /// Returns whether the `base` potentially represents `A` or `G`.
    #[must_use]
    pub fn maybe_ag(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_A | MAYBE_G) > 0
    }

    /// Returns whether the `base` potentially represents `A` or `T`.
    #[must_use]
    pub fn maybe_at(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_A | MAYBE_T) > 0
    }

    /// Returns whether the `base` potentially represents `C` or `G`.
    #[must_use]
    pub fn maybe_cg(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_C | MAYBE_G) > 0
    }

    /// Returns whether the `base` potentially represents `C` or `T`.
    #[must_use]
    pub fn maybe_ct(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_C | MAYBE_T) > 0
    }

    /// Returns whether the `base` potentially represents `G` or `T`.
    #[must_use]
    pub fn maybe_gt(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_G | MAYBE_T) > 0
    }

    /// Returns whether the `base` potentially represents `A`, `C`, or `G`.
    #[must_use]
    pub fn maybe_acg(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_A | MAYBE_C | MAYBE_G) > 0
    }

    /// Returns whether the `base` potentially represents `A`, `C`, or `T`.
    #[must_use]
    pub fn maybe_act(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_A | MAYBE_C | MAYBE_T) > 0
    }

    /// Returns whether the `base` potentially represents `A`, `G`, or `T`.
    #[must_use]
    pub fn maybe_agt(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_A | MAYBE_G | MAYBE_T) > 0
    }

    /// Returns whether the `base` potentially represents `A`, `G`, or `T`.
    #[must_use]
    pub fn maybe_cgt(base: u8) -> bool {
        Self::TABLE[base as usize] & (MAYBE_C | MAYBE_G | MAYBE_T) > 0
    }

    /// Returns whether the `base` is valid IUPAC, and hence represents either
    /// `A`, `C`, `G`, or `T`.
    #[must_use]
    pub fn maybe_acgt(base: u8) -> bool {
        Self::TABLE[base as usize] > 0
    }
}

/// The bit used to represent an IUPAC code potentially representing `A`.
const MAYBE_A: u8 = 0b0001;
/// The bit used to represent an IUPAC code potentially representing `C`.
const MAYBE_C: u8 = 0b0010;
/// The bit used to represent an IUPAC code potentially representing `G`.
const MAYBE_G: u8 = 0b0100;
/// The bit used to represent an IUPAC code potentially representing `T`.
const MAYBE_T: u8 = 0b1000;

/// The look-up table used in [`DnaDisambiguation`], where each index represents
/// a byte and each value has bits set to indicate which bases it might
/// represent.
const IUPAC_DISAMBIGUATION: [u8; 256] = {
    /// A helper function for constructing [`IUPAC_DISAMBIGUATION`], which sets
    /// the value for the uppercase and lowercase `byte`.
    const fn fill(byte: u8, to: &[u8], out: &mut [u8; 256]) {
        let lower = byte.to_ascii_lowercase();
        let upper = byte.to_ascii_uppercase();

        let mut to_encoded = 0;
        let mut i = 0;
        while i < to.len() {
            let bit = match to[i] {
                b'A' => MAYBE_A,
                b'C' => MAYBE_C,
                b'G' => MAYBE_G,
                b'T' => MAYBE_T,
                _ => panic!("Invalid byte found"),
            };
            to_encoded |= bit;
            i += 1;
        }

        out[lower as usize] = to_encoded;
        out[upper as usize] = to_encoded;
    }

    let mut out = [0; 256];

    fill(b'A', b"A", &mut out);
    fill(b'C', b"C", &mut out);
    fill(b'G', b"G", &mut out);
    fill(b'T', b"T", &mut out);
    fill(b'U', b"T", &mut out);
    fill(b'R', b"AG", &mut out);
    fill(b'Y', b"CT", &mut out);
    fill(b'S', b"GC", &mut out);
    fill(b'W', b"AT", &mut out);
    fill(b'K', b"GT", &mut out);
    fill(b'M', b"AC", &mut out);
    fill(b'B', b"CGT", &mut out);
    fill(b'D', b"AGT", &mut out);
    fill(b'H', b"ACT", &mut out);
    fill(b'V', b"ACG", &mut out);
    fill(b'N', b"ACGT", &mut out);

    out
};
