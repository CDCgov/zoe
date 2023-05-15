use atoi::atoi;
use lazy_static::lazy_static;
use regex::bytes::Regex;

#[derive(Debug, Clone)]
pub struct Cigar(Vec<u8>);

impl Cigar {
    #[must_use]
    pub fn expand_cigar(&self) -> ExpandedCigar {
        lazy_static! {
            static ref CIGAR_PARSE: Regex =
                Regex::new(r"(\d+)([MIDNSHP])").expect("REGEX cigar ID didn't compile.");
        }

        let caps = CIGAR_PARSE.captures_iter(self.0.as_slice());
        let mut expanded = Vec::new();
        for cap in caps {
            let inc: usize = if let Some(number) = atoi(&cap[1]) {
                number
            } else {
                continue;
            };
            let op = cap[2][0];
            for _ in 0..inc {
                expanded.push(op);
            }
        }
        ExpandedCigar(expanded)
    }

    #[must_use]
    pub fn match_length(&self) -> usize {
        lazy_static! {
            static ref CIGAR_MATCH_STATE: Regex =
                Regex::new(r"(\d+)[MDN]").expect("REGEX match_length didn't compile.");
        }
        let mut length: usize = 0;
        for caps in CIGAR_MATCH_STATE.captures_iter(self.0.as_slice()) {
            // Valid regex should never be default
            length += atoi::<usize>(&caps[1]).unwrap_or_default();
        }

        length
    }

    pub fn into_iter_tuple(&self) -> impl Iterator<Item = (usize, u8)> + '_ {
        lazy_static! {
            static ref CIGAR_PARSE: Regex =
                Regex::new(r"(\d+)([MIDNSHP])").expect("REGEX cigar ID didn't compile.");
        }
        let caps = CIGAR_PARSE.captures_iter(self.0.as_slice());

        caps.map(|c| (atoi::<usize>(&c[1]).unwrap_or_default(), c[2][0]))
    }

    pub fn into_iter(&self) -> impl Iterator<Item = Ciglet> + '_ {
        lazy_static! {
            static ref CIGAR_PARSE: Regex =
                Regex::new(r"(\d+)([MIDNSHP])").expect("REGEX cigar ID didn't compile.");
        }
        let caps = CIGAR_PARSE.captures_iter(self.0.as_slice());

        caps.map(|c| Ciglet {
            inc: atoi::<usize>(&c[1]).unwrap_or_default(),
            op: c[2][0],
        })
    }
}

use std::fmt;
impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))?;
        Ok(())
    }
}

/// A single increment quantifier-operation pair.
/// Don't hate the name, it could have been Cygnet or Cigarette!

#[derive(Clone, Copy)]
pub struct Ciglet {
    pub inc: usize,
    pub op: u8,
}

impl fmt::Display for Ciglet {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({},{})", self.inc, self.op as char)?;
        Ok(())
    }
}

impl From<&str> for Cigar {
    fn from(s: &str) -> Self {
        Cigar(s.as_bytes().to_owned())
    }
}

impl From<&[u8]> for Cigar {
    fn from(v: &[u8]) -> Self {
        Cigar(v.to_vec())
    }
}

impl From<ExpandedCigar> for Cigar {
    fn from(c: ExpandedCigar) -> Cigar {
        c.condense_cigar()
    }
}

#[derive(Debug)]
pub struct ExpandedCigar(Vec<u8>);

impl ExpandedCigar {
    #[inline]
    #[must_use]
    pub fn condense_cigar(self) -> Cigar {
        lazy_static! {
            static ref CIGAR_EXPANDED: Regex = Regex::new(r"([M]+|[D]+|[I]+|[H]+|[N]+|[S]+)")
                .expect("REGEX condense_cigar didn't compile.");
        }

        let mut condensed: Vec<u8> = Vec::new();
        for caps in CIGAR_EXPANDED.captures_iter(self.0.as_slice()) {
            condensed.extend_from_slice(caps[1].len().to_string().as_bytes());
            condensed.push(caps[1][0]);
        }

        Cigar(condensed)
    }
}

impl From<Cigar> for ExpandedCigar {
    fn from(c: Cigar) -> ExpandedCigar {
        c.expand_cigar()
    }
}

impl From<&str> for ExpandedCigar {
    fn from(s: &str) -> Self {
        ExpandedCigar(s.as_bytes().to_owned())
    }
}

impl From<&[u8]> for ExpandedCigar {
    fn from(v: &[u8]) -> Self {
        ExpandedCigar(v.to_vec())
    }
}
