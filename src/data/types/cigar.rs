use lazy_static::lazy_static;
use regex::bytes::Regex;
use atoi::atoi;


#[derive(Debug, Clone)]
pub struct Cigar(Vec<u8>);

impl Cigar {
    pub fn expand_cigar(&self) -> ExpandedCigar {
        lazy_static! {
            static ref CIGAR_PARSE: Regex =
                Regex::new(r"?-u:(\d+)([MIDNSHP])").expect("REGEX cigar ID didn't compile.");
        }

        let caps = CIGAR_PARSE.captures_iter(self.0.as_slice());
        let mut expanded = Vec::new();
        for cap in caps {
            let inc: usize = atoi(&cap[1]).unwrap();
            let op = cap[2][0];
            for _ in 0..inc {
                expanded.push(op);
            }
        }
        return ExpandedCigar(expanded);
    }
    pub fn match_length(&self) -> usize {
        lazy_static! {
            static ref CIGAR_MATCH_STATE: Regex =
                Regex::new(r"(\d+)[MDN]").expect("REGEX match_length didn't compile.");
        }
        let mut length: usize = 0;
        for caps in CIGAR_MATCH_STATE.captures_iter(self.0.as_slice()) {
            length += atoi::<usize>(&caps[1]).unwrap();
        }
    
        return length;
    }

    pub fn into_iter_tuple(&self) -> impl Iterator<Item=(usize,u8)> + '_ {
        lazy_static! {
            static ref CIGAR_PARSE: Regex =
                Regex::new(r"(\d+)([MIDNSHP])").expect("REGEX cigar ID didn't compile.");
        }
        let caps = CIGAR_PARSE.captures_iter(self.0.as_slice());

        return caps.map( |c| (atoi::<usize>(&c[1]).unwrap(), c[2][0]));
    }

    pub fn into_iter(&self) -> impl Iterator<Item=Ciglet> + '_ {
        lazy_static! {
            static ref CIGAR_PARSE: Regex =
                Regex::new(r"(\d+)([MIDNSHP])").expect("REGEX cigar ID didn't compile.");
        }
        let caps = CIGAR_PARSE.captures_iter(self.0.as_slice());

        return caps.map( |c| Ciglet{inc: atoi::<usize>(&c[1]).unwrap(), op: c[2][0] } );
    }
}

use std::fmt;
impl fmt::Display for Cigar {
   /*  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for v in &self.0 {
            write!(f, "\t{}", v)?;
        }
        Ok(())
    }
    */
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"{}", String::from_utf8_lossy(&self.0))?;
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
         write!(f,"({},{})", self.inc, self.op as char)?;
         Ok(())
     }
 }


/* Infinite loop, need Iter Struct to manage it
impl Iterator for Cigar {
    type Item = (usize, u8);
    fn next(&mut self) -> Option<Self::Item> {
        lazy_static! {
            static ref CIGAR_PARSE: Regex =
                Regex::new(r"(\d+)([MIDNSHP])").expect("REGEX cigar ID didn't compile.");
        }

        let mut caps = CIGAR_PARSE.captures_iter(self.0.as_slice());
        //caps.map( |c| (atoi::<usize>(&c[1]), c[2][0]));
        if let Some(cap) = caps.next() {
            let inc: usize = atoi(&cap[1]).unwrap();
            let op = cap[2][0];
            Some((inc,op))
        } else {
            None
        }
        
    }
}*/

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
    #[inline(always)]
    pub fn condense_cigar(self) -> Cigar {
        lazy_static! {
            static ref CIGAR_EXPANDED: Regex = Regex::new(r"?-u:([M]+|[D]+|[I]+|[H]+|[N]+|[S]+)")
                .expect("REGEX condense_cigar didn't compile.");
        }

        let mut condensed: Vec<u8> = Vec::new();
        for caps in CIGAR_EXPANDED.captures_iter(self.0.as_slice()) {
            condensed.extend_from_slice(caps[1].len().to_string().as_bytes());
            condensed.push(caps[1][0]);
        }

        return Cigar(condensed);
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
