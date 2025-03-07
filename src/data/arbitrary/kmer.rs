use arbitrary::{Arbitrary, Result, Unstructured};

use crate::kmer::Kmer;

use super::{GraphicAsciiByte, VecBounded};

impl<'a, const MAX_LEN: usize> Arbitrary<'a> for Kmer<MAX_LEN> {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let kmer = VecBounded::<2, MAX_LEN, GraphicAsciiByte>::arbitrary(u)?
            .0
            .into_iter()
            .map(|x| x.0)
            .collect::<Vec<_>>();

        Ok(Kmer::new(kmer))
    }
}
