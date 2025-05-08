use super::impl_deref;
use crate::data::cigar::{Cigar, Ciglet};
use arbitrary::{Arbitrary, Result, Unstructured};
use std::num::NonZeroUsize;

impl<'a> Arbitrary<'a> for Ciglet {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Ciglet {
            inc: usize::arbitrary(u)?,
            op:  u8::arbitrary(u)?,
        })
    }
}

/// A wrapper around [`Ciglet`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates operations in MIDNSHP=X and the increment is nonzero
#[derive(Debug)]
pub struct CigletOp(pub Ciglet);

impl_deref! {CigletOp, Ciglet}

impl<'a> Arbitrary<'a> for CigletOp {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        const ALPHA: &[u8] = b"MIDNSHP=X";
        Ok(CigletOp(Ciglet {
            inc: NonZeroUsize::arbitrary(u)?.get(),
            op:  ALPHA[usize::arbitrary(u)? % ALPHA.len()],
        }))
    }
}

impl<'a> Arbitrary<'a> for Cigar {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Cigar(Vec::<u8>::arbitrary(u)?))
    }
}

/// A wrapper around [`Cigar`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// only generates Cigar string with the following assumptions:
/// * The string consists of alternating numbers and operations
/// * The numbers must be between 1 and `usize::MAX`
/// * The operations must be in MIDNSHP=X
#[derive(Debug)]
pub struct CigarIncOp(pub Cigar);

impl_deref! {CigarIncOp, Cigar}

impl<'a> Arbitrary<'a> for CigarIncOp {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let mut out = Vec::new();
        for ciglet in u.arbitrary_iter::<CigletOp>()?.flatten() {
            out.extend_from_slice(ciglet.inc.to_string().as_bytes());
            out.push(ciglet.op);
        }
        Ok(CigarIncOp(Cigar(out)))
    }
}

/// A wrapper around [`Cigar`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// is similar to [`CigarIncOp`] but may contain leading zeros before numbers.
#[derive(Debug)]
pub struct CigarIncOpLeadingZero(pub Cigar);

impl_deref! {CigarIncOpLeadingZero, Cigar}

impl<'a> Arbitrary<'a> for CigarIncOpLeadingZero {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        struct CigletOpWithLen(CigletOp, usize);

        impl<'a> Arbitrary<'a> for CigletOpWithLen {
            fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
                let ciglet_op = u.arbitrary()?;
                let len = u.arbitrary_len::<u8>()?;
                Ok(CigletOpWithLen(ciglet_op, len))
            }
        }

        let mut out = Vec::new();
        for CigletOpWithLen(ciglet, leading_zeros) in u.arbitrary_iter::<CigletOpWithLen>()?.flatten() {
            out.extend(std::iter::repeat_n(b'0', leading_zeros));
            out.extend_from_slice(ciglet.inc.to_string().as_bytes());
            out.push(ciglet.op);
        }
        Ok(CigarIncOpLeadingZero(Cigar(out)))
    }
}
