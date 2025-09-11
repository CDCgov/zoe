use crate::{
    alignment::{AlignmentStates, CheckedCigar},
    data::{
        arbitrary::impl_deref,
        cigar::{Cigar, Ciglet},
    },
};
use arbitrary::{Arbitrary, Result, Unstructured};
use std::{marker::PhantomData, num::NonZeroUsize};

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

impl From<CigletOp> for Ciglet {
    #[inline]
    fn from(value: CigletOp) -> Self {
        value.0
    }
}

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

impl<'a> Arbitrary<'a> for AlignmentStates {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(AlignmentStates(Vec::<Ciglet>::arbitrary(u)?))
    }
}

/// A wrapper around [`Cigar`] such that the implementation of
/// [`Arbitrary`](https://docs.rs/arbitrary/latest/arbitrary/trait.Arbitrary.html)
/// is generated from [`Ciglet`] structs.
///
/// The type of the ciglet is specified with `C`, which could be [`Ciglet`] or
/// [`CigletOp`]. If `Z` is true, then leading zeros are added to the start of
/// some increments in the CIGAR string.
#[derive(Debug)]
pub struct CigarArbitrary<C, const Z: bool>(pub Cigar, PhantomData<C>);

impl_deref! {CigarArbitrary<C, Z>, Cigar, <C, const Z: bool>}

impl<'a, C> Arbitrary<'a> for CigarArbitrary<C, false>
where
    C: Arbitrary<'a> + Into<Ciglet>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let mut out = Vec::new();
        for ciglet in u.arbitrary_iter::<C>()?.flatten() {
            let ciglet = ciglet.into();
            out.extend_from_slice(ciglet.inc.to_string().as_bytes());
            out.push(ciglet.op);
        }
        Ok(CigarArbitrary(Cigar(out), PhantomData))
    }
}

impl<'a, C> Arbitrary<'a> for CigarArbitrary<C, true>
where
    C: Arbitrary<'a> + Into<Ciglet>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        struct CigletOpWithLen<C>(C, usize);

        impl<'a, C> Arbitrary<'a> for CigletOpWithLen<C>
        where
            C: Arbitrary<'a> + Into<Ciglet>,
        {
            fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
                let ciglet_op = u.arbitrary()?;
                let len = u.arbitrary_len::<u8>()?;
                Ok(CigletOpWithLen(ciglet_op, len))
            }
        }

        let mut out = Vec::new();
        for CigletOpWithLen(ciglet, leading_zeros) in u.arbitrary_iter::<CigletOpWithLen<C>>()?.flatten() {
            out.extend(std::iter::repeat_n(b'0', leading_zeros));
            let ciglet = ciglet.into();
            out.extend_from_slice(ciglet.inc.to_string().as_bytes());
            out.push(ciglet.op);
        }
        Ok(CigarArbitrary(Cigar(out), PhantomData))
    }
}

/// A wrapper around [`AlignmentStates`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
/// ## Parameters:
///
/// * `C`: The [`Ciglet`] type or arbitrary wrapper for generating the ciglets
/// * `D`: If true, ensures adjacent operations are different
/// * `Q`: If true, ensures the number of residues of the query consumed by the
///   alignment is not over [`usize::MAX`]
/// * `R`: If true, ensures the number of residues of the reference consumed by
///   the alignment is not over [`usize::MAX`]
#[derive(Debug)]
pub struct AlignmentStatesArbitrary<C, const D: bool, const Q: bool, const R: bool>(pub AlignmentStates, pub PhantomData<C>);

impl_deref! {AlignmentStatesArbitrary<C, D, Q, R>, AlignmentStates, <C, const D: bool, const Q: bool, const R: bool>}

impl<C, const D: bool, const Q: bool, const R: bool> From<AlignmentStatesArbitrary<C, D, Q, R>> for AlignmentStates {
    #[inline]
    fn from(value: AlignmentStatesArbitrary<C, D, Q, R>) -> Self {
        value.0
    }
}

impl<'a, C, const D: bool, const Q: bool, const R: bool> Arbitrary<'a> for AlignmentStatesArbitrary<C, D, Q, R>
where
    C: Arbitrary<'a> + Into<Ciglet>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let mut vec = Vec::<C>::arbitrary(u)?.into_iter().map(Into::into).collect::<Vec<_>>();

        // Adjust the operations so that adjacent ones are not the same, if D is
        // true
        if let Some(mut last_ciglet) = vec.first().copied()
            && D
        {
            for next_ciglet in &mut vec[1..] {
                if next_ciglet.op == last_ciglet.op {
                    next_ciglet.op = match last_ciglet.op {
                        b'M' => b'I',
                        b'I' => b'D',
                        b'D' => b'N',
                        b'N' => b'S',
                        b'S' => b'H',
                        b'H' => b'P',
                        b'P' => b'=',
                        b'=' => b'X',
                        b'X' => b'M',
                        _ => unreachable!(),
                    }
                }
                last_ciglet = *next_ciglet;
            }
        }

        let mut states = AlignmentStates(vec);

        if Q {
            ensure_no_query_len_overflow(&mut states);
        }

        if R {
            ensure_no_match_len_overflow(&mut states);
        }

        Ok(AlignmentStatesArbitrary(states, PhantomData))
    }
}

/// Removes ciglets from [`AlignmentStates`] to ensure calculating the number of
/// residues it consumes in the query does not overflow.
pub(crate) fn ensure_no_query_len_overflow(states: &mut AlignmentStates) {
    if states.num_query_consumed_checked().is_none() {
        let mut new_vec = Vec::new();
        let mut match_len = 0usize;
        for ciglet in states.iter().copied() {
            if matches!(ciglet.op, b'M' | b'I' | b'S' | b'=' | b'X') {
                match_len = match match_len.checked_add(ciglet.inc) {
                    Some(match_len) => match_len,
                    None => break,
                }
            }
            new_vec.push(ciglet);
        }
        *states = AlignmentStates(new_vec);
    }
}

/// Removes ciglets from [`AlignmentStates`] to ensure calculating the number of
/// residues it consumes in the reference does not overflow.
pub(crate) fn ensure_no_match_len_overflow(states: &mut AlignmentStates) {
    if states.num_ref_consumed_checked().is_none() {
        let mut new_vec = Vec::new();
        let mut match_len = 0usize;
        for ciglet in states.iter().copied() {
            if matches!(ciglet.op, b'M' | b'D' | b'N' | b'=' | b'X') {
                match_len = match match_len.checked_add(ciglet.inc) {
                    Some(match_len) => match_len,
                    None => break,
                }
            }
            new_vec.push(ciglet);
        }
        *states = AlignmentStates(new_vec);
    }
}
