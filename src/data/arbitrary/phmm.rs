use std::{fmt::Debug, marker::PhantomData};

use super::impl_deref;
use crate::{
    alignment::phmm::{EmissionParams, LayerParams, Phmm, PhmmState, TransitionParams},
    data::mappings::DNA_UNAMBIG_PROFILE_MAP,
    math::Float,
};
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a, T: Arbitrary<'a>, const S: usize> Arbitrary<'a> for EmissionParams<T, S> {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(EmissionParams(<[T; S]>::arbitrary(u)?))
    }
}

impl<'a, T: Arbitrary<'a>> Arbitrary<'a> for TransitionParams<T> {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(TransitionParams(<[[T; 3]; 3]>::arbitrary(u)?))
    }
}

impl<'a, T: Arbitrary<'a>, const S: usize> Arbitrary<'a> for LayerParams<T, S> {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(LayerParams {
            transition:      TransitionParams::arbitrary(u)?,
            emission_match:  EmissionParams::arbitrary(u)?,
            emission_insert: EmissionParams::arbitrary(u)?,
        })
    }
}

/// ## Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
/// * `S`: The number of parameters (size of the pHMM alphabet)
pub struct EmissionParamsArbitrary<T, F, const S: usize>(pub EmissionParams<T, S>, PhantomData<F>);

impl_deref! {EmissionParamsArbitrary<T, F, S>, EmissionParams<T, S>, <T, F, const S: usize>}

impl<'a, T, const S: usize, F> Arbitrary<'a> for EmissionParamsArbitrary<T, F, S>
where
    [F; S]: Arbitrary<'a>,
    F: Into<T>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(EmissionParamsArbitrary(
            EmissionParams(<[F; S]>::arbitrary(u)?.map(std::convert::Into::into)),
            PhantomData,
        ))
    }
}

/// ## Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
pub struct TransitionParamsArbitrary<T, F>(pub TransitionParams<T>, PhantomData<F>);

impl_deref! {TransitionParamsArbitrary<T, F>, TransitionParams<T>, <T, F>}

impl<'a, T: Float + Arbitrary<'a>, F> Arbitrary<'a> for TransitionParamsArbitrary<T, F>
where
    [[F; 3]; 3]: Arbitrary<'a>,
    F: Into<T>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(TransitionParamsArbitrary(
            TransitionParams(<[[F; 3]; 3]>::arbitrary(u)?.map(|x| x.map(std::convert::Into::into))),
            PhantomData,
        ))
    }
}

/// ## Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
/// * `S`: The size of the pHMM alphabet
pub struct LayerParamsArbitrary<T, F, const S: usize>(LayerParams<T, S>, PhantomData<F>);

impl_deref! {LayerParamsArbitrary<T, F, S>, LayerParams<T, S>, <T, F, const S: usize>}

impl<'a, T: Arbitrary<'a>, F, const S: usize> Arbitrary<'a> for LayerParamsArbitrary<T, F, S>
where
    TransitionParamsArbitrary<T, F>: Arbitrary<'a>,
    EmissionParamsArbitrary<T, F, S>: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(LayerParamsArbitrary(
            LayerParams {
                transition:      TransitionParamsArbitrary::<T, F>::arbitrary(u)?.0,
                emission_match:  EmissionParamsArbitrary::<T, F, S>::arbitrary(u)?.0,
                emission_insert: EmissionParamsArbitrary::<T, F, S>::arbitrary(u)?.0,
            },
            PhantomData,
        ))
    }
}

/// Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
/// * `R`: If true, ensure invalid transitions are set to infinity
/// * `L`: If true, ensures at least two layers are present
#[derive(Debug)]
pub struct DnaPhmm<T, F, const R: bool, const L: bool>(pub Phmm<T, 4>, PhantomData<F>);

impl_deref! {DnaPhmm<T, F, R, L>, Phmm<T, 4>, <T: Float, F, const R: bool, const L: bool>}

// Basic implementation: no correction for invalid transitions, no minimum on
// layers
impl<'a, T: Float + Arbitrary<'a>, F> Arbitrary<'a> for DnaPhmm<T, F, false, false>
where
    LayerParamsArbitrary<T, F, 4>: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(DnaPhmm(
            Phmm {
                mapping: &DNA_UNAMBIG_PROFILE_MAP,
                params:  Vec::<LayerParamsArbitrary<T, F, 4>>::arbitrary(u)?
                    .into_iter()
                    .map(|x| x.0)
                    .collect::<Vec<_>>(),
            },
            PhantomData,
        ))
    }
}

// Implementation that adds minimum on number of layers (calls off to basic
// implementation)
impl<'a, T: Float + Arbitrary<'a>, F> Arbitrary<'a> for DnaPhmm<T, F, false, true>
where
    LayerParamsArbitrary<T, F, 4>: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let first_layer = LayerParamsArbitrary::<T, F, 4>::arbitrary(u)?;
        let last_layer = LayerParamsArbitrary::<T, F, 4>::arbitrary(u)?;
        let mut phmm = DnaPhmm::<T, F, false, false>::arbitrary(u)?;
        phmm.params.insert(0, first_layer.0);
        phmm.params.push(last_layer.0);
        Ok(DnaPhmm(phmm.0, PhantomData))
    }
}

// Implementation that corrects for invalid transitions (calls off to previous
// implementation)
impl<'a, T: Float + Arbitrary<'a>, F, const L: bool> Arbitrary<'a> for DnaPhmm<T, F, true, L>
where
    DnaPhmm<T, F, false, L>: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        use PhmmState::*;

        let mut out = DnaPhmm::<T, F, false, L>::arbitrary(u)?;

        if let Some(first_layer) = out.0.params.first_mut() {
            first_layer.transition[(Delete, Delete)] = T::INFINITY;
            first_layer.transition[(Delete, Match)] = T::INFINITY;
            first_layer.transition[(Delete, Insert)] = T::INFINITY;
        }

        if let Some(last_layer) = out.0.params.last_mut() {
            last_layer.transition[(Delete, Delete)] = T::INFINITY;
            last_layer.transition[(Insert, Delete)] = T::INFINITY;
            last_layer.transition[(Match, Delete)] = T::INFINITY;
        }

        Ok(DnaPhmm(out.0, PhantomData))
    }
}
