use std::{fmt::Debug, marker::PhantomData};

use super::impl_deref;
use crate::{
    alignment::phmm::{
        Begin, CorePhmm, DomainModule, EmissionParams, End, GlobalPhmm, LayerParams, LocalModule, LocalPhmm, PhmmIndexable,
        PhmmState, SemiLocalModule, TransitionParams,
    },
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

impl<'a, T: Arbitrary<'a>, const S: usize> Arbitrary<'a> for CorePhmm<T, S> {
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let mut layers = vec![LayerParams::<T, S>::arbitrary(u)?, LayerParams::<T, S>::arbitrary(u)?];
        layers.extend(Vec::<LayerParams<T, S>>::arbitrary(u)?);
        Ok(CorePhmm::new_unchecked(layers))
    }
}

impl<'a, T: Arbitrary<'a>> Arbitrary<'a> for SemiLocalModule<T> {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self(Vec::<T>::arbitrary(u)?))
    }
}

impl<'a, T: Arbitrary<'a>, const S: usize> Arbitrary<'a> for DomainModule<T, S> {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self {
            start_to_insert:     T::arbitrary(u)?,
            insert_to_insert:    T::arbitrary(u)?,
            insert_to_end:       T::arbitrary(u)?,
            start_to_end:        T::arbitrary(u)?,
            background_emission: EmissionParams::<T, S>::arbitrary(u)?,
        })
    }
}

impl<'a, T: Arbitrary<'a>, const S: usize> Arbitrary<'a> for LocalModule<T, S> {
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self {
            external_params: SemiLocalModule::arbitrary(u)?,
            internal_params: DomainModule::arbitrary(u)?,
        })
    }
}

/// A wrapper type for [`EmissionParams`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
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

/// A wrapper type for [`TransitionParams`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
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

/// A wrapper type for [`LayerParams`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
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

/// A wrapper type for [`CorePhmm`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
/// Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
/// * `S`: The size of the pHMM alphabet
/// * `R`: If true, ensure invalid transitions are set to infinity
#[derive(Debug)]
pub struct CorePhmmArbitrary<T, F, const S: usize, const R: bool>(pub CorePhmm<T, S>, PhantomData<F>);

impl_deref! {CorePhmmArbitrary<T, F, S, R>, CorePhmm<T, S>, <T: Float, F, const S: usize, const R: bool>}

// Basic implementation: no correction for invalid transitions, no minimum on
// layers
impl<'a, T: Float + Arbitrary<'a>, F, const S: usize> Arbitrary<'a> for CorePhmmArbitrary<T, F, S, false>
where
    LayerParamsArbitrary<T, F, S>: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let mut layers = vec![
            LayerParamsArbitrary::<T, F, S>::arbitrary(u)?.0,
            LayerParamsArbitrary::<T, F, S>::arbitrary(u)?.0,
        ];
        layers.extend(Vec::<LayerParamsArbitrary<T, F, S>>::arbitrary(u)?.into_iter().map(|x| x.0));
        Ok(CorePhmmArbitrary(CorePhmm::new_unchecked(layers), PhantomData))
    }
}

// Implementation that corrects for invalid transitions (calls off to previous
// implementation)
impl<'a, T: Float + Arbitrary<'a>, F, const S: usize> Arbitrary<'a> for CorePhmmArbitrary<T, F, S, true>
where
    CorePhmmArbitrary<T, F, S, false>: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        use PhmmState::*;

        let mut core = CorePhmmArbitrary::<T, F, S, false>::arbitrary(u)?.0;

        let first_layer = core.get_layer_mut(Begin);
        first_layer.transition[(Delete, Delete)] = T::INFINITY;
        first_layer.transition[(Delete, Match)] = T::INFINITY;
        first_layer.transition[(Delete, Insert)] = T::INFINITY;

        let last_layer = core.get_layer_mut(End);
        last_layer.transition[(Delete, Delete)] = T::INFINITY;
        last_layer.transition[(Insert, Delete)] = T::INFINITY;
        last_layer.transition[(Match, Delete)] = T::INFINITY;

        Ok(CorePhmmArbitrary(core, PhantomData))
    }
}

/// A wrapper type for [`SemiLocalModule`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
/// ## Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
pub struct SemiLocalModuleArbitrary<T, F>(SemiLocalModule<T>, PhantomData<F>);

impl_deref! {SemiLocalModuleArbitrary<T, F>, SemiLocalModule<T>, <T, F>}

impl<'a, T, F> Arbitrary<'a> for SemiLocalModuleArbitrary<T, F>
where
    T: Arbitrary<'a>,
    F: Into<T> + Arbitrary<'a>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self(
            SemiLocalModule(
                Vec::<F>::arbitrary(u)?
                    .into_iter()
                    .map(std::convert::Into::into)
                    .collect::<Vec<_>>(),
            ),
            PhantomData,
        ))
    }
}

impl<'a, T, F> SemiLocalModuleArbitrary<T, F>
where
    T: Arbitrary<'a>,
    F: Into<T> + Arbitrary<'a>,
{
    /// Generates an arbitrary [`SemiLocalModule`] which is compatible with
    /// `core` (there are the correct number of pseudomatch states).
    #[allow(clippy::missing_errors_doc)]
    pub fn arbitrary_compatible<const S: usize>(u: &mut Unstructured<'a>, core: &CorePhmm<T, S>) -> Result<Self>
    where
        T: Arbitrary<'a>, {
        Ok(Self(
            SemiLocalModule(
                std::iter::from_fn(|| Some(F::arbitrary(u).map(Into::into)))
                    .take(core.num_pseudomatch())
                    .collect::<Result<Vec<_>, _>>()?,
            ),
            PhantomData,
        ))
    }
}

/// A wrapper type for [`DomainModule`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
/// ## Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
/// * `S`: The size of the pHMM alphabet
pub struct DomainModuleArbitrary<T, F, const S: usize>(DomainModule<T, S>, PhantomData<F>);

impl_deref! {DomainModuleArbitrary<T, F, S>, DomainModule<T, S>, <T, F, const S: usize>}

impl<'a, T, F, const S: usize> Arbitrary<'a> for DomainModuleArbitrary<T, F, S>
where
    T: Arbitrary<'a>,
    F: Into<T> + Arbitrary<'a>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self(
            DomainModule {
                start_to_insert:     F::arbitrary(u)?.into(),
                insert_to_insert:    F::arbitrary(u)?.into(),
                insert_to_end:       F::arbitrary(u)?.into(),
                start_to_end:        F::arbitrary(u)?.into(),
                background_emission: EmissionParamsArbitrary::<T, F, S>::arbitrary(u)?.0,
            },
            PhantomData,
        ))
    }
}

/// A wrapper type for [`LocalModule`] offering more flexibility on the
/// [`Arbitrary`] implementation.
///
/// ## Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
/// * `S`: The size of the pHMM alphabet
pub struct LocalModuleArbitrary<T, F, const S: usize>(LocalModule<T, S>, PhantomData<F>);

impl_deref! {LocalModuleArbitrary<T, F, S>, LocalModule<T, S>, <T, F, const S: usize>}

impl<'a, T, F, const S: usize> Arbitrary<'a> for LocalModuleArbitrary<T, F, S>
where
    T: Arbitrary<'a>,
    F: Into<T> + Arbitrary<'a>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self(
            LocalModule {
                external_params: SemiLocalModuleArbitrary::<T, F>::arbitrary(u)?.0,
                internal_params: DomainModuleArbitrary::<T, F, S>::arbitrary(u)?.0,
            },
            PhantomData,
        ))
    }
}

impl<'a, T, F, const S: usize> LocalModuleArbitrary<T, F, S>
where
    T: Arbitrary<'a>,
    F: Into<T> + Arbitrary<'a>,
{
    /// Generates an arbitrary [`SemiLocalModule`] which is compatible with
    /// `core` (there are the correct number of pseudomatch states).
    #[allow(clippy::missing_errors_doc)]
    pub fn arbitrary_compatible(u: &mut Unstructured<'a>, core: &CorePhmm<T, S>) -> Result<Self>
    where
        T: Arbitrary<'a>, {
        Ok(Self(
            LocalModule {
                external_params: SemiLocalModuleArbitrary::<T, F>::arbitrary_compatible(u, core)?.0,
                internal_params: DomainModuleArbitrary::<T, F, S>::arbitrary(u)?.0,
            },
            PhantomData,
        ))
    }
}

/// A wrapper type for [`GlobalPhmm`] using a DNA alphabet and offering more
/// flexibility on the [`Arbitrary`] implementation.
///
/// Parameters:
///
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
/// * `R`: If true, ensure invalid transitions are set to infinity
#[derive(Debug)]
pub struct DnaGlobalPhmm<T, F, const R: bool>(pub GlobalPhmm<T, 4>, PhantomData<F>);

impl_deref! {DnaGlobalPhmm<T, F, R>, GlobalPhmm<T, 4>, <T: Float, F, const R: bool>}

impl<'a, T: Float + Arbitrary<'a>, F, const R: bool> Arbitrary<'a> for DnaGlobalPhmm<T, F, R>
where
    CorePhmmArbitrary<T, F, 4, R>: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(DnaGlobalPhmm(
            GlobalPhmm {
                mapping: &DNA_UNAMBIG_PROFILE_MAP,
                core:    CorePhmmArbitrary::<T, F, 4, R>::arbitrary(u)?.0,
            },
            PhantomData,
        ))
    }
}

/// A wrapper type for [`LocalPhmm`] using a DNA alphabet and offering more
/// flexibility on the [`Arbitrary`] implementation.
///
/// Parameters:
/// * `T`: The floating point type for the parameters
/// * `F`: The floating point type or arbitrary wrapper for generating the
///   parameters
/// * `M`: If true, ensure there are the [`LocalModule`]s are compatible with
///   the pHMM
/// * `R`: If true, ensure invalid transitions are set to infinity
#[derive(Debug)]
pub struct DnaLocalPhmm<T, F, const M: bool, const R: bool>(pub LocalPhmm<T, 4>, pub PhantomData<F>);

impl_deref! {DnaLocalPhmm<T, F, M, R>, LocalPhmm<T, 4>, <T: Float, F, const M: bool, const R: bool>}

impl<'a, T: Float + Arbitrary<'a>, F, const R: bool> Arbitrary<'a> for DnaLocalPhmm<T, F, false, R>
where
    CorePhmmArbitrary<T, F, 4, R>: Arbitrary<'a>,
    LocalModuleArbitrary<T, F, 4>: Arbitrary<'a>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self(
            LocalPhmm {
                mapping: &DNA_UNAMBIG_PROFILE_MAP,
                core:    CorePhmmArbitrary::<T, F, 4, R>::arbitrary(u)?.0,
                begin:   LocalModuleArbitrary::<T, F, 4>::arbitrary(u)?.0,
                end:     LocalModuleArbitrary::<T, F, 4>::arbitrary(u)?.0,
            },
            PhantomData,
        ))
    }
}

impl<'a, T: Float + Arbitrary<'a>, F, const R: bool> Arbitrary<'a> for DnaLocalPhmm<T, F, true, R>
where
    T: Arbitrary<'a>,
    F: Into<T> + Arbitrary<'a>,
    CorePhmmArbitrary<T, F, 4, R>: Arbitrary<'a>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let core = CorePhmmArbitrary::<T, F, 4, R>::arbitrary(u)?.0;
        let begin = LocalModuleArbitrary::<T, F, 4>::arbitrary_compatible(u, &core)?.0;
        let end = LocalModuleArbitrary::<T, F, 4>::arbitrary_compatible(u, &core)?.0;

        Ok(Self(
            LocalPhmm {
                mapping: &DNA_UNAMBIG_PROFILE_MAP,
                core,
                begin,
                end,
            },
            PhantomData,
        ))
    }
}
