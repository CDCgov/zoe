//! Arbitrary implementations and specification structs for profile Hidden
//! Markov Models.

use crate::{
    alignment::phmm::{
        CorePhmm, DomainPhmm, EmissionParams, GetLayer, GlobalPhmm, LayerParams, LocalPhmm, PhmmState, SemiLocalPhmm,
        TransitionParams,
        indexing::{Begin, End, PhmmIndexable},
        modules::{DomainModule, LocalModule, SemiLocalModule},
    },
    data::{
        arbitrary::{ArbitrarySpecs, ArraySpecs, FloatSpecs, VecSpecs},
        mappings::DNA_UNAMBIG_PROFILE_MAP,
    },
    math::Float,
};
use arbitrary::{Arbitrary, Result, Unstructured};

impl<'a, T, const S: usize> Arbitrary<'a> for EmissionParams<T, S>
where
    T: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(EmissionParams::from_array(<[T; S]>::arbitrary(u)?))
    }
}

impl<'a, T> Arbitrary<'a> for TransitionParams<T>
where
    T: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(TransitionParams(<[[T; 3]; 3]>::arbitrary(u)?))
    }
}

impl<'a, T, const S: usize> Arbitrary<'a> for LayerParams<T, S>
where
    T: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(LayerParams {
            transition:      TransitionParams::arbitrary(u)?,
            emission_match:  EmissionParams::arbitrary(u)?,
            emission_insert: EmissionParams::arbitrary(u)?,
        })
    }
}

impl<'a, T, const S: usize> Arbitrary<'a> for CorePhmm<T, S>
where
    T: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let mut layers = vec![LayerParams::<T, S>::arbitrary(u)?, LayerParams::<T, S>::arbitrary(u)?];
        layers.extend(Vec::<LayerParams<T, S>>::arbitrary(u)?);
        Ok(CorePhmm::new_unchecked(layers))
    }
}

impl<'a, T> Arbitrary<'a> for SemiLocalModule<T>
where
    T: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self(Vec::<T>::arbitrary(u)?))
    }
}

impl<'a, T, const S: usize> Arbitrary<'a> for DomainModule<T, S>
where
    T: Arbitrary<'a>,
{
    #[inline]
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

impl<'a, T, const S: usize> Arbitrary<'a> for LocalModule<T, S>
where
    T: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self {
            external_params: SemiLocalModule::arbitrary(u)?,
            internal_params: DomainModule::arbitrary(u)?,
        })
    }
}

/// Specifications for generating an arbitrary [`EmissionParams`] struct.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct EmissionParamsSpecs<T, const S: usize> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,
}

impl<'a, T, const S: usize> ArbitrarySpecs<'a> for EmissionParamsSpecs<T, S>
where
    T: Float + Arbitrary<'a>,
{
    type Output = EmissionParams<T, S>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let specs = ArraySpecs {
            element_specs: self.float_specs,
        };

        Ok(EmissionParams::from_array(specs.make_arbitrary(u)?))
    }
}

/// Specifications for generating an arbitrary [`TransitionParamsSpecs`] struct.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
struct TransitionParamsSpecs<T> {
    /// The specifications for generating the floating point parameters.
    float_specs: FloatSpecs<T>,
}

impl<'a, T> ArbitrarySpecs<'a> for TransitionParamsSpecs<T>
where
    T: Float + Arbitrary<'a>,
{
    type Output = TransitionParams<T>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let specs = ArraySpecs {
            element_specs: ArraySpecs {
                element_specs: self.float_specs,
            },
        };

        Ok(TransitionParams(specs.make_arbitrary(u)?))
    }
}

/// Specifications for generating an arbitrary [`LayerParams`] struct.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
struct LayerParamsSpecs<T, const S: usize> {
    /// The specifications for generating the floating point parameters.
    float_specs: FloatSpecs<T>,
}

impl<'a, T, const S: usize> ArbitrarySpecs<'a> for LayerParamsSpecs<T, S>
where
    T: Float + Arbitrary<'a>,
{
    type Output = LayerParams<T, S>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let transition_specs = TransitionParamsSpecs {
            float_specs: self.float_specs,
        };

        let emission_specs = EmissionParamsSpecs {
            float_specs: self.float_specs,
        };

        Ok(LayerParams {
            transition:      transition_specs.make_arbitrary(u)?,
            emission_match:  emission_specs.make_arbitrary(u)?,
            emission_insert: emission_specs.make_arbitrary(u)?,
        })
    }
}

/// Specifications for generating an arbitrary [`CorePhmm`].
#[cfg_attr(feature = "alignment-diagnostics", visibility::make(pub))]
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub(crate) struct CorePhmmSpecs<T, const S: usize> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,

    /// Whether to disallow invalid transitions for the first and last layers by
    /// setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,
}

impl<'a, T, const S: usize> ArbitrarySpecs<'a> for CorePhmmSpecs<T, S>
where
    T: Float + Arbitrary<'a>,
{
    type Output = CorePhmm<T, S>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        use PhmmState::*;

        let specs = VecSpecs {
            element_specs: LayerParamsSpecs {
                float_specs: self.float_specs,
            },
            min_len:       2,
            len:           None,
            max_len:       usize::MAX,
        };

        let mut core = CorePhmm::new_unchecked(specs.make_arbitrary(u)?);

        if self.disallow_invalid {
            let first_layer = core.get_layer_mut(Begin);
            first_layer.transition[(Delete, Delete)] = T::INFINITY;
            first_layer.transition[(Delete, Match)] = T::INFINITY;
            first_layer.transition[(Delete, Insert)] = T::INFINITY;

            let last_layer = core.get_layer_mut(End);
            last_layer.transition[(Delete, Delete)] = T::INFINITY;
            last_layer.transition[(Insert, Delete)] = T::INFINITY;
            last_layer.transition[(Match, Delete)] = T::INFINITY;
        }

        Ok(core)
    }
}

/// Specifications for generating an arbitrary [`SemiLocalModule`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct SemiLocalModuleSpecs<T> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,

    /// The number of match states (including BEGIN and END) to include.
    pub num_pseudomatch: Option<usize>,
}

impl<'a, T> ArbitrarySpecs<'a> for SemiLocalModuleSpecs<T>
where
    T: Float + Arbitrary<'a>,
{
    type Output = SemiLocalModule<T>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let specs = VecSpecs {
            element_specs: self.float_specs,
            min_len:       0,
            len:           self.num_pseudomatch,
            max_len:       usize::MAX,
        };

        specs.make_arbitrary(u).map(SemiLocalModule)
    }
}

/// Specifications for generating an arbitrary [`DomainModule`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DomainModuleSpecs<T, const S: usize> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,
}

impl<'a, T, const S: usize> ArbitrarySpecs<'a> for DomainModuleSpecs<T, S>
where
    T: Float + Arbitrary<'a>,
{
    type Output = DomainModule<T, S>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let emission_specs = EmissionParamsSpecs {
            float_specs: self.float_specs,
        };

        Ok(DomainModule {
            start_to_insert:     self.float_specs.make_arbitrary(u)?,
            insert_to_insert:    self.float_specs.make_arbitrary(u)?,
            insert_to_end:       self.float_specs.make_arbitrary(u)?,
            start_to_end:        self.float_specs.make_arbitrary(u)?,
            background_emission: emission_specs.make_arbitrary(u)?,
        })
    }
}

/// Specification for generating an arbitrary [`LocalModule`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct LocalModuleSpecs<T, const S: usize> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,

    /// The number of match states (including BEGIN and END) to include.
    pub num_pseudomatch: Option<usize>,
}

impl<'a, T, const S: usize> ArbitrarySpecs<'a> for LocalModuleSpecs<T, S>
where
    T: Float + Arbitrary<'a>,
{
    type Output = LocalModule<T, S>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let semilocal_specs = SemiLocalModuleSpecs {
            float_specs:     self.float_specs,
            num_pseudomatch: self.num_pseudomatch,
        };

        let domain_specs = DomainModuleSpecs {
            float_specs: self.float_specs,
        };

        Ok(LocalModule {
            external_params: semilocal_specs.make_arbitrary(u)?,
            internal_params: domain_specs.make_arbitrary(u)?,
        })
    }
}

/// Specifications for generating an arbitrary [`GlobalPhmm`] with a DNA
/// alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DnaGlobalPhmmSpecs<T> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,

    /// Whether to disallow invalid transitions for the first and last layers by
    /// setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,
}

impl<'a, T> ArbitrarySpecs<'a> for DnaGlobalPhmmSpecs<T>
where
    T: Float + Arbitrary<'a>,
{
    type Output = GlobalPhmm<T, 4>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let specs = CorePhmmSpecs {
            float_specs:      self.float_specs,
            disallow_invalid: self.disallow_invalid,
        };

        Ok(GlobalPhmm::new(&DNA_UNAMBIG_PROFILE_MAP, specs.make_arbitrary(u)?))
    }
}

/// Specifications for generating an arbitrary [`LocalPhmm`] with a DNA
/// alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DnaLocalPhmmSpecs<T> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,

    /// Whether to disallow invalid transitions for the first and last layers by
    /// setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,

    /// Whether to ensure that the modules have a compatible size with the core
    /// pHMM.
    pub compatible_modules: bool,
}

impl<'a, T> ArbitrarySpecs<'a> for DnaLocalPhmmSpecs<T>
where
    T: Float + Arbitrary<'a>,
{
    type Output = LocalPhmm<T, 4>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let core_specs = CorePhmmSpecs {
            float_specs:      self.float_specs,
            disallow_invalid: self.disallow_invalid,
        };

        let core = core_specs.make_arbitrary(u)?;

        let module_specs = LocalModuleSpecs {
            float_specs:     self.float_specs,
            num_pseudomatch: self.compatible_modules.then_some(core.num_pseudomatch()),
        };

        let begin = module_specs.make_arbitrary(u)?;
        let end = module_specs.make_arbitrary(u)?;

        Ok(LocalPhmm::new(&DNA_UNAMBIG_PROFILE_MAP, core, begin, end))
    }
}

/// Specifications for generating an arbitrary [`DomainPhmm`] with a DNA
/// alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DnaDomainPhmmSpecs<T> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,

    /// Whether to disallow invalid transitions for the first and last layers by
    /// setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,
}

impl<'a, T> ArbitrarySpecs<'a> for DnaDomainPhmmSpecs<T>
where
    T: Float + Arbitrary<'a>,
{
    type Output = DomainPhmm<T, 4>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let core_specs = CorePhmmSpecs {
            float_specs:      self.float_specs,
            disallow_invalid: self.disallow_invalid,
        };

        let core = core_specs.make_arbitrary(u)?;

        let module_specs = DomainModuleSpecs {
            float_specs: self.float_specs,
        };

        let begin = module_specs.make_arbitrary(u)?;
        let end = module_specs.make_arbitrary(u)?;

        Ok(DomainPhmm::new(&DNA_UNAMBIG_PROFILE_MAP, core, begin, end))
    }
}

/// Specifications for generating an arbitrary [`SemiLocalPhmm`] with a DNA
/// alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DnaSemiLocalPhmmSpecs<T> {
    /// The specifications for generating the floating point parameters.
    pub float_specs: FloatSpecs<T>,

    /// Whether to disallow invalid transitions for the first and last layers by
    /// setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,

    /// Whether to ensure that the modules have a compatible size with the core
    /// pHMM.
    pub compatible_modules: bool,
}

impl<'a, T> ArbitrarySpecs<'a> for DnaSemiLocalPhmmSpecs<T>
where
    T: Float + Arbitrary<'a>,
{
    type Output = SemiLocalPhmm<T, 4>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let core_specs = CorePhmmSpecs {
            float_specs:      self.float_specs,
            disallow_invalid: self.disallow_invalid,
        };

        let core = core_specs.make_arbitrary(u)?;

        let module_specs = SemiLocalModuleSpecs {
            float_specs:     self.float_specs,
            num_pseudomatch: self.compatible_modules.then_some(core.num_pseudomatch()),
        };

        let begin = module_specs.make_arbitrary(u)?;
        let end = module_specs.make_arbitrary(u)?;

        Ok(SemiLocalPhmm::new(&DNA_UNAMBIG_PROFILE_MAP, core, begin, end))
    }
}
