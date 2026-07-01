//! Arbitrary implementations and specification structs for profile Hidden
//! Markov Models.

use crate::{
    alignment::phmm::{
        DomainPhmm, EmissionParams, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
        indexing::PhmmIndexable,
        modules::{DomainModule, LocalModule, SemiLocalModule},
    },
    data::{
        arbitrary::{ArbitrarySpecs, VecSpecs},
        mappings::DNA_UNAMBIG_PROFILE_MAP,
    },
};
use arbitrary::{Arbitrary, Result, Unstructured};

#[cfg(not(feature = "alignment-diagnostics"))]
#[doc(auto_cfg(hide(feature, values("alignment-diagnostics"))))]
mod components;
#[cfg(not(feature = "alignment-diagnostics"))]
#[doc(auto_cfg(hide(feature, values("alignment-diagnostics"))))]
pub use components::EmissionParamsSpecs;
#[cfg(not(feature = "alignment-diagnostics"))]
#[doc(auto_cfg(hide(feature, values("alignment-diagnostics"))))]
pub(crate) use components::*;

#[cfg(feature = "alignment-diagnostics")]
pub(crate) mod components;
#[cfg(feature = "alignment-diagnostics")]
pub use components::*;

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
            semilocal_params: SemiLocalModule::arbitrary(u)?,
            domain_params:    DomainModule::arbitrary(u)?,
        })
    }
}

/// Specifications for generating an arbitrary [`SemiLocalModule`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct SemiLocalModuleSpecs<K> {
    /// The specifications for generating the parameters.
    pub param_specs: K,

    /// The number of match states (including BEGIN and END) to include.
    pub num_pseudomatch: Option<usize>,
}

impl<'a, K> ArbitrarySpecs<'a> for SemiLocalModuleSpecs<K>
where
    K: ArbitrarySpecs<'a, Output: Default> + Copy,
{
    type Output = SemiLocalModule<K::Output>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let specs = VecSpecs {
            element_specs: self.param_specs,
            min_len:       0,
            len:           self.num_pseudomatch,
            max_len:       usize::MAX,
        };

        specs.make_arbitrary(u).map(SemiLocalModule)
    }
}

/// Specifications for generating an arbitrary [`DomainModule`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DomainModuleSpecs<K, const S: usize> {
    /// The specifications for generating the parameters.
    pub param_specs: K,
}

impl<'a, K, const S: usize> ArbitrarySpecs<'a> for DomainModuleSpecs<K, S>
where
    K: ArbitrarySpecs<'a, Output: Default> + Copy,
{
    type Output = DomainModule<K::Output, S>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let emission_specs = EmissionParamsSpecs {
            param_specs: self.param_specs,
        };

        Ok(DomainModule {
            start_to_insert:     self.param_specs.make_arbitrary(u)?,
            insert_to_insert:    self.param_specs.make_arbitrary(u)?,
            insert_to_end:       self.param_specs.make_arbitrary(u)?,
            start_to_end:        self.param_specs.make_arbitrary(u)?,
            background_emission: emission_specs.make_arbitrary(u)?,
        })
    }
}

/// Specification for generating an arbitrary [`LocalModule`].
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct LocalModuleSpecs<K, const S: usize> {
    /// The specifications for generating the floating point parameters.
    pub param_specs: K,

    /// The number of match states (including BEGIN and END) to include.
    pub num_pseudomatch: Option<usize>,
}

impl<'a, K, const S: usize> ArbitrarySpecs<'a> for LocalModuleSpecs<K, S>
where
    K: ArbitrarySpecs<'a, Output: Default> + Copy,
{
    type Output = LocalModule<K::Output, S>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let semilocal_specs = SemiLocalModuleSpecs {
            param_specs:     self.param_specs,
            num_pseudomatch: self.num_pseudomatch,
        };

        let domain_specs = DomainModuleSpecs {
            param_specs: self.param_specs,
        };

        Ok(LocalModule {
            semilocal_params: semilocal_specs.make_arbitrary(u)?,
            domain_params:    domain_specs.make_arbitrary(u)?,
        })
    }
}

/// Specifications for generating an arbitrary [`GlobalPhmm`] with a DNA
/// alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DnaGlobalPhmmSpecs<K> {
    /// The specifications for generating the parameters.
    pub param_specs: K,

    /// Whether to disallow invalid transitions/emissions for the first and last
    /// layers by setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,
}

impl<'a, K> ArbitrarySpecs<'a> for DnaGlobalPhmmSpecs<K>
where
    K: ArbitrarySpecs<'a, Output: PhmmNumber> + Copy,
{
    type Output = GlobalPhmm<K::Output, 4>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let specs = CorePhmmSpecs {
            param_specs:      self.param_specs,
            disallow_invalid: self.disallow_invalid,
        };

        Ok(GlobalPhmm {
            mapping: &DNA_UNAMBIG_PROFILE_MAP,
            core:    specs.make_arbitrary(u)?,
        })
    }
}

/// Specifications for generating an arbitrary [`LocalPhmm`] with a DNA
/// alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DnaLocalPhmmSpecs<K> {
    /// The specifications for generating the parameters.
    pub param_specs: K,

    /// Whether to disallow invalid transitions/emissions for the first and last
    /// layers by setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,

    /// Whether to ensure that the modules have a compatible size with the core
    /// pHMM.
    pub compatible_modules: bool,
}

impl<'a, K> ArbitrarySpecs<'a> for DnaLocalPhmmSpecs<K>
where
    K: ArbitrarySpecs<'a, Output: PhmmNumber> + Copy,
{
    type Output = LocalPhmm<K::Output, 4>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let core_specs = CorePhmmSpecs {
            param_specs:      self.param_specs,
            disallow_invalid: self.disallow_invalid,
        };

        let core = core_specs.make_arbitrary(u)?;

        let module_specs = LocalModuleSpecs {
            param_specs:     self.param_specs,
            num_pseudomatch: self.compatible_modules.then_some(core.num_pseudomatch()),
        };

        let begin = module_specs.make_arbitrary(u)?;
        let end = module_specs.make_arbitrary(u)?;

        Ok(LocalPhmm {
            mapping: &DNA_UNAMBIG_PROFILE_MAP,
            core,
            begin,
            end,
        })
    }
}

/// Specifications for generating an arbitrary [`DomainPhmm`] with a DNA
/// alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DnaDomainPhmmSpecs<K> {
    /// The specifications for generating the floating point parameters.
    pub param_specs: K,

    /// Whether to disallow invalid transitions/emissions for the first and last
    /// layers by setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,
}

impl<'a, K> ArbitrarySpecs<'a> for DnaDomainPhmmSpecs<K>
where
    K: ArbitrarySpecs<'a, Output: PhmmNumber> + Copy,
{
    type Output = DomainPhmm<K::Output, 4>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let core_specs = CorePhmmSpecs {
            param_specs:      self.param_specs,
            disallow_invalid: self.disallow_invalid,
        };

        let core = core_specs.make_arbitrary(u)?;

        let module_specs = DomainModuleSpecs {
            param_specs: self.param_specs,
        };

        let begin = module_specs.make_arbitrary(u)?;
        let end = module_specs.make_arbitrary(u)?;

        Ok(DomainPhmm {
            mapping: &DNA_UNAMBIG_PROFILE_MAP,
            core,
            begin,
            end,
        })
    }
}

/// Specifications for generating an arbitrary [`SemiLocalPhmm`] with a DNA
/// alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct DnaSemiLocalPhmmSpecs<K> {
    /// The specifications for generating the parameters.
    pub param_specs: K,

    /// Whether to disallow invalid transitions/emissions for the first and last
    /// layers by setting the parameters to infinity (probability zero).
    pub disallow_invalid: bool,

    /// Whether to ensure that the modules have a compatible size with the core
    /// pHMM.
    pub compatible_modules: bool,
}

impl<'a, K> ArbitrarySpecs<'a> for DnaSemiLocalPhmmSpecs<K>
where
    K: ArbitrarySpecs<'a, Output: PhmmNumber> + Copy,
{
    type Output = SemiLocalPhmm<K::Output, 4>;

    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        let core_specs = CorePhmmSpecs {
            param_specs:      self.param_specs,
            disallow_invalid: self.disallow_invalid,
        };

        let core = core_specs.make_arbitrary(u)?;

        let module_specs = SemiLocalModuleSpecs {
            param_specs:     self.param_specs,
            num_pseudomatch: self.compatible_modules.then_some(core.num_pseudomatch()),
        };

        let begin = module_specs.make_arbitrary(u)?;
        let end = module_specs.make_arbitrary(u)?;

        Ok(SemiLocalPhmm {
            mapping: &DNA_UNAMBIG_PROFILE_MAP,
            core,
            begin,
            end,
        })
    }
}
