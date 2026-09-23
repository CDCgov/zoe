//! Arbitrary implementations and specification structs for profile Hidden
//! Markov Models.

use crate::{
    alignment::phmm::{
        DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
        at_least_two::VecAtLeast2,
        components::EmissionParams,
        indexing::PhmmLen,
        modules::{DomainModule, LocalModule, SemiLocalModule},
    },
    data::{
        arbitrary::{
            ArbitrarySpecs, NoConstraintSpecs, VecSpecs,
            components::{CorePhmmSpecs, EmissionParamsSpecs},
        },
        mappings::DNA_UNAMBIG_PROFILE_MAP,
    },
    iter_utils::ProcessResultsExt,
};
use arbitrary::{Arbitrary, Result, Unstructured};

pub mod components;

impl<'a, T> Arbitrary<'a> for VecAtLeast2<T>
where
    T: Arbitrary<'a>,
{
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        let specs = VecSpecs {
            element_specs: NoConstraintSpecs::default(),
            min_len: 2,
            ..Default::default()
        };

        let vec = specs.make_arbitrary(u)?;
        Ok(vec.try_into().unwrap())
    }
}

/// Specifications for generating an arbitrary [`VecAtLeast2`].
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct VecAtLeast2Specs<S> {
    /// The specifications for generating each element of the [`VecAtLeast2`].
    pub element_specs: S,

    /// The minimum length of the [`VecAtLeast2`].
    ///
    /// This must be less than or equal to `max_len`, and will be clamped to be
    /// above 2 during usage.
    pub min_len: usize,

    /// The exact length of the [`VecAtLeast2`] to generate.
    ///
    /// If set, this ignores the `min_len` and `max_len` fields. This must be at
    /// least 2.
    pub len: Option<usize>,

    /// The maximum length of the [`VecAtLeast2`].
    ///
    /// This must be greater than or equal to `min_len`. This must be at least
    /// 2.
    pub max_len: usize,
}

impl<S> Default for VecAtLeast2Specs<S>
where
    S: Default,
{
    fn default() -> Self {
        Self {
            element_specs: S::default(),
            min_len:       2,
            len:           None,
            max_len:       usize::MAX,
        }
    }
}

impl<'a, S> ArbitrarySpecs<'a> for VecAtLeast2Specs<S>
where
    S: ArbitrarySpecs<'a>,
{
    type Output = VecAtLeast2<S::Output>;

    /// Generates an arbitrary [`VecAtLeast2`] conforming to the given
    /// specifications.
    ///
    /// ## Errors
    ///
    /// Any errors from the underlying [`arbitrary`] calls are propagated.
    ///
    /// ## Panics
    ///
    /// `min_len` must be less than or equal to `max_len`, and `len` and
    /// `max_len` cannot be less than 2.
    ///
    /// [`arbitrary`]: arbitrary::Arbitrary::arbitrary
    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        assert!(
            self.max_len >= 2,
            "The max_len field must be at least 2 for VecAtLeast2Specs (found {})",
            self.max_len
        );

        let min_len = self.min_len.max(2);

        let Some(len_range) = self.max_len.checked_sub(min_len) else {
            panic!(
                "The min_len field must be less than or equal to the max_len field for VecAtLeast2Specs (found min_len={min_len} and max_len={max_len}",
                max_len = self.max_len
            );
        };

        let vec = if let Some(len) = self.len {
            assert!(
                len >= 2,
                "The len field must be at least 2 for VecAtLeast2Specs (found {len})"
            );

            std::iter::repeat_with(|| self.element_specs.make_arbitrary(u))
                .take(len)
                .collect::<Result<Vec<_>>>()?
        } else {
            let start = std::iter::repeat_with(|| self.element_specs.make_arbitrary(u)).take(min_len);

            let mut out = start.collect::<Result<Vec<_>>>()?;

            let remaining = self.element_specs.make_arbitrary_iter(u).take(len_range);

            remaining.process_results(|iter| {
                out.extend(iter);
            })?;

            out
        };

        Ok(vec.try_into().unwrap())
    }
}

impl<'a, T> Arbitrary<'a> for SemiLocalModule<T>
where
    T: Arbitrary<'a>,
{
    #[inline]
    fn arbitrary(u: &mut Unstructured<'a>) -> Result<Self> {
        Ok(Self(VecAtLeast2::arbitrary(u)?))
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
        let specs = VecAtLeast2Specs {
            element_specs: self.param_specs,
            min_len:       2,
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
    type Output = GlobalPhmm<'static, K::Output, 4>;

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
    type Output = LocalPhmm<'static, K::Output, 4>;

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
    type Output = DomainPhmm<'static, K::Output, 4>;

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
    type Output = SemiLocalPhmm<'static, K::Output, 4>;

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
