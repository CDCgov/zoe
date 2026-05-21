//! Structs defining the modules that appear may appear at the beginning and/or
//! end of a pHMM in order to change the alignment mode.
//!
//! Traditionally, a pHMM performs global alignment by default (as in
//! [`GlobalPhmm`]). However, by adding extra states/transitions to the
//! beginning/end of a pHMM, it can allow other types of alignment by
//! permitting:
//!
//! - Extra residues to be consumed at the tail ends of a pHMM with differing
//!   scoring criteria
//! - The alignment to skip to an arbitrary layer in the pHMM, rather than
//!   passing through the BEGIN state
//! - The alignment to exit early from a layer in the pHMM, rather than continue
//!   to the END state
//!
//! [`GlobalPhmm`]: super::GlobalPhmm

use crate::{
    alignment::phmm::{
        CorePhmm, EmissionParams, PhmmNumber,
        indexing::{PhmmIndex, PhmmIndexable, QueryIndex, QueryIndexable},
    },
    data::mappings::ByteIndexMap,
};

/// A module placed before or after a pHMM to support semilocal alignment
/// (aligning the full sequence to a portion of the model).
///
/// This module can jump to/from an arbitrary match statement in the pHMM. Each
/// match state (and the BEGIN state and END state) have a transition parameter
/// associated with them.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SemiLocalModule<T>(pub(crate) Vec<T>);

impl<T> SemiLocalModule<T> {
    /// Returns the parameters as a slice.
    #[inline]
    #[must_use]
    #[cfg(feature = "alignment-diagnostics")]
    pub fn as_slice(&self) -> &[T] {
        &self.0
    }
}

impl<T: PhmmNumber> SemiLocalModule<T> {
    /// Constructs a [`SemiLocalModule`] from a slice of parameters.
    #[cfg(feature = "alignment-diagnostics")]
    pub fn from_slice(params: &[T]) -> Self {
        Self(params.to_vec())
    }

    /// Constructs a [`SemiLocalModule`] where all transitions are free.
    ///
    /// This will not produce a true probabilistic model where the probabilities
    /// sum to 1. However, it can be used for simplicity and ease of
    /// interpretation.
    #[inline]
    #[must_use]
    pub(crate) fn no_penalty<const S: usize>(core: &CorePhmm<T, S>) -> Self {
        Self(vec![T::ZERO; core.num_pseudomatch()])
    }

    /// Gets the score for entering directly into layer `j` (when this module is
    /// placed at the beginning of the [`CorePhmm`] or exiting early from layer
    /// `j` (when this module is placed at the end of the [`CorePhmm`]).
    #[inline]
    #[must_use]
    pub(crate) fn get_score(&self, j: impl PhmmIndex) -> T {
        self.0[self.get_dp_index(j)]
    }
}

/// A module paced before or after a pHMM to support domain alignment (aligning
/// a subsequence to the full model).
///
/// This module can skip arbitrarily many residues from the beginning/end of the
/// query sequence.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct DomainModule<T, const S: usize> {
    /// The transition parameter from the start of the module to the insert
    /// state
    pub start_to_insert:     T,
    /// The transition parameter for staying in the insert state
    pub insert_to_insert:    T,
    /// The transition parameter for exiting the insert state and going to the
    /// end of the module
    pub insert_to_end:       T,
    /// The transition parameter for going from the start to the end of the
    /// module without any insertions (no residues in the sequence are skipped)
    pub start_to_end:        T,
    /// The emission parameters for all inserted/skipped residues
    pub background_emission: EmissionParams<T, S>,
}

impl<T: PhmmNumber, const S: usize> DomainModule<T, S> {
    /// Precomputes the transition probabilities for skipping any number of
    /// residues from the beginning of `seq`.
    ///
    /// This should only be called on modules placed at the beginning of the
    /// pHMM. See [`PrecomputedDomainModule`] for more details.
    pub(crate) fn precompute_begin_mod<Q: AsRef<[u8]>>(
        &self, seq: Q, mapping: &ByteIndexMap<S>,
    ) -> PrecomputedDomainModule<T, S> {
        let seq = seq.as_ref();
        let mut domain_params = vec![self.start_to_insert + self.insert_to_end; seq.len() + 1];
        domain_params[0] = self.start_to_end;

        let mut emissions_sum = T::ZERO;
        for (i, x_idx) in seq.iter().map(|x| mapping.to_index(*x)).enumerate() {
            emissions_sum += self.background_emission[x_idx];
            // Special casing needed in case insert_to_insert is infinite,
            // causing a NAN to appear when multiplied by 0
            let insert_to_insert = if i > 0 {
                T::cast_from(i) * self.insert_to_insert
            } else {
                T::ZERO
            };
            domain_params[i + 1] += emissions_sum + insert_to_insert;
        }

        PrecomputedDomainModule(domain_params)
    }

    /// Precomputes the transition probabilities for skipping any number of
    /// residues from the end of `seq`.
    ///
    /// This should only be called on modules placed at the end of the pHMM. See
    /// [`PrecomputedDomainModule`] for more details.
    pub(crate) fn precompute_end_mod<Q: AsRef<[u8]>>(
        &self, seq: Q, mapping: &ByteIndexMap<S>,
    ) -> PrecomputedDomainModule<T, S> {
        let seq = seq.as_ref();
        let mut domain_params = vec![self.start_to_insert + self.insert_to_end; seq.len() + 1];
        domain_params[seq.len()] = self.start_to_end;

        let mut emissions_sum = T::ZERO;
        for (i, x_idx) in seq.iter().rev().map(|x| mapping.to_index(*x)).enumerate() {
            emissions_sum += self.background_emission[x_idx];
            // Special casing needed in case insert_to_insert is infinite,
            // causing a NAN to appear when multiplied by 0
            let insert_to_insert = if i > 0 {
                T::cast_from(i) * self.insert_to_insert
            } else {
                T::ZERO
            };
            domain_params[seq.len() - i - 1] += emissions_sum + insert_to_insert;
        }

        PrecomputedDomainModule(domain_params)
    }
}

/// Precomputed parameters from a [`DomainModule`] for use in a Viterbi
/// alignment.
///
/// When a [`CorePhmm`] includes a [`DomainModule`] at the beginning, the
/// transition probabilities of skipping $i$ residues from the start of the
/// query require $O(i)$ time. Over all $i$, this results in $O(n^2)$ time,
/// where $n$ is the sequence length.
///
/// To avoid this, [`PrecomputedDomainModule`] carries out all the computations
/// beforehand in $O(n)$ time.
pub(crate) struct PrecomputedDomainModule<T, const S: usize>(pub(crate) Vec<T>);

impl<T: Copy, const S: usize> PrecomputedDomainModule<T, S> {
    /// Gets the score for skipping the first `i` residues in the query (when
    /// this module is placed at the beginning of the [`CorePhmm`]) or skipping
    /// the last `i` residues in the query (when this module is placed at the
    /// end of the [`CorePhmm`]).
    pub(crate) fn get_score(&self, i: impl QueryIndex) -> T {
        self.0[self.get_dp_index(i)]
    }
}

/// A module placed before or after a pHMM to support local alignment (aligning
/// a subsequence to a portion of the model).
///
/// This module can skip arbitrarily many residues from the beginning/end of the
/// query sequence, and it can also jump to/from an arbitrary match statement in
/// the pHMM. Internally, this is a combination of a [`SemiLocalModule`] (for
/// supporting arbitrary jumps to/from the pHMM) and a [`DomainModule`] (for
/// skipping residues in the query).
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct LocalModule<T, const S: usize> {
    /// The parameters enabling semilocal behavior (skipping layers from the
    /// beginning or end of the model).
    ///
    /// These are either the transition parameters into the match states of the
    /// pHMM (if the module is at the beginning) or the transition parameters
    /// out of the match states of the pHMM (if the module is at the end).
    pub semilocal_params: SemiLocalModule<T>,
    /// The parameters enabling domain behavior (skipping residues at the start
    /// and end of the query).
    pub domain_params:    DomainModule<T, S>,
}

impl<T: PhmmNumber, const S: usize> LocalModule<T, S> {
    /// Constructs a [`LocalModule`] where all transitions are free, and the
    /// only penalty is the emission parameters for any skipped bases in the
    /// query.
    ///
    /// This will not produce a true probabilistic model where the probabilities
    /// sum to 1. However, it can be used for simplicity and ease of
    /// interpretation.
    #[inline]
    #[must_use]
    pub(crate) fn no_penalty(core: &CorePhmm<T, S>, background_emission: EmissionParams<T, S>) -> Self {
        Self {
            semilocal_params: SemiLocalModule::no_penalty(core),
            domain_params:    DomainModule::no_penalty(background_emission),
        }
    }
}

/// Precomputed parameters from a [`LocalModule`] for use in a Viterbi
/// alignment.
///
/// When a [`CorePhmm`] includes a [`LocalModule`] at the beginning, every match
/// state computation is influenced by the [`LocalModule`], due to the
/// possibility of skipping directly to that layer. Since this would result in
/// repeated computations, the [`PrecomputedLocalModule`] carries out all such
/// computations beforehand. These are the same computations performed by
/// [`PrecomputedDomainModule`].
pub(crate) struct PrecomputedLocalModule<'a, T, const S: usize> {
    /// The transition parameters from the last state of the module to match
    /// state j in the core pHMM.
    pub(crate) semilocal_params: &'a SemiLocalModule<T>,
    /// The transition parameters for consuming i bases within the module.
    pub(crate) domain_params:    PrecomputedDomainModule<T, S>,
}

impl<T: PhmmNumber, const S: usize> LocalModule<T, S> {
    /// Precomputes the transition probabilities for skipping any number of
    /// residues from the beginning of `seq`.
    ///
    /// This should only be called on modules placed at the beginning of the
    /// pHMM. See [`PrecomputedLocalModule`] for more details.
    pub(crate) fn precompute_begin_mod<Q: AsRef<[u8]>>(
        &self, seq: Q, mapping: &ByteIndexMap<S>,
    ) -> PrecomputedLocalModule<'_, T, S> {
        PrecomputedLocalModule {
            domain_params:    self.domain_params.precompute_begin_mod(seq, mapping),
            semilocal_params: &self.semilocal_params,
        }
    }

    /// Precomputes the transition probabilities for skipping any number of
    /// residues from the end of `seq`.
    ///
    /// This should only be called on modules placed at the end of the pHMM. See
    /// [`PrecomputedLocalModule`] for more details.
    pub(crate) fn precompute_end_mod<Q: AsRef<[u8]>>(
        &self, seq: Q, mapping: &ByteIndexMap<S>,
    ) -> PrecomputedLocalModule<'_, T, S> {
        PrecomputedLocalModule {
            domain_params:    self.domain_params.precompute_end_mod(seq, mapping),
            semilocal_params: &self.semilocal_params,
        }
    }
}

impl<T: PhmmNumber, const S: usize> PrecomputedLocalModule<'_, T, S> {
    /// Gets the score for skipping `i` residues in the query and
    /// entering/exiting layer `j`.
    ///
    /// Specifically, the score is for:
    ///
    /// - Skipping the first `i` residues in the query and then entering
    ///   directly into layer `j` (when this module is placed at the beginning
    ///   of the [`CorePhmm`])
    /// - Exiting early from layer `j` then skipping the last `i` residues in
    ///   the query (when this module is placed at the end of the [`CorePhmm`])
    pub(crate) fn get_score(&self, i: impl QueryIndex, j: impl PhmmIndex) -> T {
        self.domain_params.get_score(i) + self.semilocal_params.get_score(j)
    }
}

impl<T: PhmmNumber, const S: usize> DomainModule<T, S> {
    /// Constructs a [`DomainModule`] where all transitions are free, and the
    /// only penalty is the emission parameters for any skipped bases in the
    /// query.
    ///
    /// This will not produce a true probabilistic model where the probabilities
    /// sum to 1. However, it can be used for simplicity and ease of
    /// interpretation.
    #[inline]
    #[must_use]
    pub fn no_penalty(background_emission: EmissionParams<T, S>) -> Self {
        Self {
            start_to_insert: T::ZERO,
            insert_to_insert: T::ZERO,
            insert_to_end: T::ZERO,
            start_to_end: T::ZERO,
            background_emission,
        }
    }
}

/// A trait unifying the different modules providing semilocal-style parameters.
#[allow(dead_code)]
pub(crate) trait SemiLocalParams<T> {
    /// Returns a reference to the semilocal parameters of the module.
    #[must_use]
    fn semilocal_params(&self) -> &SemiLocalModule<T>;
}

impl<T> SemiLocalParams<T> for SemiLocalModule<T> {
    fn semilocal_params(&self) -> &SemiLocalModule<T> {
        self
    }
}

impl<T, const S: usize> SemiLocalParams<T> for LocalModule<T, S> {
    fn semilocal_params(&self) -> &SemiLocalModule<T> {
        &self.semilocal_params
    }
}

/// A trait unifying the different modules providing domain-style parameters.
#[allow(dead_code)]
pub(crate) trait DomainParams<T, const S: usize> {
    /// Returns a reference to the domain parameters of the module.
    #[must_use]
    fn domain_params(&self) -> &DomainModule<T, S>;
}

impl<T, const S: usize> DomainParams<T, S> for DomainModule<T, S> {
    fn domain_params(&self) -> &DomainModule<T, S> {
        self
    }
}

impl<T, const S: usize> DomainParams<T, S> for LocalModule<T, S> {
    fn domain_params(&self) -> &DomainModule<T, S> {
        &self.domain_params
    }
}
