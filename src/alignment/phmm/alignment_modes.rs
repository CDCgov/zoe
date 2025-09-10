//! Structs defining the modules that appear may appear at the beginning and/or
//! end of a [`CorePhmm`] in order to change the alignment mode

use crate::{
    alignment::phmm::{
        CorePhmm, EmissionParams, QueryIndex, QueryIndexable,
        indexing::{PhmmIndex, PhmmIndexable},
    },
    data::ByteIndexMap,
    math::Float,
};

/// A module placed before or after a [`CorePhmm`] to support local alignment
/// (aligning a subsequence to a portion of the model).
///
/// This module can skip arbitrarily many residues from the beginning/end of the
/// query sequence, and it can also jump to/from an arbitrary match statement in
/// the pHMM. Internally, this is a combination of a [`SemiLocalModule`] (for
/// supporting arbitrary jumps to/from the pHMM) and a [`DomainModule`] (for
/// skipping residues in the query).
#[derive(Debug)]
pub struct LocalModule<T, const S: usize> {
    /// The parameters for connecting this module with the match states of the
    /// pHMM.
    ///
    /// These are either the transition parameters into the match states of the
    /// pHMM (if the module is at the beginning) or the transition parameters
    /// out of the match states of the pHMM (if the module is at the end).
    pub external_params: SemiLocalModule<T>,
    /// The parameters within this module.
    ///
    /// These control how skipped residues in the query are scored.
    pub internal_params: DomainModule<T, S>,
}

/// A module placed before or after a [`CorePhmm`] to support semilocal
/// alignment (aligning the full sequence to a portion of the model).
///
/// This module can jump to/from an arbitrary match statement in the pHMM. Each
/// match state (and the BEGIN state and END state) have a transition parameter
/// associated with them.
#[derive(Debug)]
pub struct SemiLocalModule<T>(pub Vec<T>);

/// A module paced before or after a [`CorePhmm`] to support domain alignment
/// (aligning a subsequence to the full model).
///
/// This module can skip arbitrarily many residues from the beginning/end of the
/// query sequence.
#[derive(Debug)]
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

impl<T: Float, const S: usize> DomainModule<T, S> {
    /// Given the residues which were inserted, lazily compute the score for a
    /// [`DomainModule`] at the beginning of a pHMM (rather than using
    /// [`PrecomputedDomainModule`]).
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    pub(crate) fn get_begin_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        if inserted.is_empty() {
            self.start_to_end
        } else {
            // Special casing needed in case insert_to_insert is infinite,
            // causing a NAN to appear when multiplied by 0
            let insert_to_insert = if inserted.len() > 1 {
                T::cast_from(inserted.len() - 1) * self.insert_to_insert
            } else {
                T::ZERO
            };

            self.start_to_insert
                + self.insert_to_end
                + (inserted
                    .iter()
                    .map(|x| self.background_emission[mapping.to_index(*x)])
                    .fold(T::ZERO, |acc, elem| acc + elem)
                    + insert_to_insert)
        }
    }

    /// Given the residues which were inserted, lazily compute the score for a
    /// [`DomainModule`] at the end of a pHMM (rather than using
    /// [`PrecomputedDomainModule`]).
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    pub(crate) fn get_end_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        if inserted.is_empty() {
            self.start_to_end
        } else {
            // Special casing needed in case insert_to_insert is infinite,
            // causing a NAN to appear when multiplied by 0
            let insert_to_insert = if inserted.len() > 1 {
                T::cast_from(inserted.len() - 1) * self.insert_to_insert
            } else {
                T::ZERO
            };

            self.start_to_insert
                + self.insert_to_end
                + (inserted
                    .iter()
                    .rev()
                    .map(|x| self.background_emission[mapping.to_index(*x)])
                    .fold(T::ZERO, |acc, elem| acc + elem)
                    + insert_to_insert)
        }
    }

    /// Precomputes the transition probabilities for skipping any number of
    /// residues from the beginning of `seq`.
    ///
    /// This should only be called on modules placed at the beginning of the
    /// pHMM. See [`PrecomputedDomainModule`] for more details.
    fn precompute_begin_mod<Q: AsRef<[u8]>>(&self, seq: Q, mapping: &ByteIndexMap<S>) -> PrecomputedDomainModule<T, S> {
        let seq = seq.as_ref();
        let mut internal_params = vec![self.start_to_insert + self.insert_to_end; seq.len() + 1];
        internal_params[0] = self.start_to_end;

        let mut emissions_sum = T::ZERO;
        for (i, x_idx) in seq.iter().map(|x| mapping.to_index(*x)).enumerate() {
            emissions_sum += self.background_emission[x_idx];
            internal_params[i + 1] += emissions_sum + T::cast_from(i) * self.insert_to_insert;
        }

        PrecomputedDomainModule(internal_params)
    }

    /// Precomputes the transition probabilities for skipping any number of
    /// residues from the end of `seq`.
    ///
    /// This should only be called on modules placed at the end of the pHMM. See
    /// [`PrecomputedDomainModule`] for more details.
    fn precompute_end_mod<Q: AsRef<[u8]>>(&self, seq: Q, mapping: &ByteIndexMap<S>) -> PrecomputedDomainModule<T, S> {
        let seq = seq.as_ref();
        let mut internal_params = vec![self.start_to_insert + self.insert_to_end; seq.len() + 1];
        internal_params[seq.len()] = self.start_to_end;

        let mut emissions_sum = T::ZERO;
        for (i, x_idx) in seq.iter().rev().map(|x| mapping.to_index(*x)).enumerate() {
            emissions_sum += self.background_emission[x_idx];
            internal_params[seq.len() - i - 1] += emissions_sum + T::cast_from(i) * self.insert_to_insert;
        }

        PrecomputedDomainModule(internal_params)
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
    fn get_score(&self, i: impl QueryIndex) -> T {
        self.0[self.get_dp_index(i)]
    }
}

impl<T: Float, const S: usize> LocalModule<T, S> {
    /// Constructs a [`LocalModule`] where all transitions are free, and the
    /// only penalty is the emission parameters for any skipped bases in the
    /// query.
    ///
    /// This will not produce a true probabilistic model where the probabilities
    /// sum to 1. However, it can be used for simplicity and ease of
    /// interpretation.
    #[inline]
    #[must_use]
    pub fn no_penalty(core: &CorePhmm<T, S>, background_emission: EmissionParams<T, S>) -> Self {
        Self {
            external_params: SemiLocalModule::no_penalty(core),
            internal_params: DomainModule::no_penalty(background_emission),
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
    pub(crate) external_params: &'a SemiLocalModule<T>,
    /// The transition parameters for consuming i bases within the module.
    pub(crate) internal_params: PrecomputedDomainModule<T, S>,
}

impl<T: Float, const S: usize> LocalModule<T, S> {
    /// Precomputes the transition probabilities for skipping any number of
    /// residues from the beginning of `seq`.
    ///
    /// This should only be called on modules placed at the beginning of the
    /// pHMM. See [`PrecomputedLocalModule`] for more details.
    pub(crate) fn precompute_begin_mod<Q: AsRef<[u8]>>(
        &self, seq: Q, mapping: &ByteIndexMap<S>,
    ) -> PrecomputedLocalModule<'_, T, S> {
        PrecomputedLocalModule {
            internal_params: self.internal_params.precompute_begin_mod(seq, mapping),
            external_params: &self.external_params,
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
            internal_params: self.internal_params.precompute_end_mod(seq, mapping),
            external_params: &self.external_params,
        }
    }
}

impl<T: Float, const S: usize> PrecomputedLocalModule<'_, T, S> {
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
        self.internal_params.get_score(i) + self.external_params.get_score(j)
    }
}

impl<T: Float, const S: usize> DomainModule<T, S> {
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

impl<T: Float> SemiLocalModule<T> {
    /// Constructs a [`SemiLocalModule`] where all transitions are free.
    ///
    /// This will not produce a true probabilistic model where the probabilities
    /// sum to 1. However, it can be used for simplicity and ease of
    /// interpretation.
    #[inline]
    #[must_use]
    pub fn no_penalty<const S: usize>(core: &CorePhmm<T, S>) -> Self {
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
