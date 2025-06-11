//! Structs defining the modules that appear may appear at the beginning and/or
//! end of a [`CorePhmm`] in order to change the alignment mode

use crate::{
    alignment::phmm::{
        CorePhmm, EmissionParams,
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
    /// `DomainPrecomputed`).
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    pub(crate) fn get_begin_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        if inserted.is_empty() {
            self.start_to_end
        } else {
            self.start_to_insert
                + self.insert_to_end
                + (inserted
                    .iter()
                    .map(|x| self.background_emission[mapping.to_index(*x)])
                    .fold(T::ZERO, |acc, elem| acc + elem)
                    + T::cast_from(inserted.len() - 1) * self.insert_to_insert)
        }
    }

    /// Given the residues which were inserted, lazily compute the score for a
    /// [`DomainModule`] at the end of a pHMM (rather than using
    /// `DomainPrecomputed`).
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    pub(crate) fn get_end_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        if inserted.is_empty() {
            self.start_to_end
        } else {
            self.start_to_insert
                + self.insert_to_end
                + (inserted
                    .iter()
                    .rev()
                    .map(|x| self.background_emission[mapping.to_index(*x)])
                    .fold(T::ZERO, |acc, elem| acc + elem)
                    + T::cast_from(inserted.len() - 1) * self.insert_to_insert)
        }
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

    /// Gets the transition parameter stored within a [`SemiLocalModule`] for
    /// transitioning from/to a given [`PhmmIndex`].
    #[inline]
    #[must_use]
    pub(crate) fn get_score(&self, index: impl PhmmIndex) -> T {
        self.0[index.get_phmm_dp_index(self)]
    }
}
