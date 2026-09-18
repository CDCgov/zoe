//! Implementations for sampling sequences/alignments probabilistically from
//! pHMMs.

use crate::{
    alignment::{
        Alignment,
        phmm::{
            DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
            state::{PhmmState, PhmmStateArr, PhmmStateOrModule, PhmmStateOrModuleArr},
            traverse::score_from_path::{ScoreVisitor, WithScore},
        },
    },
    math::sample_one_weighted,
};
use rand::{
    Rng,
    distr::{
        uniform::SampleUniform,
        weighted::{Error, Weight},
    },
};

mod error;
mod visitor;

pub use error::*;
pub use visitor::*;

impl<X> PhmmStateArr<X> {
    /// Samples a state from the [`PhmmStateArr`], using the values as weights.
    ///
    /// ## Errors
    ///
    /// See [`sample_one_weighted`].
    pub fn sample<R>(&self, rng: &mut R) -> Result<PhmmState, Error>
    where
        X: Weight + SampleUniform + PartialOrd,
        R: Rng + ?Sized, {
        let idx = sample_one_weighted(rng, &self.0)?;
        Ok(PhmmState::VARIANTS[idx])
    }
}

impl<X> PhmmStateOrModuleArr<X> {
    /// Samples a state from the [`PhmmStateOrModuleArr`], using the values as
    /// weights.
    ///
    /// ## Errors
    ///
    /// See [`sample_one_weighted`].
    pub fn sample<R>(&self, rng: &mut R) -> Result<PhmmStateOrModule, Error>
    where
        X: Weight + SampleUniform + PartialOrd,
        R: Rng + ?Sized, {
        let idx = sample_one_weighted(rng, &self.0)?;
        Ok(PhmmStateOrModule::VARIANTS[idx])
    }
}

/// A pHMM-sampled sequence alongside its alignment information.
#[derive(Debug)]
pub struct SampledSequence<T> {
    /// The sampled query sequence.
    pub sequence:  Vec<u8>,
    /// The alignment used when generating the sequence.
    ///
    /// Note that this may not be the optimal alignment yielded by the Viterbi
    /// algorithm.
    pub alignment: Alignment<T>,
}

impl<T: PhmmNumber, const S: usize> GlobalPhmm<T, S> {
    /// Samples a sequence and corresponding [`Alignment`] from the
    /// [`GlobalPhmm`].
    ///
    /// To avoid pathological cases, insertion sizes are truncated to 1000.
    ///
    /// ## Errors
    ///
    /// Any errors while sampling (e.g., bad parameters, no path with non-zero
    /// probability found, or an invalid model) are returned as
    /// [`SamplingError`].
    pub fn sample(&self, rng: &mut impl Rng) -> Result<SampledSequence<T>, SamplingError<T, S>>
    where
        T: Clone + 'static, {
        let visitor = ScoreVisitor::new(SampleVisitor::new(rng));
        let output = self.traverse(visitor)?;

        let WithScore {
            output: SampledSequence { sequence, alignment },
            score,
        } = output;

        let alignment = Alignment {
            score,
            ref_range: alignment.ref_range,
            query_range: alignment.query_range,
            states: alignment.states,
            ref_len: alignment.ref_len,
            query_len: alignment.query_len,
        };

        Ok(SampledSequence { sequence, alignment })
    }
}

impl<T: PhmmNumber, const S: usize> DomainPhmm<T, S> {
    /// Samples a sequence and corresponding [`Alignment`] from the
    /// [`DomainPhmm`].
    ///
    /// To avoid pathological cases, insertion sizes are truncated to 1000.
    ///
    /// ## Errors
    ///
    /// Any errors while sampling (e.g., bad parameters, no path with non-zero
    /// probability found, or an invalid model) are returned as
    /// [`SamplingError`].
    pub fn sample(&self, rng: &mut impl Rng) -> Result<SampledSequence<T>, SamplingError<T, S>>
    where
        T: Clone + 'static, {
        let visitor = ScoreVisitor::new(SampleVisitor::new(rng));
        let output = self.traverse(visitor)?;

        let WithScore {
            output: SampledSequence { sequence, alignment },
            score,
        } = output;

        let alignment = Alignment {
            score,
            ref_range: alignment.ref_range,
            query_range: alignment.query_range,
            states: alignment.states,
            ref_len: alignment.ref_len,
            query_len: alignment.query_len,
        };

        Ok(SampledSequence { sequence, alignment })
    }
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<T, S> {
    /// Samples a sequence and corresponding [`Alignment`] from the
    /// [`SemiLocalPhmm`].
    ///
    /// To avoid pathological cases, insertion sizes are truncated to 1000.
    ///
    /// ## Errors
    ///
    /// Any errors while sampling (e.g., bad parameters, no path with non-zero
    /// probability found, or an invalid model) are returned as
    /// [`SamplingError`].
    pub fn sample(&self, rng: &mut impl Rng) -> Result<SampledSequence<T>, SamplingError<T, S>>
    where
        T: Clone + 'static, {
        let output = self.traverse(ScoreVisitor::new(SampleVisitor::new(rng)))?;

        let WithScore {
            output: SampledSequence { sequence, alignment },
            score,
        } = output;

        let alignment = Alignment {
            score,
            ref_range: alignment.ref_range,
            query_range: alignment.query_range,
            states: alignment.states,
            ref_len: alignment.ref_len,
            query_len: alignment.query_len,
        };

        Ok(SampledSequence { sequence, alignment })
    }
}

impl<T: PhmmNumber, const S: usize> LocalPhmm<T, S> {
    /// Samples a sequence and corresponding [`Alignment`] from the
    /// [`LocalPhmm`].
    ///
    /// To avoid pathological cases, insertion sizes are truncated to 1000.
    ///
    /// ## Errors
    ///
    /// Any errors while sampling (e.g., bad parameters, no path with non-zero
    /// probability found, or an invalid model) are returned as
    /// [`SamplingError`].
    pub fn sample(&self, rng: &mut impl Rng) -> Result<SampledSequence<T>, SamplingError<T, S>>
    where
        T: Clone + 'static, {
        let output = self.traverse(ScoreVisitor::new(SampleVisitor::new(rng)))?;

        let WithScore {
            output: SampledSequence { sequence, alignment },
            score,
        } = output;

        let alignment = Alignment {
            score,
            ref_range: alignment.ref_range,
            query_range: alignment.query_range,
            states: alignment.states,
            ref_len: alignment.ref_len,
            query_len: alignment.query_len,
        };

        Ok(SampledSequence { sequence, alignment })
    }
}
