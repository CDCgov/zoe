use crate::alignment::{
    Alignment,
    phmm::{
        DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
        traverse::score_from_path::{ScoreVisitor, WithScore},
    },
};
use rand::Rng;

mod error;
mod visitor;

pub use error::*;
pub use visitor::*;

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
    /// ## Validity
    ///
    /// If the transition probabilities out of an insert state are zero, this
    /// implementation may enter an infinite loop. If they are too small, the
    /// sampling may run slowly.
    ///
    /// The pHMM parameters must be correctly specified at the beginning and end
    /// of the model, so that invalid transitions have a probability of 0.
    /// *Zoe*'s parsers automatically handle this.
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
    /// ## Validity
    ///
    /// If the transition probabilities out of an insert state are zero, this
    /// implementation may enter an infinite loop. If they are too small, the
    /// sampling may run slowly.
    ///
    /// The pHMM parameters must be correctly specified at the beginning and end
    /// of the model, so that invalid transitions have a probability of 0.
    /// *Zoe*'s parsers automatically handle this.
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
    /// ## Validity
    ///
    /// If the transition probabilities out of an insert state are zero, this
    /// implementation may enter an infinite loop. If they are too small, the
    /// sampling may run slowly.
    ///
    /// The pHMM parameters must be correctly specified at the beginning and end
    /// of the model, so that invalid transitions have a probability of 0.
    /// *Zoe*'s parsers automatically handle this.
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
    /// ## Validity
    ///
    /// If the transition probabilities out of an insert state are zero, this
    /// implementation may enter an infinite loop. If they are too small, the
    /// sampling may run slowly.
    ///
    /// The pHMM parameters must be correctly specified at the beginning and end
    /// of the model, so that invalid transitions have a probability of 0.
    /// *Zoe*'s parsers automatically handle this.
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
