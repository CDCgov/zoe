use rand::{
    Rng,
    distr::{
        uniform::{SampleUniform, UniformSampler},
        weighted::{Error, Weight},
    },
};
use std::borrow::Borrow;

/// An allocation-free, one-off alternative to [`WeightedIndex`] for selecting
/// between options with different probabilities.
///
/// ## Parameters
///
/// - `N`: The number of options being selected between.
/// - `X`: The type specifying the weights.
/// - `R`: The random number generator used to pick the option.
///
/// ## Errors
///
/// - [`InvalidWeight`] if any weight is not at least zero (or is `NaN`).
/// - [`Overflow`] if the sum of the weights overflows.
/// - [`InsufficientNonZero`] if the sum of the weights is zero.
///
/// ## Panics
///
/// Panics if the number of options `N` is 0.
///
/// ## Acknowledgements
///
/// This contains a modified version of code from `rand`.
///
/// [`WeightedIndex`]: rand::distr::weighted::WeightedIndex
/// [`InvalidWeight`]: Error::InvalidWeight
/// [`Overflow`]: Error::Overflow
/// [`InsufficientNonZero`]: Error::InsufficientNonZero
#[allow(clippy::neg_cmp_op_on_partial_ord, reason = "same code as rand")]
pub fn sample_one_weighted<const N: usize, X, R>(rng: &mut R, weights: &[X; N]) -> Result<usize, Error>
where
    X: Weight + SampleUniform + PartialOrd,
    R: Rng + ?Sized, {
    const { assert!(N > 0, "At least one option must be possible when using sample_weighted") }

    let mut total_weight = weights[0].borrow().clone();

    let zero = X::ZERO;
    if !(total_weight >= zero) {
        return Err(Error::InvalidWeight);
    }

    // Only the first N-1 values are used
    let mut cumulative_weights = weights.clone();

    for (i, w) in weights[1..].iter().enumerate() {
        // Note that `!(w >= x)` is not equivalent to `w < x` for partially
        // ordered types due to NaNs which are equal to nothing.
        if !(w.borrow() >= &zero) {
            return Err(Error::InvalidWeight);
        }
        cumulative_weights[i] = total_weight.clone();

        if let Err(()) = total_weight.checked_add_assign(w.borrow()) {
            return Err(Error::Overflow);
        }
    }

    if total_weight == zero {
        return Err(Error::InsufficientNonZero);
    }
    let distr = X::Sampler::new(zero, total_weight.clone()).unwrap();

    let chosen_weight = distr.sample(rng);
    // Find the first item which has a weight *higher* than the chosen weight.
    Ok(cumulative_weights[..N - 1].partition_point(|w| w <= &chosen_weight))
}
