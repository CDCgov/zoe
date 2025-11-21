//! ## Subsampling
//!
//! *Zoe* provides the ability to efficiently and randomly downsample iterators
//! (with and without a known size). The following methods are provided:
//!
//! - Vitter's [**Method
//!   D**](https://www.ittc.ku.edu/~jsv/Papers/Vit87.RandomSampling.pdf):
//!   [`downsample_known_size`] or [`SkipSampler::new`]
//!
//!   This method operates by iterating through the collection, and calculating
//!   random values used to skip over multiple elements as it iterates through
//!   the collection. This method is very efficient, but requires knowing the
//!   length of the collection.
//!
//! - **Reservoir sampling** (or Vitter's [Method
//!   L](https://dl.acm.org/doi/pdf/10.1145/198429.198435)):
//!   [`downsample_reservoir`]
//!
//!   This method uses a reservoir to hold downsampled elements, and moves
//!   through an iterator of elements, randomly inserting elements into the
//!   reservoir. Reservoir sampling requires allocating all the elements in the
//!   reservoir, but has the advantage of not requiring knowledge of the length
//!   of the iterator.
//!
//! - **Bernoulli Sampling**: [`downsample_bernoulli`]
//!
//!   This method keeps each element with a given probability. Unlike reservoir
//!   sampling or *Method D*, the number of yielded elements may fluctuate.
//!
//! All methods require a random number generator, for which *Zoe* uses
//! `rand_xoshiro`'s [`Xoshiro256StarStar`].
//!
//! ### Example
//!
//! The following example shows downsampling of a `Vec<usize>`. Since the
//! iterator over `Vec<usize>` implements [`ExactSizeIterator`], the
//! [`downsample_known_size`] method can be used. Otherwise,
//! [`SkipSampler::new`] would have to be manually called.
//!
//! ```
//! # use rand_xoshiro::{Xoshiro256StarStar, rand_core::SeedableRng};
//! use zoe::iter_utils::sampling::DownsampleKnownSize;
//!
//! let target: usize = 4;
//!
//! // Hint: use from_os_rng() normally.
//! let seed = 42;
//! let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
//!
//! let nums: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
//!
//! let downsampled = nums.iter().downsample_known_size(&mut rng, target).unwrap();
//!
//! assert_eq!(downsampled.count(), target);
//! ```
//!
//! The following example shows downsampling of a `Vec<usize>`, but first, a
//! filter is applied, which leaves an iterator of unknown size, requiring use
//! of reservoir sampling.
//!
//! ```
//! # use rand_xoshiro::{Xoshiro256StarStar, rand_core::SeedableRng};
//! use zoe::iter_utils::sampling::downsample_reservoir;
//!
//! let target: usize = 4;
//!
//! // Hint: use from_os_rng() normally.
//! let seed = 42;
//! let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
//!
//! let nums: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
//! let filtered = nums.iter().filter(|num| *num % 2 == 0);
//!
//! let downsampled = downsample_reservoir(filtered, &mut rng, target);
//! assert_eq!(downsampled.len(), target);
//! ```
//!
//! <div class="warning note">
//!
//! **Note**
//!
//! You must enable the *rand* feature in your `Cargo.toml` to use these
//! functions. This feature is enabled by default.
//!
//! </div>
//!
//! [`downsample_known_size`]: DownsampleKnownSize::downsample_known_size
//! [`downsample_bernoulli`]: DownsampleBernoulli::downsample_bernoulli

#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
use rand::{Rng, seq::IndexedMutRandom};
use rand_xoshiro::Xoshiro256StarStar;

/// An iterator that efficiently downsampler from another iterator using *Method
/// D*.
pub struct SkipSampler<'a, I: Iterator> {
    iterator:             I,
    rng:                  &'a mut Xoshiro256StarStar,
    remaining_samples:    usize,
    remaining_population: usize,
    use_method_a:         bool,
    vprime:               f32,
    max_start_index:      usize,
    rem_samples_inv:      f32,
    threshold:            usize,
}

impl<'a, I: Iterator> SkipSampler<'a, I> {
    const ALPHAINV: usize = 13;

    /// Creates a new [`SkipSampler`] iterator that wraps and randomly
    /// subsamples another iterator of known size using *Method D*.
    ///
    /// ## Errors
    ///
    /// The sample target size must be smaller than the total population.
    pub fn new(iterator: I, target: usize, total_items: usize, rng: &'a mut Xoshiro256StarStar) -> std::io::Result<Self> {
        if target > total_items {
            return Err(std::io::Error::other(format!(
                "The sample target size {target} should not be larger than the population size {total_items}!"
            )));
        }

        let remaining_samples = target;
        let remaining_population = total_items;
        let max_start_index = remaining_population - remaining_samples + 1;
        let rem_samples_inv = (remaining_samples as f32).recip();
        let vprime = rng.random::<f32>().powf(rem_samples_inv);

        Ok(SkipSampler {
            iterator,
            remaining_samples,
            remaining_population,
            max_start_index,
            rem_samples_inv,
            vprime,
            threshold: Self::ALPHAINV * target,
            use_method_a: false,
            rng,
        })
    }

    fn next_skip(&mut self) -> Option<usize> {
        if self.remaining_samples == 0 {
            return None;
        }

        if self.use_method_a {
            return Some(self.next_skip_method_a());
        }

        if self.threshold >= self.remaining_population || self.remaining_samples == 1 {
            self.use_method_a = true;
            return Some(self.next_skip_method_a());
        }

        let rem_samples_min1_inv = ((self.remaining_samples - 1) as f32).recip();

        loop {
            let mut x = self.remaining_population as f32 * (1.0 - self.vprime);
            let mut skip = x as usize;

            while skip >= self.max_start_index {
                self.vprime = self.rng.random::<f32>().powf(self.rem_samples_inv);
                x = self.remaining_population as f32 * (1.0 - self.vprime);
                skip = x as usize;
            }

            let y1 = (self.rng.random::<f32>() * self.remaining_population as f32 / self.max_start_index as f32)
                .powf(rem_samples_min1_inv);
            self.vprime = y1
                * (1.0 - x / self.remaining_population as f32)
                * (self.max_start_index as f32 / (self.max_start_index - skip) as f32);

            if self.vprime <= 1.0 {
                return Some(self.yield_skip(skip, rem_samples_min1_inv));
            }

            let mut y2 = 1.0;
            let (mut bottom, limit) = if self.remaining_samples - 1 > skip {
                (
                    self.remaining_population - self.remaining_samples,
                    self.remaining_population - skip,
                )
            } else {
                (self.remaining_population - skip - 1, self.max_start_index)
            };

            for top in (limit..self.remaining_population).rev() {
                y2 *= top as f32 / bottom as f32;
                bottom -= 1;
            }

            if self.remaining_population as f32 / (self.remaining_population as f32 - x)
                >= y1 * y2.powf(rem_samples_min1_inv)
            {
                self.vprime = self.rng.random::<f32>().powf(rem_samples_min1_inv);
                return Some(self.yield_skip(skip, rem_samples_min1_inv));
            }

            self.vprime = self.rng.random::<f32>().powf(self.rem_samples_inv);
        }
    }

    /// Given a skip to return, updates the remaining population and sample, and
    /// update any cached values.
    fn yield_skip(&mut self, skip: usize, rem_samples_min1_inv: f32) -> usize {
        // Decrement due to skip, and due to the item yielded
        self.remaining_population -= skip + 1;
        self.remaining_samples -= 1;

        // Update cached values
        self.rem_samples_inv = rem_samples_min1_inv;
        self.max_start_index -= skip;
        self.threshold -= Self::ALPHAINV;
        skip
    }

    /// Similar to [`next_skip`], but uses *Algorithm A* which is faster when
    /// `remaining_samples` is large relative to `remaining_population`.
    ///
    /// Rather than using rejection sampling (*Method D*) which has overhead,
    /// *Algorithm A* uses inverse transform sampling. This requires checking
    /// each value of `skips` from 0 until the final value returned, so when
    /// `remaining_samples` is large (and hence skips tend to be small), this
    /// may be more efficient.
    ///
    /// [`next_skip`]: SkipSampler::next_skip
    fn next_skip_method_a(&mut self) -> usize {
        if self.remaining_samples > 1 {
            let mut top = self.remaining_population - self.remaining_samples;
            let mut skip = 0;
            let variate = self.rng.random::<f32>();

            // Find largest quotient of falling factorials (quot) greater than
            // random variate
            let mut quot = top as f32 / self.remaining_population as f32;
            while quot > variate {
                skip += 1;
                top -= 1;
                self.remaining_population -= 1;
                quot *= top as f32 / self.remaining_population as f32;
            }

            self.remaining_population -= 1;
            self.remaining_samples -= 1;
            skip
        } else {
            // The last remaining sample is picked uniformly from those
            // remaining
            let skip = self.rng.random_range(0..self.remaining_population);
            self.remaining_population -= skip + 1;
            self.remaining_samples = 0;
            skip
        }
    }
}

impl<I: Iterator> Iterator for SkipSampler<'_, I> {
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        let skip = self.next_skip()?;
        self.iterator.nth(skip)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        // Lower bound of 0 is due to case where target is larger than the
        // original iterator size
        (0, Some(self.remaining_samples))
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        let mut inner_n = self.next_skip()?;
        for _ in 0..n {
            // Add 1 to consume the element that would've been yielded, then add
            // the next skip
            inner_n += 1 + self.next_skip()?;
        }
        self.iterator.nth(inner_n)
    }
}

/// An iterator providing sampling using the Bernoulli method.
pub struct BernoulliSampler<'a, I: Iterator> {
    iterator: I,
    prob:     f32,
    rng:      &'a mut Xoshiro256StarStar,
}

impl<I> Iterator for BernoulliSampler<'_, I>
where
    I: Iterator,
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let item = self.iterator.next()?;
            if self.rng.random::<f32>() < self.prob {
                return Some(item);
            }
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, self.iterator.size_hint().1)
    }
}

/// Performs reservoir sampling on an iterator.
///
/// This returns an owned [`Vec`] that is a randomly downsampled subset of an
/// iterator (typically of unknown length). If the length of the iterator is
/// known, consider using [`downsample_known_size`].
///
/// [`downsample_known_size`]: DownsampleKnownSize::downsample_known_size
pub fn downsample_reservoir<I, T>(iter: I, rng: &mut Xoshiro256StarStar, target: usize) -> Vec<T>
where
    I: Iterator<Item = T>, {
    let mut reservoir = Vec::with_capacity(target);
    let n_sample = target;

    // Initial calculation of W and S
    let mut r = rng.random::<f32>();
    let mut w = (r.ln() / n_sample as f32).exp();
    r = rng.random::<f32>();
    let mut s = (r.ln() / (1.0 - w).ln()).floor() as usize;

    for (i, sample) in iter.enumerate() {
        // Initialize the reservoir
        if i < n_sample {
            reservoir.push(sample);
        }
        // s=0 case, finished skipping, make swap
        else if s == 0 {
            // Insert into a random position in reservoir
            if let Some(slot) = reservoir.choose_mut(rng) {
                *slot = sample;
            }

            // Recalculate S and W
            r = rng.random::<f32>();
            w *= (r.ln() / n_sample as f32).exp();
            r = rng.random::<f32>();
            s = (r.ln() / (1.0 - w).ln()).floor() as usize;
        }
        // s>0 case, keep skipping
        else {
            s -= 1;
        }
    }
    reservoir
}

/// Extension trait for [`Iterator`] that allows downsampling via the Bernoulli
/// method.
pub trait DownsampleBernoulli: Iterator + Sized {
    /// Downsamples the iterator using the Bernoulli method (keeping each item
    /// with probability `prob`).
    #[inline]
    fn downsample_bernoulli(self, prob: f32, rng: &mut Xoshiro256StarStar) -> BernoulliSampler<'_, Self> {
        BernoulliSampler {
            iterator: self,
            prob,
            rng,
        }
    }
}

impl<I: Iterator> DownsampleBernoulli for I {}

/// Extension trait for [`ExactSizeIterator`] that allows for *Method D*
/// sampling to be conducted on an iterator of known size.
pub trait DownsampleKnownSize: ExactSizeIterator + Sized {
    /// Creates a new [`SkipSampler`] iterator that wraps and randomly
    /// subsamples an [`ExactSizeIterator`] using *Method D*.
    ///
    /// If your iterator does not implement [`ExactSizeIterator`], then
    /// [`SkipSampler::new`] will need to be called directly.
    ///
    /// ## Errors
    ///
    /// The sample target size must be smaller than the total population.
    fn downsample_known_size(self, rng: &mut Xoshiro256StarStar, target: usize) -> std::io::Result<SkipSampler<'_, Self>>;
}

impl<I> DownsampleKnownSize for I
where
    I: ExactSizeIterator,
{
    #[inline]
    fn downsample_known_size(self, rng: &mut Xoshiro256StarStar, target: usize) -> std::io::Result<SkipSampler<'_, I>> {
        let total_items = self.len();
        SkipSampler::new(self, target, total_items, rng)
    }
}
