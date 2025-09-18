#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;

const ALPHAINV: usize = 13;

pub struct MethodDSampler<'a, I: Iterator> {
    reader: I,
    rng: &'a mut Xoshiro256StarStar,
    remaining_samples: usize,
    remaining_population: usize,
    alphainv: usize,
    use_method_a: bool,
    vprime: f32,
    max_start_index: usize,
    rem_samples_inv: f32,
    threshold: usize,
}

impl<'a, I: Iterator> MethodDSampler<'a, I> {
    /// Create a new [`MethodDSampler`] iterator that wraps and randomly
    /// subsamples another iterator of known size.
    ///
    /// # Errors
    ///
    /// The sample target size must be smaller than the total population,
    /// otherwise, an error is returned.
    pub fn new(reader: I, target: usize, total_items: usize, rng: &'a mut Xoshiro256StarStar) -> std::io::Result<Self> {
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

        Ok(MethodDSampler {
            reader,
            remaining_samples,
            remaining_population,
            max_start_index,
            rem_samples_inv,
            vprime,
            threshold: ALPHAINV * target,
            use_method_a: false,
            alphainv: ALPHAINV,
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

    /// Given a skip to return, update the remaining population and sample, and
    /// update any cached values.
    ///
    /// `None` is never returned; an option is used for convenience in
    /// [`next_skip`].
    ///
    /// [`next_skip`]: MethodDSampler::next_skip
    fn yield_skip(&mut self, skip: usize, rem_samples_min1_inv: f32) -> usize {
        // Decrement due to skip, and due to the item yielded
        self.remaining_population -= skip + 1;
        self.remaining_samples -= 1;

        // Update cached values
        self.rem_samples_inv = rem_samples_min1_inv;
        self.max_start_index -= skip;
        self.threshold -= self.alphainv;
        skip
    }

    /// Similar to [`next_skip`], but uses Algorithm A which is faster when
    /// `remaining_samples` is large relative to `remaining_population`.
    ///
    /// Rather than using rejection sampling (Algorithm D) which has overhead,
    /// Algorithm A uses inverse transform sampling. This requires checking each
    /// value of `skips` from 0 until the final value returned, so when
    /// `remaining_samples` is large (and hence skips tend to be small), this
    /// may be more efficient.
    ///
    /// [`next_skip`]: MethodDSampler::next_skip
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

impl<T, I> Iterator for MethodDSampler<'_, I>
where
    I: Iterator<Item = T>,
{
    type Item = I::Item;
    /// Calculates the number of skips and then skips it for the iterator inside
    fn next(&mut self) -> Option<Self::Item> {
        // makes sure to validate all remaining reads even if the target number
        // of reads is reached
        let skip = self.next_skip()?;
        self.reader.nth(skip)
    }
}