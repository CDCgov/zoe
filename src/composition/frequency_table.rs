/// Trait for creating frequency tables from bytes slices and performing other
/// operations on said table.
pub trait FrequencyTable: AsRef<[u8]> {
    const MIN: u8 = 0;
    const MAX: u8 = 255;

    /// Generate a frequency table for the given range.
    #[inline]
    #[must_use]
    fn get_frequency_table(&self) -> ([usize; 256], usize) {
        let raw_data = self.as_ref();
        let mut count_table = [0; 256];
        let mut valid_items = 0;

        for &num in raw_data {
            count_table[num as usize] += 1;
            if Self::MIN <= num && num <= Self::MAX {
                valid_items += 1;
            }
        }
        (count_table, valid_items)
    }

    /// Generate a frequency table but assume all values are in range.
    #[inline]
    #[must_use]
    fn get_frequency_table_unchecked(&self) -> ([usize; 256], usize) {
        let raw_data = self.as_ref();
        let mut count_table = [0; 256];

        for &num in raw_data {
            count_table[num as usize] += 1;
        }

        (count_table, raw_data.len())
    }

    /// Calculate the median value in the frequency table.
    #[inline]
    #[must_use]
    fn get_median(&self) -> Option<f32> {
        if self.as_ref().is_empty() {
            return None;
        }

        let (count_table, num_items) = Self::get_frequency_table(self);
        median_from_table(&count_table, num_items, Self::MIN)
    }

    /// Calculate the median value in the frequency table but assume all values
    /// are in `MIN..=MAX`.
    #[inline]
    #[must_use]
    fn get_median_unchecked(&self) -> Option<f32> {
        let data = self.as_ref();
        if data.is_empty() {
            return None;
        } else if data.len() == 1 {
            return data.first().map(|&b| b.into());
        }

        let (count_table, num_items) = Self::get_frequency_table_unchecked(self);
        median_from_table(&count_table, num_items, Self::MIN)
    }

    // Calculate the minimum, median, and maximum values from the frequency
    // table. This version is more performant than running min, max, and median
    // independently.
    #[inline]
    #[must_use]
    fn get_min_med_max(&self) -> Option<(u8, f32, u8)> {
        if self.as_ref().is_empty() {
            return None;
        }

        let (count_table, num_items) = Self::get_frequency_table(self);

        let mut min = None;
        for i in Self::MIN..=Self::MAX {
            let count = count_table[i as usize];
            if count > 0 {
                min = Some(i);
                break;
            }
        }

        let mut max = None;
        for i in (Self::MIN..=Self::MAX).rev() {
            let count = count_table[i as usize];
            if count > 0 {
                max = Some(i);
                break;
            }
        }

        Some((min?, median_from_table(&count_table, num_items, Self::MIN)?, max?))
    }

    // Calculate the minimum, median, and maximum values from the frequency
    // table. Assumes values are aleady in range for `MIN..=MAX`.
    #[inline]
    #[must_use]
    fn get_min_med_max_unchecked(&self) -> Option<(u8, f32, u8)> {
        let data = self.as_ref();
        if data.is_empty() {
            return None;
        } else if data.len() == 1 {
            return Some((data[0], data[0].into(), data[0]));
        }

        let (count_table, num_items) = Self::get_frequency_table_unchecked(self);

        let mut min = None;
        for i in Self::MIN..=Self::MAX {
            let count = count_table[i as usize];
            if count > 0 {
                min = Some(i);
                break;
            }
        }

        let mut max = None;
        for i in (Self::MIN..=Self::MAX).rev() {
            let count = count_table[i as usize];
            if count > 0 {
                max = Some(i);
                break;
            }
        }

        Some((min?, median_from_table(&count_table, num_items, Self::MIN)?, max?))
    }
}

/// Private function to get the median from a frequency table.
#[must_use]
#[inline]
fn median_from_table(count_table: &[usize; 256], num_items: usize, min: u8) -> Option<f32> {
    if num_items == 0 {
        return None;
    }

    // Consume half of our target total, rounding up
    let threshold = num_items.div_ceil(2);
    let mut index: u16 = min.into();
    let mut total = 0;

    while total < threshold {
        total += count_table[index as usize];
        index += 1;
    }
    index -= 1;

    // total > Floor(N/2), total ≥ Ceil(N/2)
    let threshold = num_items / 2;
    if total > threshold {
        Some(f32::from(index))

    // Else: Ceil(N/2) ≤ total ≤ Floor(N/2) where Floor(N/2) ≤ Ceil(N/2)
    // ⇒ Ceil(N/2) ≤ Floor(N/2) ⇒ Floor(N/2) = Ceil(N/2) = N/2 (even)
    // ⇒ total = N/2
    } else {
        let mut next_index = index + 1;
        while count_table[next_index as usize] == 0 {
            next_index += 1;
        }
        Some(f32::from(index + next_index) / 2.0)
    }
}

#[cfg(test)]
mod test {
    #![allow(clippy::type_complexity)]
    use super::*;

    struct BigFake(Vec<u8>);

    impl FrequencyTable for BigFake {
        const MIN: u8 = 0;
        const MAX: u8 = 255;
    }

    impl AsRef<[u8]> for BigFake {
        fn as_ref(&self) -> &[u8] {
            &self.0
        }
    }

    struct SmallFake(Vec<u8>);

    impl FrequencyTable for SmallFake {
        const MIN: u8 = 222;
        const MAX: u8 = 240;
    }

    impl AsRef<[u8]> for SmallFake {
        fn as_ref(&self) -> &[u8] {
            &self.0
        }
    }

    #[test]
    fn median_full() {
        let median_data: Vec<(Vec<u8>, Option<f32>)> = vec![
            (vec![0, 255, 1, 3], Some(2.0)),
            (vec![], None),
            (vec![24, 129], Some(76.5)),
            (vec![1, 255, 3], Some(3.0)),
            (vec![0, 0], Some(0.0)),
            (vec![0], Some(0.0)),
            (vec![255], Some(255.0)),
        ];

        for (v, expected) in median_data {
            let bf: BigFake = BigFake(v);
            assert_eq!(bf.get_median(), expected);
        }
    }

    #[test]
    fn median_small() {
        let median_data: Vec<(Vec<u8>, Option<f32>)> = vec![
            (vec![222, 240, 230, 233], Some(231.5)),
            (vec![], None),
            (vec![224, 229], Some(226.5)),
            (vec![231, 235, 230], Some(231.0)),
            (vec![222, 240], Some(231.0)),
            (vec![222], Some(222.0)),
            (vec![240], Some(240.0)),
        ];

        for (v, expected) in median_data {
            let bf: SmallFake = SmallFake(v);
            assert_eq!(bf.get_median(), expected);
        }
    }

    #[test]
    fn median_data_too_large() {
        let f = SmallFake(vec![222, 245, 230, 233]);
        assert_eq!(f.get_median(), Some(230.0));
    }

    #[test]
    fn median_data_too_small() {
        let f = SmallFake(vec![222, 200, 230, 233]);
        assert_eq!(f.get_median(), Some(230.0));
    }

    #[test]
    fn med_min_max_full() {
        let all_data: Vec<(Vec<u8>, Option<(u8, f32, u8)>)> = vec![
            (vec![0, 255, 1, 3], Some((0, 2.0, 255))),
            (vec![], None),
            (vec![24, 129], Some((24, 76.5, 129))),
            (vec![1, 255, 3], Some((1, 3.0, 255))),
            (vec![0, 0], Some((0, 0.0, 0))),
            (vec![0], Some((0, 0.0, 0))),
            (vec![255], Some((255, 255.0, 255))),
        ];

        for (v, expected) in all_data {
            let bf: BigFake = BigFake(v);
            assert_eq!(bf.get_min_med_max(), expected);
        }
    }

    #[test]
    fn med_min_max_small() {
        let median_data: Vec<(Vec<u8>, Option<(u8, f32, u8)>)> = vec![
            (vec![222, 240, 230, 233], Some((222, 231.5, 240))),
            (vec![], None),
            (vec![224, 229], Some((224, 226.5, 229))),
            (vec![231, 235, 230], Some((230, 231.0, 235))),
            (vec![222, 240], Some((222, 231.0, 240))),
            (vec![222], Some((222, 222.0, 222))),
            (vec![240], Some((240, 240.0, 240))),
        ];

        for (v, expected) in median_data {
            let bf: SmallFake = SmallFake(v);
            assert_eq!(bf.get_min_med_max(), expected);
        }
    }

    #[test]
    fn med_min_max_too_small() {
        let f = SmallFake(vec![222, 200, 230, 233]);
        assert_eq!(f.get_median(), Some(230.0));
    }

    #[test]
    fn med_min_max_too_large() {
        let f = SmallFake(vec![222, 245, 230, 233]);
        assert_eq!(f.get_median(), Some(230.0));
    }
}
