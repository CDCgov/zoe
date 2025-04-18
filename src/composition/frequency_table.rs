/// Represents the frequency of u8 bytes within a sequence.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct FrequencyTable<const MIN: u8 = 0, const MAX: u8 = 255> {
    table:     [usize; 256],
    num_items: usize,
}

impl<const MIN: u8, const MAX: u8> FrequencyTable<MIN, MAX> {
    /// Initialize a new empty frequency table (with counts as zero).
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Generate a frequency table.
    #[inline]
    #[must_use]
    pub fn new_from<Q: AsRef<[u8]> + ?Sized>(seq: &Q) -> Self {
        let mut out = Self::new();
        out.tally_from(seq);
        out
    }

    /// Generate a frequency table, but assume all values are in range.
    #[inline]
    #[must_use]
    pub fn new_from_unchecked<Q: AsRef<[u8]> + ?Sized>(seq: &Q) -> Self {
        let mut out = Self::new();
        out.tally_from_unchecked(seq);
        out
    }

    /// Tally the counts from a sequence into an already existing frequency
    /// table.
    #[inline]
    pub fn tally_from<Q: AsRef<[u8]> + ?Sized>(&mut self, seq: &Q) {
        self.num_items += tally_from::<MIN, MAX, _>(seq, &mut self.table);
    }

    /// Tally the counts from a sequence into an already existing frequency
    /// table, but assume all values are in range.
    #[inline]
    pub fn tally_from_unchecked<Q: AsRef<[u8]> + ?Sized>(&mut self, seq: &Q) {
        self.num_items += tally_from_unchecked(seq, &mut self.table);
    }

    /// Calculate the median value in the frequency table.
    #[inline]
    #[must_use]
    pub fn get_median(&self) -> Option<f32> {
        get_median::<MIN>(&self.table, self.num_items)
    }

    /// Calculate the minimum, median, and maximum values from the frequency
    /// table. This version is more performant than running min, max, and median
    /// independently.
    #[inline]
    #[must_use]
    pub fn get_min_med_max(&self) -> Option<(u8, f32, u8)> {
        get_min_med_max::<MIN, MAX>(&self.table, self.num_items)
    }

    /// Returns the array of counts stored within the frequency table.
    #[inline]
    #[must_use]
    pub fn into_inner(self) -> [usize; 256] {
        self.table
    }
}

#[inline]
#[must_use]
pub(crate) fn tally_from<const MIN: u8, const MAX: u8, Q: AsRef<[u8]> + ?Sized>(seq: &Q, table: &mut [usize; 256]) -> usize {
    let mut num_items = 0;
    for &num in seq.as_ref() {
        table[num as usize] += 1;
        if (MIN..=MAX).contains(&num) {
            num_items += 1;
        }
    }
    num_items
}

#[inline]
#[must_use]
pub(crate) fn tally_from_unchecked<Q: AsRef<[u8]> + ?Sized>(seq: &Q, table: &mut [usize; 256]) -> usize {
    for &num in seq.as_ref() {
        table[num as usize] += 1;
    }
    seq.as_ref().len()
}

#[inline]
#[must_use]
pub(crate) fn get_median<const MIN: u8>(table: &[usize; 256], num_items: usize) -> Option<f32> {
    if num_items == 0 {
        return None;
    }

    // Consume half of our target total, rounding up
    let threshold = num_items.div_ceil(2);
    let mut index: u16 = MIN.into();
    let mut total = table[index as usize];

    while total < threshold {
        index += 1;
        total += table[index as usize];
    }

    // total > Floor(N/2), total ≥ Ceil(N/2)
    let threshold = num_items / 2;
    if total > threshold {
        Some(f32::from(index))

    // Else: Ceil(N/2) ≤ total ≤ Floor(N/2) where Floor(N/2) ≤ Ceil(N/2)
    // ⇒ Ceil(N/2) ≤ Floor(N/2) ⇒ Floor(N/2) = Ceil(N/2) = N/2 (even)
    // ⇒ total = N/2
    } else {
        let mut next_index = index + 1;
        while table[next_index as usize] == 0 {
            next_index += 1;
        }
        Some(f32::from(index + next_index) / 2.0)
    }
}

#[inline]
#[must_use]
pub(crate) fn get_min_med_max<const MIN: u8, const MAX: u8>(
    table: &[usize; 256], num_items: usize,
) -> Option<(u8, f32, u8)> {
    let mut min = None;
    for i in MIN..=MAX {
        let count = table[i as usize];
        if count > 0 {
            min = Some(i);
            break;
        }
    }

    let mut max = None;
    for i in (MIN..=MAX).rev() {
        let count = table[i as usize];
        if count > 0 {
            max = Some(i);
            break;
        }
    }

    // TODO: Reduce to a single pass through data
    Some((min?, get_median::<MIN>(table, num_items)?, max?))
}

impl<const MIN: u8, const MAX: u8> Default for FrequencyTable<MIN, MAX> {
    #[inline]
    fn default() -> Self {
        FrequencyTable {
            table:     [0; 256],
            num_items: 0,
        }
    }
}

#[cfg(test)]
mod test {
    #![allow(clippy::type_complexity)]
    use super::*;

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
            let freq_table = FrequencyTable::<0, 255>::new_from(&v);
            assert_eq!(freq_table.get_median(), expected);
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
            let freq_table = FrequencyTable::<222, 240>::new_from(&v);
            assert_eq!(freq_table.get_median(), expected);
        }
    }

    #[test]
    fn median_data_too_large() {
        let freq_table = FrequencyTable::<222, 240>::new_from(&[222, 245, 230, 233]);
        assert_eq!(freq_table.get_median(), Some(230.0));
    }

    #[test]
    fn median_data_too_small() {
        let freq_table = FrequencyTable::<222, 240>::new_from(&[222, 200, 230, 233]);
        assert_eq!(freq_table.get_median(), Some(230.0));
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
            let freq_table = FrequencyTable::<0, 255>::new_from(&v);
            assert_eq!(freq_table.get_min_med_max(), expected);
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
            let freq_table = FrequencyTable::<222, 240>::new_from(&v);
            assert_eq!(freq_table.get_min_med_max(), expected);
        }
    }

    #[test]
    fn med_min_max_too_small() {
        let freq_table = FrequencyTable::<222, 240>::new_from(&[222, 200, 230, 233]);
        assert_eq!(freq_table.get_median(), Some(230.0));
    }

    #[test]
    fn med_min_max_too_large() {
        let freq_table = FrequencyTable::<222, 240>::new_from(&[222, 245, 230, 233]);
        assert_eq!(freq_table.get_median(), Some(230.0));
    }
}
