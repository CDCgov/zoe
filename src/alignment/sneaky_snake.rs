/// A pre-alignment filter based upon [Sneaky Snake] (1) which will reject
/// queries exceeding a given edit distance from the reference.
///
/// The `threshold` is defined as a decimal between 0.0 and 1.0, representing a
/// proportion of the query length. For example, a threshold of 0.2 for a query
/// of length 100 will reject any query with more than 20 edits (insertions,
/// deletions, or substitutions) compared to the reference.
///
/// Returns `Some(true)` if number of edits between query and reference is
/// within the provided threshold. Returns `Some(false)` if number of edits
/// between query and reference exceeds the provided threshold. Returns `None`
/// if inputs are invalid.
///
/// ## Limitations
///
/// - The provided edit distance threshold must be greater than or equal to the
///   difference in length between the reference and the query. Returns `None`
///   otherwise.
/// - Although Sneaky Snake is designed for sequences of the same length, this
///   implementation is able to handle sequences of slightly different lengths
///   (following the above requirement that the length difference is less than
///   the edit distance threshold). The excess residues at either end of the
///   longer sequence are not penalized. More specifically, if the reference is
///   $\Delta$ residues longer than the query, the algorithm supports up to
///   $\Delta$ unpenalized residues at either end of the reference. Note that
///   this is not the same as semiglobal edit distance, since only $\Delta$
///   unpenalized residues are supported instead of arbitrarily many. Rather,
///   this algorithm is a relaxation of global edit distance, designed to not
///   overly penalize sequences with different lengths.
/// - The filter calculates an approximated edit distance value that is very
///   close to the actual edit distance. Its calculated edit distance is always
///   less than or equal to the actual edit distance and then compares to the
///   provided threshold.
///
/// ## Time Complexity
///
/// $N$ is the length of the shorter sequence and $E$ is the edit distance
/// threshold.
///
/// - When the true number of edits is close to the threshold, the worst case
///   time complexity approaches $ \mathcal{O}(N * E) $.
///
/// - When the sequences are nearly identical and when enough of the sequence
///   has been examined and cannot possibly get filtered out with the remaining
///   length, the best case time complexity approaches $ \Omega(N - E) $.
///
/// - Alternatively, when the sequences differ significantly at the beginning
///   causing early termination, the best case time complexity approaches $
///   \Omega(E^2) $.
///
/// ## Example
///
/// ```
/// # use zoe::alignment::sneaky_snake;
/// let reference = b"GGTGCAGAGCTC";
/// let query = b"GGTGAGAGTTGT";
/// let threshold = 0.25;
///
/// let result = sneaky_snake(reference, query, threshold);
/// assert_eq!(result, Some(true));
/// ```
///
/// ## Module Citations
///
/// 1. Alser, Mohammed; Shahroodi, Taha; Gómez-Luna, Juan; Alkan, Can; Mutlu,
///    Onur (2020). "SneakySnake: a fast and accurate universal genome
///    pre-alignment filter for CPUs, GPUs and FPGAs". *Bioinformatics*. 36
///    (22–23): 5282–5290.
///
/// [Sneaky Snake]: https://doi.org/10.1093/bioinformatics/btaa1015
#[must_use]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::doc_markdown
)]
pub fn sneaky_snake(reference: &[u8], query: &[u8], threshold: f32) -> Option<bool> {
    if !(0. ..=1.).contains(&threshold) {
        return None;
    }

    let edit_thresh = ((query.len() as f32) * threshold).round() as usize;

    if reference.is_empty() {
        return Some(query.len() <= edit_thresh);
    }
    if query.is_empty() {
        return Some(reference.len() <= edit_thresh);
    }
    if edit_thresh == query.len() {
        return Some(true);
    }
    if edit_thresh == 0 {
        return Some(reference == query);
    }

    // choose shorter string
    let (s1, s2, len_diff) = if reference.len() > query.len() {
        (query, reference, reference.len() - query.len())
    } else {
        (reference, query, query.len() - reference.len())
    };

    if len_diff > edit_thresh {
        return None;
    }

    let window = 2 * edit_thresh + 1;
    let diffpad_len = len_diff / 2;
    let mut obstacles = 0;
    let mut checkpoint = 0;

    // the sneaky snake will now solve the chip maze
    while checkpoint < s1.len() && obstacles <= edit_thresh && s1.len() - checkpoint > edit_thresh - obstacles {
        let mut last_col = checkpoint;
        // Microbenchmarks suggest that iterating over the rows sequentially is
        // in general faster than controlling the order of row access so middle
        // rows are checked first.
        for row in 0..window {
            for (col, s1_char) in s1.iter().enumerate().skip(checkpoint) {
                if let Some(s2_idx) = (col + row + diffpad_len).checked_sub(edit_thresh)
                    && s2.get(s2_idx) == Some(s1_char)
                {
                    if col == s1.len() - 1 || s1.len() - col - 1 <= edit_thresh - obstacles {
                        return Some(true);
                    }
                } else {
                    last_col = last_col.max(col);
                    break;
                }
            }
        }

        checkpoint = last_col + 1;
        obstacles += 1;
    }

    Some(obstacles <= edit_thresh)
}
