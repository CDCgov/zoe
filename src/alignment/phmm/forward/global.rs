//! The implementation of the forward algorithm for [`GlobalPhmm`].

use crate::alignment::phmm::{
    GlobalPhmm, PhmmNumber,
    forward::{update_delete, update_insert},
    indexing::{GetLayer, GetMapping},
    state::PhmmState::{Delete, Insert, Match},
    sum_probs,
};

impl<T: PhmmNumber, const S: usize> GlobalPhmm<'_, T, S> {
    /// Computes the probability of a sequence arising from the pHMM.
    ///
    /// The probability $P(x)$ is returned in negative log space as
    /// $-\ln(P(x))$. This means that lower values are more probable. To convert
    /// this into a probability, use [`PhmmNumber::to_prob`].
    ///
    /// See [here](crate::alignment::phmm#alignment-algorithms) for more
    /// details.
    pub fn forward<Q: AsRef<[u8]>>(&self, seq: Q) -> T {
        let seq = seq.as_ref();

        let (end, layers) = self.layers().split_last();

        let query_dim = seq.len() + 1;

        let mut f_m = vec![T::INFINITY; query_dim];
        f_m[0] = T::ZERO;
        let mut f_i = vec![T::INFINITY; query_dim];
        let mut f_d = vec![T::INFINITY; query_dim];

        for layer in layers {
            let mut cur_m = f_m[0];

            for (i, x_idx) in seq.iter().map(|x| self.mapping().to_index(*x)).enumerate() {
                let match_val = cur_m;
                let delete_val = f_d[i];
                let insert_val = f_i[i];

                let match_score = layer.emission_match[x_idx]
                    + sum_probs([
                        match_val + layer.transition[(Match, Match)],
                        delete_val + layer.transition[(Delete, Match)],
                        insert_val + layer.transition[(Insert, Match)],
                    ]);

                cur_m = std::mem::replace(&mut f_m[i + 1], match_score);

                f_i[i + 1] = update_insert(layer, x_idx, match_val, delete_val, insert_val);
                f_d[i] = update_delete(layer, match_val, delete_val, insert_val);
            }

            let i = seq.len();
            f_d[i] = update_delete(layer, cur_m, f_d[i], f_i[i]);

            // Cannot enter a match state after BEGIN without consuming a base
            f_m[0] = T::INFINITY;
        }

        // Calculate last column of insertion table and previous states for END
        // state
        for (i, x_idx) in seq.iter().map(|x| self.mapping().to_index(*x)).enumerate() {
            f_i[i + 1] = update_insert(end, x_idx, f_m[i], f_d[i], f_i[i]);
        }

        let i = seq.len();

        sum_probs([
            f_m[i] + end.transition[(Match, Match)],
            f_d[i] + end.transition[(Delete, Match)],
            f_i[i] + end.transition[(Insert, Match)],
        ])
    }
}
