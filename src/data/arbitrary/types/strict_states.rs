use crate::{
    alignment::{AlignOp, AlignmentStates},
    data::{
        arbitrary::{ArbitrarySpecs, UnstructuredExt},
        cigar::Ciglet,
    },
};
use arbitrary::{Result, Unstructured};
use std::num::{NonZero, NonZeroUsize};

/// Specifications for generating an arbitrary [`AlignmentStates`] subject to
/// the conditions that the operations are only `M`, `D`, or `I` and the
/// increments of each [`Ciglet`] are nonzero. This allows for increased control
/// over the generated alignment, such as enforcing an exact total query length
/// or reference length.
#[derive(Clone, Default, Debug)]
pub struct StrictAlignmentStatesSpecs {
    // TODO: Could add avoid_repeated_ops which could condense after formation,
    // but then it might be possible to overflow usize during merging, so
    // another path for rearranging ciglets would be necessary. And this could
    // also fail, leading to unsolvable situation.
    /// Requires the number of residues consumed in the query to equal the
    /// given value.
    pub query_inc: Option<usize>,

    /// Requires the number of residues consumed in the reference to equal the
    /// given value.
    pub ref_inc: Option<usize>,

    /// Requires the sum of all increments to equal the given value.
    pub total_inc: Option<usize>,

    /// Ensures the number of residues consumed in the query is at most the
    /// given value.
    pub max_query_inc: Option<usize>,

    /// Ensures the number of residues consumed in the reference is at most the
    /// given value.
    pub max_ref_inc: Option<usize>,

    /// Ensures the sum of all increments is at most the given value.
    pub max_total_inc: Option<usize>,
}

/// Determines the `M`, `I`, and `D` increments for fixed query and reference
/// increments.
///
/// ## Panics
///
/// This will panic unless the following are true (for the case where
/// `max_total_inc` is set):
///
/// - `query_inc <= max_total_inc`
/// - `ref_inc <= max_total_inc`
fn strict_increments_from_query_ref(
    u: &mut Unstructured<'_>, query_inc: usize, ref_inc: usize, max_total_inc: Option<usize>,
) -> Result<(usize, usize, usize)> {
    let max_match_inc = query_inc.min(ref_inc);

    // Let excess_indel_inc measure how much the increments strays from a
    // match_inc-maximizing set of increments. Specifically, each unit of
    // excess_indel_inc represents replacing one match increment with a pair of
    // insertion and deletion increments.

    // The maximum value for this is max_match_inc by definition when
    // max_total_inc is None. Otherwise, it is also bounded above so that
    // match_inc+insert_inc+delete_inc <= max_total_inc. Substituting with
    // definitions below and rearranging gives excess_indel_inc <= max_total_inc
    // + max_match_inc - query_inc - ref_inc, which is excess_indel_inc <=
    // max_total_inc - query_inc.max(ref_inc)
    let max_excess_indel_inc = if let Some(max_total_inc) = max_total_inc {
        assert!(query_inc <= max_total_inc, "Must have query_inc <= max_total_inc");
        assert!(ref_inc <= max_total_inc, "Must have ref_inc <= max_total_inc");
        max_match_inc.min(max_total_inc - query_inc.max(ref_inc))
    } else {
        max_match_inc
    };

    let excess_indel_inc = u.int_in_range(0..=max_excess_indel_inc)?;
    let match_inc = max_match_inc - excess_indel_inc;
    let insert_inc = query_inc - match_inc;
    let delete_inc = ref_inc - match_inc;

    Ok((match_inc, insert_inc, delete_inc))
}

/// Determines the `M`, `I`, and `D` increments for fixed query and total
/// increments.
///
/// ## Panics
///
/// This function panics if:
///
/// - `query_inc > total_inc`
/// - `max_ref_inc` is set and `total_inc - query_inc > max_ref_inc`
fn strict_increments_from_query_total(
    u: &mut Unstructured<'_>, query_inc: usize, total_inc: usize, max_ref_inc: Option<usize>,
) -> Result<(usize, usize, usize)> {
    let delete_inc = total_inc.checked_sub(query_inc).expect("Must have query_inc <= total_inc");
    let min_insert_inc = max_ref_inc.map_or(0, |max_ref_inc| total_inc.saturating_sub(max_ref_inc));

    assert!(min_insert_inc <= query_inc, "Must have total_inc - query_inc <= max_ref_inc");

    let insert_inc = u.int_in_range(min_insert_inc..=query_inc)?;
    let match_inc = query_inc - insert_inc;

    Ok((match_inc, insert_inc, delete_inc))
}

/// Determines the `M`, `I`, and `D` increments for fixed reference and total
/// increments.
///
/// ## Panics
///
/// This function panics if:
///
/// - `ref_inc > total_inc`
/// - `max_query_inc` is set and `total_inc - ref_inc > max_query_inc`
fn strict_increments_from_ref_total(
    u: &mut Unstructured<'_>, ref_inc: usize, total_inc: usize, max_query_inc: Option<usize>,
) -> Result<(usize, usize, usize)> {
    let insert_inc = total_inc.checked_sub(ref_inc).expect("Must have ref_inc <= total_inc");
    let min_delete_inc = max_query_inc.map_or(0, |max_query_inc| total_inc.saturating_sub(max_query_inc));

    assert!(min_delete_inc <= ref_inc, "Must have total_inc - ref_inc <= max_query_inc");

    let delete_inc = u.int_in_range(min_delete_inc..=ref_inc)?;
    let match_inc = ref_inc - delete_inc;

    Ok((match_inc, insert_inc, delete_inc))
}

/// Determines the `M`, `I`, and possibly `D` increments for a fixed query
/// increment.
///
/// ## Panics
///
/// This function panics if `query_inc > max_total_inc` when `max_total_inc` is
/// set.
fn strict_increments_from_query(
    u: &mut Unstructured<'_>, query_inc: usize, max_ref_inc: Option<usize>, max_total_inc: Option<usize>,
) -> Result<(usize, usize, Option<usize>)> {
    if let Some(max_total_inc) = max_total_inc {
        assert!(query_inc <= max_total_inc, "Must have query_inc <= max_total_inc");
    }

    // Based on match_inc + delete_inc <= max_ref_inc, we derive insert_inc >=
    // query_inc - max_ref_inc + delete_inc, and since delete_inc has no
    // minimum, we get a relaxed insert_inc >= query_inc - max_ref_inc
    let min_insert_inc = max_ref_inc.map_or(0, |max_ref_inc| query_inc.saturating_sub(max_ref_inc));
    let insert_inc = u.int_in_range(min_insert_inc..=query_inc)?;
    let match_inc = query_inc - insert_inc;

    let delete_inc = if max_ref_inc.is_some() || max_total_inc.is_some() {
        let mut max_delete_inc = usize::MAX;
        if let Some(max_ref_inc) = max_ref_inc {
            max_delete_inc = max_delete_inc.min(max_ref_inc - match_inc);
        }
        if let Some(max_total_inc) = max_total_inc {
            max_delete_inc = max_delete_inc.min(max_total_inc - query_inc);
        }
        Some(u.int_in_range(0..=max_delete_inc)?)
    } else {
        None
    };

    Ok((match_inc, insert_inc, delete_inc))
}

/// Determines the `M`, possibly `I`, and `D` increments for a fixed reference
/// increment.
///
/// ## Panics
///
/// This function panics if `ref_inc > max_total_inc` when `max_total_inc` is
/// set.
fn strict_increments_from_ref(
    u: &mut Unstructured<'_>, ref_inc: usize, max_query_inc: Option<usize>, max_total_inc: Option<usize>,
) -> Result<(usize, Option<usize>, usize)> {
    if let Some(max_total_inc) = max_total_inc {
        assert!(ref_inc <= max_total_inc, "Must have ref_inc <= max_total_inc");
    }

    // Based on match_inc + insert_inc <= max_query_inc, we derive delete_inc >=
    // ref_inc - max_query_inc + insert_inc, and since insert_inc has no
    // minimum, we get a relaxed delete_inc >= ref_inc - max_query_inc
    let min_delete_inc = max_query_inc.map_or(0, |max_query_inc| ref_inc.saturating_sub(max_query_inc));
    let delete_inc = u.int_in_range(min_delete_inc..=ref_inc)?;
    let match_inc = ref_inc - delete_inc;

    let insert_inc = if max_query_inc.is_some() || max_total_inc.is_some() {
        let mut max_insert_inc = usize::MAX;
        if let Some(max_query_inc) = max_query_inc {
            max_insert_inc = max_insert_inc.min(max_query_inc - match_inc);
        }
        if let Some(max_total_inc) = max_total_inc {
            max_insert_inc = max_insert_inc.min(max_total_inc - ref_inc);
        }
        Some(u.int_in_range(0..=max_insert_inc)?)
    } else {
        None
    };

    Ok((match_inc, insert_inc, delete_inc))
}

/// Determines the `M`, `I`, and `D` increments for a fixed total increment.
///
/// ## Panics
///
/// Panics if `max_ref_inc` and `max_query_inc` are both set but their sum is
/// less than `total_inc`.
fn strict_increments_from_total(
    u: &mut Unstructured<'_>, total_inc: usize, max_query_inc: Option<usize>, max_ref_inc: Option<usize>,
) -> Result<(usize, usize, usize)> {
    let min_insert_inc = max_ref_inc.map_or(0, |max_ref_inc| total_inc.saturating_sub(max_ref_inc));
    let min_delete_inc = max_query_inc.map_or(0, |max_query_inc| total_inc.saturating_sub(max_query_inc));

    assert!(
        min_insert_inc <= total_inc - min_delete_inc,
        "Must have max_ref_inc + max_query_inc >= total_inc"
    );

    let delete_inc = u.int_in_range(min_delete_inc..=total_inc - min_insert_inc)?;
    let insert_inc = u.int_in_range(min_insert_inc..=total_inc - delete_inc)?;
    let match_inc = total_inc - delete_inc - insert_inc;

    Ok((match_inc, insert_inc, delete_inc))
}

/// Possibly determines the `M`, `I`, and `D` increments when only maximums are
/// given.
fn strict_increments_from_maximums(
    u: &mut Unstructured<'_>, max_query_inc: Option<usize>, max_ref_inc: Option<usize>, max_total_inc: Option<usize>,
) -> Result<(Option<usize>, Option<usize>, Option<usize>)> {
    let increments = match (max_query_inc, max_ref_inc, max_total_inc) {
        (None, None, None) => (None, None, None),
        (Some(max_query_inc), None, None) => {
            let query_inc = u.int_in_range(0..=max_query_inc)?;
            // Validity: max_total_inc is not set
            let (match_inc, insert_inc, delete_inc) = strict_increments_from_query(u, query_inc, None, None)?;
            (Some(match_inc), Some(insert_inc), delete_inc)
        }
        (None, Some(max_ref_inc), None) => {
            let ref_inc = u.int_in_range(0..=max_ref_inc)?;
            // Validity: max_total_inc is not set
            let (match_inc, insert_inc, delete_inc) = strict_increments_from_ref(u, ref_inc, None, None)?;
            (Some(match_inc), insert_inc, Some(delete_inc))
        }
        (None, None, Some(max_total_inc)) => {
            let total_inc = u.int_in_range(0..=max_total_inc)?;
            // Validity: max_query_inc and max_ref_inc are not set
            let (match_inc, insert_inc, delete_inc) = strict_increments_from_total(u, total_inc, None, None)?;
            (Some(match_inc), Some(insert_inc), Some(delete_inc))
        }
        (Some(max_query_inc), Some(max_ref_inc), None) => {
            let query_inc = u.int_in_range(0..=max_query_inc)?;
            let ref_inc = u.int_in_range(0..=max_ref_inc)?;
            // Validity: max_total_inc is not set
            let (match_inc, insert_inc, delete_inc) = strict_increments_from_query_ref(u, query_inc, ref_inc, None)?;
            (Some(match_inc), Some(insert_inc), Some(delete_inc))
        }
        (Some(max_query_inc), None, Some(max_total_inc)) => {
            let total_inc = u.int_in_range(0..=max_total_inc)?;
            let query_inc = u.int_in_range(0..=max_query_inc.min(total_inc))?;
            // Validity: max_ref_inc is not set, and a panic is desired if
            // query_inc > total_inc
            let (match_inc, insert_inc, delete_inc) = strict_increments_from_query_total(u, query_inc, total_inc, None)?;
            (Some(match_inc), Some(insert_inc), Some(delete_inc))
        }
        (None, Some(max_ref_inc), Some(max_total_inc)) => {
            let total_inc = u.int_in_range(0..=max_total_inc)?;
            let ref_inc = u.int_in_range(0..=max_ref_inc.min(total_inc))?;
            // Validity: max_query_inc is not set, and a panic is desired if
            // ref_inc > total_inc
            let (match_inc, insert_inc, delete_inc) = strict_increments_from_ref_total(u, ref_inc, total_inc, None)?;
            (Some(match_inc), Some(insert_inc), Some(delete_inc))
        }
        (Some(max_query_inc), Some(max_ref_inc), Some(max_total_inc)) => {
            let total_inc = u.int_in_range(0..=max_total_inc.min(max_query_inc.saturating_add(max_ref_inc)))?;
            let query_inc = u.int_in_range(total_inc.saturating_sub(max_ref_inc)..=max_query_inc.min(total_inc))?;
            let ref_inc = u.int_in_range(total_inc - query_inc..=max_ref_inc.min(total_inc))?;

            let delete_inc = total_inc - query_inc;
            let insert_inc = total_inc - ref_inc;
            let match_inc = query_inc - insert_inc;
            (Some(match_inc), Some(insert_inc), Some(delete_inc))
        }
    };

    Ok(increments)
}

/// Splits a known operation total into arbitrary positive ciglet increments,
/// or generates arbitrary positive increments when there is no known total.
fn strict_ciglet_increments(u: &mut Unstructured<'_>, total_inc: Option<usize>) -> Result<Vec<usize>> {
    if let Some(total_inc) = total_inc {
        if total_inc == 0 {
            return Ok(Vec::new());
        }

        let num_ciglets = u.arbitrary_len::<AlignOp>()?.clamp(1, total_inc);
        u.integer_composition(total_inc, num_ciglets)
    } else {
        u.arbitrary_iter::<NonZeroUsize>()?.map(|res| res.map(NonZero::get)).collect()
    }
}

impl<'a> ArbitrarySpecs<'a> for StrictAlignmentStatesSpecs {
    type Output = AlignmentStates;

    /// Generates an arbitrary [`AlignmentStates`] conforming to the given
    /// specifications, as well as all operations being only `M`, `D`, or `I`
    /// and all increments of each [`Ciglet`] being nonzero.
    ///
    /// ## Errors
    ///
    /// Any errors from the underlying [`arbitrary`] calls are propagated.
    ///
    /// ## Panics
    ///
    /// This panics if the exact and maximum increment fields are inconsistent,
    /// or if the three exact increment fields cannot be fulfilled.
    ///
    /// [`arbitrary`]: arbitrary::Arbitrary::arbitrary
    #[inline]
    fn make_arbitrary(&self, u: &mut Unstructured<'a>) -> Result<Self::Output> {
        // Validate that exact increments satisfy the corresponding max
        // increments if set
        if let Some(query_inc) = self.query_inc
            && let Some(max_query_inc) = self.max_query_inc
            && query_inc > max_query_inc
        {
            panic!("Must have query_inc <= max_query_inc")
        }
        if let Some(ref_inc) = self.ref_inc
            && let Some(max_ref_inc) = self.max_ref_inc
            && ref_inc > max_ref_inc
        {
            panic!("Must have ref_inc <= max_ref_inc")
        }
        if let Some(total_inc) = self.total_inc
            && let Some(max_total_inc) = self.max_total_inc
            && total_inc > max_total_inc
        {
            panic!("Must have total_inc <= max_total_inc")
        }

        let (match_inc, insert_inc, delete_inc) = match (self.query_inc, self.ref_inc, self.total_inc) {
            (Some(query_inc), Some(ref_inc), Some(total_inc)) => {
                let delete_inc = total_inc.checked_sub(query_inc).expect("Must have query_inc <= total_inc");
                let insert_inc = total_inc.checked_sub(ref_inc).expect("Must have ref_inc <= total_inc");
                let match_inc = query_inc
                    .checked_sub(insert_inc)
                    .expect("Must have query_inc + ref_inc >= total_inc");

                (Some(match_inc), Some(insert_inc), Some(delete_inc))
            }
            (Some(query_inc), Some(ref_inc), None) => {
                let (match_inc, insert_inc, delete_inc) =
                    strict_increments_from_query_ref(u, query_inc, ref_inc, self.max_total_inc)?;

                (Some(match_inc), Some(insert_inc), Some(delete_inc))
            }
            (Some(query_inc), None, Some(total_inc)) => {
                let (match_inc, insert_inc, delete_inc) =
                    strict_increments_from_query_total(u, query_inc, total_inc, self.max_ref_inc)?;

                (Some(match_inc), Some(insert_inc), Some(delete_inc))
            }
            (None, Some(ref_inc), Some(total_inc)) => {
                let (match_inc, insert_inc, delete_inc) =
                    strict_increments_from_ref_total(u, ref_inc, total_inc, self.max_query_inc)?;

                (Some(match_inc), Some(insert_inc), Some(delete_inc))
            }
            (Some(query_inc), None, None) => {
                let (match_inc, insert_inc, delete_inc) =
                    strict_increments_from_query(u, query_inc, self.max_ref_inc, self.max_total_inc)?;
                (Some(match_inc), Some(insert_inc), delete_inc)
            }
            (None, Some(ref_inc), None) => {
                let (match_inc, insert_inc, delete_inc) =
                    strict_increments_from_ref(u, ref_inc, self.max_query_inc, self.max_total_inc)?;
                (Some(match_inc), insert_inc, Some(delete_inc))
            }
            (None, None, Some(total_inc)) => {
                let (match_inc, insert_inc, delete_inc) =
                    strict_increments_from_total(u, total_inc, self.max_query_inc, self.max_ref_inc)?;

                (Some(match_inc), Some(insert_inc), Some(delete_inc))
            }
            (None, None, None) => {
                strict_increments_from_maximums(u, self.max_query_inc, self.max_ref_inc, self.max_total_inc)?
            }
        };

        let match_incs = strict_ciglet_increments(u, match_inc)?;
        let insert_incs = strict_ciglet_increments(u, insert_inc)?;
        let delete_incs = strict_ciglet_increments(u, delete_inc)?;

        let match_and_ins_states = u.interleave(
            match_incs.into_iter().map(|inc| Ciglet { inc, op: AlignOp::M }),
            insert_incs.into_iter().map(|inc| Ciglet { inc, op: AlignOp::I }),
        )?;

        let states = u.interleave(
            match_and_ins_states.into_iter(),
            delete_incs.into_iter().map(|inc| Ciglet { inc, op: AlignOp::D }),
        )?;

        Ok(AlignmentStates(states))
    }
}
