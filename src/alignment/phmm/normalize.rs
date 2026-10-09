use crate::alignment::phmm::{
    DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
    at_least_two::VecAtLeast2,
    components::{CorePhmm, EmissionParams, LayerParams, TransitionParams},
    indexing::GetLayerMut,
    modules::{DomainModule, LocalModule, SemiLocalModule},
    relative_sum_probs, relative_sum_probs_slice,
    state::{PhmmState, PhmmStateArr, PhmmStateOrModule, PhmmStateOrModuleArr},
};

impl<T: PhmmNumber> PhmmStateArr<T> {
    /// Normalizes the parameters so that they correspond to a probability
    /// distribution (their corresponding probabilities add to 1).
    ///
    /// If all probabilities are 0, then normalization does not occur.
    fn normalize(&mut self) {
        self.0 = normalize(self.0);
    }
}

impl<T: PhmmNumber> PhmmStateOrModuleArr<T> {
    /// Normalizes the parameters so that they correspond to a probability
    /// distribution (their corresponding probabilities add to 1).
    ///
    /// If all probabilities are 0, then normalization does not occur.
    fn normalize(&mut self) {
        self.0 = normalize(self.0);
    }
}

impl<T: PhmmNumber, const S: usize> EmissionParams<T, S> {
    /// Normalizes the parameters so that they correspond to a probability
    /// distribution (their corresponding probabilities add to 1).
    ///
    /// If all probabilities are 0, then normalization does not occur.
    fn normalize(&mut self) {
        *self = EmissionParams::from_array(normalize(*self.as_array()));
    }
}

impl<T: PhmmNumber> TransitionParams<T> {
    /// Normalizes the parameters so that all sets of [`exiting_params`]
    /// correspond to a probability distribution (their corresponding
    /// probabilities add to 1).
    ///
    /// If all probabilities for any [`exiting_params`] are 0, then
    /// normalization does not occur.
    ///
    /// [`exiting_params`]: TransitionParams::exiting_params
    fn normalize(&mut self) {
        for state in PhmmState::VARIANTS {
            let mut vals = self.exiting_params(state);
            vals.normalize();
            self.set_exiting_params(state, vals);
        }
    }

    /// Normalizes the parameters so that all sets of [`exiting_params`]
    /// correspond to a probability distribution (their corresponding
    /// probabilities add to 1), and so that exiting the pHMM from a match state
    /// is considered.
    ///
    /// If all probabilities for any set are 0, then normalization does not
    /// occur for that set.
    ///
    /// [`exiting_params`]: TransitionParams::exiting_params
    fn normalize_with_exit(&mut self, exit_param: &mut T) {
        let mut vals = self.exiting_params(PhmmState::Match).with_module_val(*exit_param);
        vals.normalize();
        *exit_param = vals[PhmmStateOrModule::Module];
        self.set_exiting_params(PhmmState::Match, vals.into());

        let mut vals = self.exiting_params(PhmmState::Insert);
        vals.normalize();
        self.set_exiting_params(PhmmState::Insert, vals);

        let mut vals = self.exiting_params(PhmmState::Delete);
        vals.normalize();
        self.set_exiting_params(PhmmState::Delete, vals);
    }
}

impl<T: PhmmNumber, const S: usize> LayerParams<T, S> {
    /// Normalizes all sets of probabilities within the [`LayerParams`] so that
    /// they correspond to probability distributions (their corresponding
    /// probabilities add to 1). This assumes that exiting early from the match
    /// state is not allowed.
    ///
    /// If all probabilities in any set are 0, then normalization does not
    /// occur for that set.
    fn normalize_no_exit(&mut self) {
        self.transition.normalize();
        self.emission_match.normalize();
        self.emission_insert.normalize();
    }

    /// Normalizes all sets of probabilities within the [`LayerParams`] so that
    /// they correspond to probability distributions (their corresponding
    /// probabilities add to 1). The probability for exiting early from the
    /// match state is passed as `exit_param`.
    ///
    /// If all probabilities in any set are 0, then normalization does not
    /// occur for that set.
    fn normalize_with_exit(&mut self, exit_param: &mut T) {
        self.transition.normalize_with_exit(exit_param);
        self.emission_match.normalize();
        self.emission_insert.normalize();
    }
}

impl<T: PhmmNumber, const S: usize> DomainModule<T, S> {
    /// Normalizes all sets of probabilities within the [`DomainModule`] so that
    /// they correspond to probability distributions (their corresponding
    /// probabilities add to 1).
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    fn normalize(&mut self) {
        self.background_emission.normalize();

        let [new_start_to_insert, new_start_to_end] = normalize([self.start_to_insert, self.start_to_end]);
        self.start_to_insert = new_start_to_insert;
        self.start_to_end = new_start_to_end;

        let [new_insert_to_insert, new_insert_to_end] = normalize([self.insert_to_insert, self.insert_to_end]);
        self.insert_to_insert = new_insert_to_insert;
        self.insert_to_end = new_insert_to_end;
    }
}

impl<T: PhmmNumber> SemiLocalModule<T> {
    /// Normalizes the probabilities for entering the core pHMM from the
    /// [`SemiLocalModule`] so that they correspond to a probability
    /// distribution (the corresponding probabilities add to 1).
    ///
    /// If all probabilities are 0, then normalization does not occur.
    fn normalize_start(&mut self) {
        normalize_slice(self.as_slice_mut());
    }
}

impl<T: PhmmNumber, const S: usize> LocalModule<T, S> {
    /// Normalizes all sets of probabilities within the [`LocalModule`] at the
    /// start of a pHMM so that they correspond to probability distributions
    /// (their corresponding probabilities add to 1).
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    fn normalize_start(&mut self) {
        self.domain_params.normalize();
        self.semilocal_params.normalize_start();
    }

    /// Normalizes all sets of probabilities within the [`LocalModule`] at the
    /// end of a pHMM so that they correspond to probability distributions
    /// (their corresponding probabilities add to 1).
    ///
    /// The semilocal parameters are normalized in
    /// [`CorePhmm::normalize_with_exit`], not in this method.
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    fn normalize_end(&mut self) {
        self.domain_params.normalize();
    }
}

impl<T: PhmmNumber, const S: usize> CorePhmm<T, S> {
    /// Normalizes all sets of probabilities within the [`CorePhmm`] so that
    /// they correspond to probability distributions (their corresponding
    /// probabilities add to 1). This assumes that exiting early from the match
    /// states is not allowed.
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    fn normalize_no_exit(&mut self) {
        for layer in self.layers_mut() {
            layer.normalize_no_exit();
        }
    }

    /// Normalizes all sets of probabilities within the [`CorePhmm`] so that
    /// they correspond to probability distributions (their corresponding
    /// probabilities add to 1). The probabilities for exiting early from each
    /// layer is given by `exit_params`.
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    fn normalize_with_exit(&mut self, exit_params: &mut VecAtLeast2<T>) {
        for (layer, exit_param) in self.layers_mut().iter_mut().zip(exit_params.iter_mut()) {
            layer.normalize_with_exit(exit_param);
        }

        // Set exit parameter from End state to probability 1
        *exit_params.last_mut() = T::ZERO;
    }
}

impl<T: PhmmNumber, const S: usize> GlobalPhmm<'_, T, S> {
    /// Normalizes all sets of probabilities within the [`GlobalPhmm`] so that
    /// they correspond to probability distributions (their corresponding
    /// probabilities add to 1).
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    pub fn normalize(&mut self) {
        self.core.normalize_no_exit();
    }
}

impl<T: PhmmNumber, const S: usize> LocalPhmm<'_, T, S> {
    /// Normalizes all sets of probabilities within the [`LocalPhmm`] so that
    /// they correspond to probability distributions (their corresponding
    /// probabilities add to 1).
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    pub fn normalize(&mut self) {
        self.core.normalize_with_exit(&mut self.end.semilocal_params.0);
        self.begin.normalize_start();
        self.end.normalize_end();
    }
}

impl<T: PhmmNumber, const S: usize> DomainPhmm<'_, T, S> {
    /// Normalizes all sets of probabilities within the [`DomainPhmm`] so that
    /// they correspond to probability distributions (their corresponding
    /// probabilities add to 1). This includes editing the modules at the start
    /// and end of the pHMM.
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    pub fn normalize(&mut self) {
        self.core.normalize_no_exit();
        self.begin.normalize();
        self.end.normalize();
    }
}

impl<T: PhmmNumber, const S: usize> SemiLocalPhmm<'_, T, S> {
    /// Normalizes all sets of probabilities within the [`SemiLocalPhmm`] so
    /// that they correspond to probability distributions (their corresponding
    /// probabilities add to 1). This includes editing the modules at the start
    /// and end of the pHMM.
    ///
    /// If all probabilities in any set are 0, then normalization does not occur
    /// for that set.
    pub fn normalize(&mut self) {
        self.core.normalize_with_exit(&mut self.end.0);
        self.begin.normalize_start();
    }
}

/// Normalizes an array of [`PhmmNumber`] so that they correspond to a
/// probability distribution (their corresponding probabilities add to 1). This
/// includes editing the modules at the start and end of the pHMM.
///
/// If all probabilities are 0, then normalization does not occur. Any
/// probabilities above 1 are clamped to 1 before normalization.
#[must_use]
fn normalize<T: PhmmNumber, const N: usize>(mut vals: [T; N]) -> [T; N] {
    // Clamping isn't expected to be necessary based on our parser, but it
    // doesn't hurt being careful! Normalization isn't called in a hot loop
    for val in &mut vals {
        *val = val.clamp_valid_prob();
    }

    let (relative_sum, minimum) = relative_sum_probs(vals);

    let sum_probs = relative_sum + minimum;

    // Check if the sum of the probabilities is 0, which assuming validity
    // requirements are met, this means all probabilities are 0 and
    // normalization shouldn't occur.
    if sum_probs == T::INFINITY {
        return vals;
    }

    for val in &mut vals {
        // val-sum_probs = (val-minimum)-relative_sum helps to reduce
        // cancellation and improve numeric stability. relative_sum is bounded
        // between -ln(N) and 0, so it is better to apply it after subtracting
        // off minimum (and thus normalizing the scale a bit)
        *val = (*val - minimum) - relative_sum;
    }
    vals
}

/// Similar to [`normalize`] but without a known compile-time size, and
/// operating in-place.
///
/// Perfer [`normalize`] for better edge-case handling when possible.
fn normalize_slice<T: PhmmNumber>(vals: &mut [T]) {
    // Clamping isn't expected to be necessary based on our parser, but it
    // doesn't hurt being careful! Normalization isn't called in a hot loop
    for val in vals.iter_mut() {
        *val = val.clamp_valid_prob();
    }

    let (relative_sum, minimum) = relative_sum_probs_slice(vals);

    let sum_probs = relative_sum + minimum;

    // Check if the sum of the probabilities is 0, which assuming validity
    // requirements are met, this means all probabilities are 0 and
    // normalization shouldn't occur.
    if sum_probs == T::INFINITY {
        return;
    }

    for val in vals {
        // val-sum_probs = (val-minimum)-relative_sum helps to reduce
        // cancellation and improve numeric stability. relative_sum is bounded
        // between -ln(N) and 0, so it is better to apply it after subtracting
        // off minimum (and thus normalizing the scale a bit)
        *val = (*val - minimum) - relative_sum;
    }
}
