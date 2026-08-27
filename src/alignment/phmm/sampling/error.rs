use crate::{
    alignment::phmm::{
        InvalidModelError, PhmmNumber,
        indexing::DpIndex,
        state::{PhmmState, PhmmStateOrExit, PhmmStateOrModule},
        traverse::{EndInsert, EndInsertExit, ModuleLocation},
    },
    data::{ByteIndexMap, err::GetCode},
};
use rand::seq::WeightError;
use std::{
    error::Error,
    fmt::{Debug, Display},
};

/// An error with the weights in a pHMM causing sampling to fail.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum ParamWeightError {
    /// No path could be found, due to all paths having probability 0.
    NoPathFound,
    /// An overflow occurred when sampling from the parameters.
    Overflow,
    /// One of the probabilities was invalid (e.g., NaN or negative).
    InvalidWeight,
    /// Some other unexpected error occurred.
    ///
    /// This is introduced for forward-compatibility with `rand`.
    Other(WeightError),
}

impl From<WeightError> for ParamWeightError {
    fn from(value: WeightError) -> ParamWeightError {
        match value {
            WeightError::InvalidWeight => ParamWeightError::InvalidWeight,
            WeightError::InsufficientNonZero => ParamWeightError::NoPathFound,
            WeightError::Overflow => ParamWeightError::Overflow,
            other => ParamWeightError::Other(other),
        }
    }
}

impl Display for ParamWeightError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParamWeightError::NoPathFound => write!(f, "No path could be found."),
            ParamWeightError::Overflow => {
                write!(f, "An overflow occurred when trying to sample the available parameters.")
            }
            ParamWeightError::InvalidWeight => write!(f, "One of the parameters was invalid."),
            ParamWeightError::Other(_) => write!(f, "An unexpected error when sampling occurred"),
        }
    }
}

impl Error for ParamWeightError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ParamWeightError::Other(e) => Some(e),
            _ => None,
        }
    }
}

impl GetCode for ParamWeightError {}

/// A parameter alongside a label, for use in [`ParamSamplingError`].
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct LabeledParam<T, L> {
    /// The label describing the parameter option (e.g., the state being
    /// transitioned into, or the symbol being emitted).
    pub label: L,
    /// The parameter indicating the likelihood of that choice.
    pub param: T,
}

/// A generic error for when sampling one of a small constant number of choices
/// fails.
///
/// This is used to represent failures to sample a transition or an emission.
/// This adds an additional level of context above [`ParamWeightError`] by
/// including the choices and parameters.
///
/// ## Parameters
///
/// - `T`: The parameter type.
/// - `L`: The label type (the options being selected).
/// - `N`: The number of options.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct ParamSamplingError<T, L, const N: usize> {
    /// The underlying problem with the weights.
    pub source: ParamWeightError,
    /// The labeled parameters.
    pub params: [LabeledParam<T, L>; N],
}

impl<T, L> LabeledParam<T, L>
where
    T: Copy,
    L: Copy,
{
    /// A helper function to build an array of [`LabeledParam`]
    fn new_arr<const N: usize>(labels: [L; N], params: [T; N]) -> [LabeledParam<T, L>; N] {
        std::array::from_fn(|i| LabeledParam {
            label: labels[i],
            param: params[i],
        })
    }
}

impl<T> ParamSamplingError<T, PhmmState, 3>
where
    T: Copy,
{
    /// Creates a new [`ParamSamplingError`] for sampling a [`PhmmState`].
    ///
    /// The parameters should be in the same order as [`PhmmState::VARIANTS`].
    pub(crate) fn new_state(e: WeightError, params: [T; 3]) -> Self {
        let source = ParamWeightError::from(e);
        let params = LabeledParam::new_arr(PhmmState::VARIANTS, params);
        Self { source, params }
    }
}

impl<T> ParamSamplingError<T, PhmmStateOrExit, 4>
where
    T: Copy,
{
    /// Creates a new [`ParamSamplingError`] for sampling a pHMM state or
    /// exiting.
    ///
    /// The parameters should be in the same order as
    /// [`PhmmStateOrModule::VARIANTS`], with exiting to the module being the
    /// last parameter.
    pub(crate) fn new_state_or_exit(e: WeightError, params: [T; 4]) -> Self {
        let source = ParamWeightError::from(e);
        let params = LabeledParam::new_arr(PhmmStateOrModule::VARIANTS.map(PhmmStateOrModule::display_exit), params);
        Self { source, params }
    }
}

impl<T> ParamSamplingError<T, EndInsert, 2>
where
    T: Copy,
{
    /// Creates a new [`ParamSamplingError`] for sampling whether to enter the
    /// END state or the final insert state from the last layer.
    pub(crate) fn new_end_or_insert(e: WeightError, end_param: T, insert_param: T) -> Self {
        let source = ParamWeightError::from(e);
        let params = LabeledParam::new_arr([EndInsert::End, EndInsert::Insert], [end_param, insert_param]);
        Self { source, params }
    }
}

impl<T> ParamSamplingError<T, EndInsertExit, 3>
where
    T: Copy,
{
    /// Creates a new [`ParamSamplingError`] for sampling whether to enter the
    /// END state, the final insert state, or exit early from the last match
    /// state.
    pub(crate) fn new_end_insert_or_exit(e: WeightError, end_param: T, insert_param: T, exit_param: T) -> Self {
        let source = ParamWeightError::from(e);
        let params = LabeledParam::new_arr(
            [EndInsertExit::End, EndInsertExit::Insert, EndInsertExit::Exit],
            [end_param, insert_param, exit_param],
        );
        Self { source, params }
    }
}

impl<T, const S: usize> ParamSamplingError<T, char, S>
where
    T: Copy,
{
    /// Creates a new [`ParamSamplingError`] for sampling an emission from an
    /// alphabet.
    ///
    /// The parameters should be in the same order as `map.byte_keys()`.
    pub(crate) fn new_emission(e: WeightError, params: &[T; S], map: &ByteIndexMap<S>) -> Self {
        let source = ParamWeightError::from(e);
        let params = LabeledParam::new_arr(map.byte_keys().map(|byte| byte as char), *params);
        Self { source, params }
    }
}

impl<T, L, const N: usize> Display for ParamSamplingError<T, L, N>
where
    T: PhmmNumber,
    L: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Could not sample from parameters:")?;

        for LabeledParam { label, param } in &self.params {
            let prob = param.to_prob::<f32>();
            write!(f, "\n{label}: param={param}, prob={prob}")?;
        }

        Ok(())
    }
}

impl<T, L, const N: usize> Error for ParamSamplingError<T, L, N>
where
    T: PhmmNumber + Debug,
    L: Display + Debug,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.source)
    }
}

impl<T, L, const N: usize> GetCode for ParamSamplingError<T, L, N> {
    fn get_code(&self) -> i32 {
        self.source.get_code()
    }
}

/// An error representing a sampling error where the only non-zero probability
/// transition is a self-loop, causing an infinite loop during sampling.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct ParamTrapError<T, L, const N: usize> {
    /// The labeled parameters.
    pub params: [LabeledParam<T, L>; N],
}

impl<T> ParamTrapError<T, PhmmState, 3>
where
    T: Copy,
{
    /// Creates a new [`ParamTrapError`] for sampling a [`PhmmState`].
    ///
    /// The parameters should be in the same order as [`PhmmState::VARIANTS`].
    pub(crate) fn new_state(params: [T; 3]) -> Self {
        let params = LabeledParam::new_arr(PhmmState::VARIANTS, params);
        Self { params }
    }
}

impl<T> ParamTrapError<T, EndInsert, 2>
where
    T: Copy,
{
    /// Creates a new [`ParamTrapError`] for sampling whether to enter an INSERT
    /// state or END state.
    pub(crate) fn new_end_or_insert(end_param: T, insert_param: T) -> Self {
        let params = LabeledParam::new_arr([EndInsert::End, EndInsert::Insert], [end_param, insert_param]);
        Self { params }
    }
}

impl<T, L, const N: usize> Display for ParamTrapError<T, L, N>
where
    T: PhmmNumber,
    L: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "The only possible transition is a self-loop, causing an infinite loop during sampling:"
        )?;

        for LabeledParam { label, param } in &self.params {
            let prob = param.to_prob::<f32>();
            write!(f, "\n{label}: param={param}, prob={prob}")?;
        }

        Ok(())
    }
}

impl<T, L, const N: usize> Error for ParamTrapError<T, L, N>
where
    T: PhmmNumber + Debug,
    L: Display + Debug,
{
}

impl<T, L, const N: usize> GetCode for ParamTrapError<T, L, N> {}

/// A sampling error caused by bad transition parameters.
///
/// In the error stack, this is a transparent enum adding no additional context.
/// It's display implementation immediately forwards to the contained error
/// variant.
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum TransitionParamError<T, L, const N: usize> {
    InvalidWeights(ParamSamplingError<T, L, N>),
    Trap(ParamTrapError<T, L, N>),
}

impl<T, L, const N: usize> Display for TransitionParamError<T, L, N>
where
    T: PhmmNumber,
    L: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TransitionParamError::InvalidWeights(err) => write!(f, "{err}"),
            TransitionParamError::Trap(err) => write!(f, "{err}"),
        }
    }
}

impl<T, L, const N: usize> Error for TransitionParamError<T, L, N>
where
    T: PhmmNumber + Debug,
    L: Display + Debug,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            TransitionParamError::InvalidWeights(err) => err.source(),
            TransitionParamError::Trap(err) => err.source(),
        }
    }
}

impl<T, L, const N: usize> GetCode for TransitionParamError<T, L, N> {
    fn get_code(&self) -> i32 {
        match self {
            TransitionParamError::InvalidWeights(err) => err.get_code(),
            TransitionParamError::Trap(err) => err.get_code(),
        }
    }
}

impl<T, L, const N: usize> From<ParamSamplingError<T, L, N>> for TransitionParamError<T, L, N> {
    fn from(value: ParamSamplingError<T, L, N>) -> Self {
        TransitionParamError::InvalidWeights(value)
    }
}

impl<T, L, const N: usize> From<ParamTrapError<T, L, N>> for TransitionParamError<T, L, N> {
    fn from(value: ParamTrapError<T, L, N>) -> Self {
        TransitionParamError::Trap(value)
    }
}

/// An enum unifying the different state-sampling errors that could occur.
///
/// For ease of error handling, [`TransitionParamError`] is used only in cases
/// where a self-loop is possible (i.e., it might be an insert state being
/// transitioned out of). Otherwise, the basic [`ParamSamplingError`] is used.
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum SampleFromStateErrorKind<T> {
    /// A sampling error while choosing whether to transition to a match
    /// state, insert state, or delete state.
    State(TransitionParamError<T, PhmmState, 3>),
    /// A sampling error while choosing whether to transition to a match
    /// state, insert state, delete state, or exit the pHMM early.
    StateOrExit(ParamSamplingError<T, PhmmStateOrExit, 4>),
    /// A sampling error while choosing whether to transition from the last
    /// layer to either the END state or the final insert state.
    EndOrInsert(TransitionParamError<T, EndInsert, 2>),
    /// A sampling error while choosing whether to transition from the last
    /// match layer to either the END state, the final insert state, or exiting
    /// directly.
    EndInsertOrExit(ParamSamplingError<T, EndInsertExit, 3>),
}

/// A sampling error while choosing the next state to enter.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SampleFromStateError<T> {
    /// The cause of the sampling error.
    pub kind:    SampleFromStateErrorKind<T>,
    /// The state that is being transitioned out of.
    pub exiting: PhmmState,
}

impl<T> Display for SampleFromStateError<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Failed to sample a transition out of state {}", self.exiting)
    }
}

impl<T> Error for SampleFromStateError<T>
where
    T: PhmmNumber + Debug + 'static,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.kind {
            SampleFromStateErrorKind::State(e) => Some(e),
            SampleFromStateErrorKind::StateOrExit(e) => Some(e),
            SampleFromStateErrorKind::EndOrInsert(e) => Some(e),
            SampleFromStateErrorKind::EndInsertOrExit(e) => Some(e),
        }
    }
}

impl<T> GetCode for SampleFromStateError<T> {
    fn get_code(&self) -> i32 {
        match &self.kind {
            SampleFromStateErrorKind::State(e) => e.get_code(),
            SampleFromStateErrorKind::StateOrExit(e) => e.get_code(),
            SampleFromStateErrorKind::EndOrInsert(e) => e.get_code(),
            SampleFromStateErrorKind::EndInsertOrExit(e) => e.get_code(),
        }
    }
}

/// An enum unifying the different sampling errors that could occur within a
/// layer of the core pHMM.
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum SamplingLayerErrorKind<T, const S: usize> {
    /// A sampling error while choosing the next state to enter.
    FromState(SampleFromStateError<T>),
    /// A sampling error while choosing a symbol to emit from a match state.
    EmissionMatch(ParamSamplingError<T, char, S>),
    /// A sampling error while choosing a symbol to emit from an insert
    /// state.
    EmissionInsert(ParamSamplingError<T, char, S>),
}

/// A sampling error that occurred within a layer of the core pHMM.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct LayerSamplingError<T, const S: usize> {
    /// The cause of the sampling error.
    pub kind:  SamplingLayerErrorKind<T, S>,
    /// The layer that the sampling is originating in.
    pub layer: DpIndex,
}

impl<T, const S: usize> Display for LayerSamplingError<T, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Failed to perform sampling in layer: {}", self.layer.0)
    }
}

impl<T, const S: usize> Error for LayerSamplingError<T, S>
where
    T: PhmmNumber + Debug + 'static,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.kind {
            SamplingLayerErrorKind::FromState(e) => Some(e),
            SamplingLayerErrorKind::EmissionMatch(e) | SamplingLayerErrorKind::EmissionInsert(e) => Some(e),
        }
    }
}

impl<T, const S: usize> GetCode for LayerSamplingError<T, S> {
    fn get_code(&self) -> i32 {
        match &self.kind {
            SamplingLayerErrorKind::FromState(e) => e.get_code(),
            SamplingLayerErrorKind::EmissionMatch(e) | SamplingLayerErrorKind::EmissionInsert(e) => e.get_code(),
        }
    }
}

/// A sampling error that occurred when choosing whether to enter the insert
/// state of a [`DomainModule`] or skip to the end.
///
/// [`DomainModule`]: crate::alignment::phmm::modules::DomainModule
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct DomainEnterInsertError<T>(pub ParamSamplingError<T, EndInsert, 2>);

impl<T> Display for DomainEnterInsertError<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Failed to sample whether to enter the insert state or transition to the end of the module"
        )
    }
}

impl<T> Error for DomainEnterInsertError<T>
where
    T: PhmmNumber + Debug + 'static,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.0)
    }
}

impl<T> GetCode for DomainEnterInsertError<T> {
    fn get_code(&self) -> i32 {
        self.0.get_code()
    }
}

/// A sampling error that occurred when choosing whether to remain in the insert
/// state of a [`DomainModule`] or exit to the end.
///
/// [`DomainModule`]: crate::alignment::phmm::modules::DomainModule
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct DomainExitInsertError<T>(pub TransitionParamError<T, EndInsert, 2>);

impl<T> Display for DomainExitInsertError<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Failed to sample whether to remain in the insert state or exit to the end of the module"
        )
    }
}

impl<T> Error for DomainExitInsertError<T>
where
    T: PhmmNumber + Debug + 'static,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.0)
    }
}

impl<T> GetCode for DomainExitInsertError<T> {
    fn get_code(&self) -> i32 {
        self.0.get_code()
    }
}

/// A sampling error that occurred when choosing the layer of the core pHMM to
/// enter from the [`SemiLocalModule`] at the start of the pHMM.
///
/// [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SemiLocalEnterCoreError<T> {
    /// The parameters contained in the [`SemiLocalModule`], describing the
    /// likelihood of entering each layer of the core pHMM.
    ///
    /// This is not displayed in the error message since it can be very large.
    ///
    /// [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule
    pub params: Vec<T>,
    /// The cause of the sampling error.
    pub source: ParamWeightError,
}

impl<T> Display for SemiLocalEnterCoreError<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Failed to sample a layer to enter from the semilocal module at the beginning of the pHMM"
        )
    }
}

impl<T> Error for SemiLocalEnterCoreError<T>
where
    T: Debug,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.source)
    }
}

impl<T> GetCode for SemiLocalEnterCoreError<T> {
    fn get_code(&self) -> i32 {
        self.source.get_code()
    }
}

/// An enum unifying the different sampling errors that could occur within a
/// module at the start or end of a pHMM.
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum ModuleSamplingErrorKind<T, const S: usize> {
    /// A sampling error when choosing whether to enter the insert state of a
    /// [`DomainModule`].
    ///
    /// [`DomainModule`]: crate::alignment::phmm::modules::DomainModule
    DomainEnterInsert(DomainEnterInsertError<T>),
    /// A sampling error when choosing whether to exit the insert state of a
    /// [`DomainModule`].
    ///
    /// [`DomainModule`]: crate::alignment::phmm::modules::DomainModule
    DomainExitInsert(DomainExitInsertError<T>),
    /// A sampling error when choosing the symbol to emit within a
    /// [`DomainModule`].
    ///
    /// [`DomainModule`]: crate::alignment::phmm::modules::DomainModule
    DomainSampleEmissionsError(ParamSamplingError<T, char, S>),
    /// A sampling error when choosing the layer to enter from the
    /// [`SemiLocalModule`] at the beginning of a [`SemiLocalPhmm`].
    ///
    /// [`SemiLocalModule`]: crate::alignment::phmm::modules::SemiLocalModule
    /// [`SemiLocalPhmm`]: crate::alignment::phmm::SemiLocalPhmm
    SemiLocalEnterCore(SemiLocalEnterCoreError<T>),
}

/// A sampling error that occurred within a module at the start or end of a
/// pHMM.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct ModuleSamplingError<T, const S: usize> {
    /// The cause of the sampling error.
    pub kind: ModuleSamplingErrorKind<T, S>,
    /// The location of the module.
    pub loc:  ModuleLocation,
}

impl<T, const S: usize> Display for ModuleSamplingError<T, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let loc = match self.loc {
            ModuleLocation::Begin => "beginning",
            ModuleLocation::End => "end",
        };

        write!(f, "Failed to sample from the module at the {loc} of the pHMM")
    }
}

impl<T, const S: usize> Error for ModuleSamplingError<T, S>
where
    T: PhmmNumber + Debug + 'static,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.kind {
            ModuleSamplingErrorKind::DomainEnterInsert(err) => Some(err),
            ModuleSamplingErrorKind::DomainExitInsert(err) => Some(err),
            ModuleSamplingErrorKind::DomainSampleEmissionsError(err) => Some(err),
            ModuleSamplingErrorKind::SemiLocalEnterCore(err) => Some(err),
        }
    }
}

impl<T, const S: usize> GetCode for ModuleSamplingError<T, S> {
    fn get_code(&self) -> i32 {
        match &self.kind {
            ModuleSamplingErrorKind::DomainEnterInsert(err) => err.get_code(),
            ModuleSamplingErrorKind::DomainExitInsert(err) => err.get_code(),
            ModuleSamplingErrorKind::DomainSampleEmissionsError(err) => err.get_code(),
            ModuleSamplingErrorKind::SemiLocalEnterCore(err) => err.get_code(),
        }
    }
}

/// Any error that can occur while sampling from a pHMM.
///
/// In the error stack, this is a transparent enum adding no additional context.
/// It's display implementation immediately forwards to the contained error
/// variant.
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum SamplingError<T, const S: usize> {
    /// An error caused by an invalid model.
    InvalidModel(InvalidModelError),
    /// An error occurring while sampling within a layer of the core pHMM.
    Layer(LayerSamplingError<T, S>),
    /// An error occurring while sampling within a module at the start or
    /// end of the pHMM.
    Module(ModuleSamplingError<T, S>),
}

impl<T, const S: usize> Display for SamplingError<T, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // This error adds no additional display message... instead immediately
        // defer to inner error
        match self {
            SamplingError::InvalidModel(err) => write!(f, "{err}"),
            SamplingError::Layer(err) => write!(f, "{err}"),
            SamplingError::Module(err) => write!(f, "{err}"),
        }
    }
}

impl<T, const S: usize> Error for SamplingError<T, S>
where
    T: PhmmNumber + Debug + 'static,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        // Since we display the contained err, skip it in the backtrace
        match self {
            SamplingError::InvalidModel(err) => err.source(),
            SamplingError::Layer(err) => err.source(),
            SamplingError::Module(err) => err.source(),
        }
    }
}

impl<T, const S: usize> From<InvalidModelError> for SamplingError<T, S> {
    fn from(value: InvalidModelError) -> Self {
        SamplingError::InvalidModel(value)
    }
}

impl<T, const S: usize> From<LayerSamplingError<T, S>> for SamplingError<T, S> {
    fn from(value: LayerSamplingError<T, S>) -> Self {
        SamplingError::Layer(value)
    }
}

impl<T, const S: usize> From<ModuleSamplingError<T, S>> for SamplingError<T, S> {
    fn from(value: ModuleSamplingError<T, S>) -> Self {
        SamplingError::Module(value)
    }
}
