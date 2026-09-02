//! An API for traversing pHMMs, useful for custom tasks or diagnostic
//! processes.
//!
//! Within *Zoe*, this API supports rescoring and sampling, and it can also be
//! used for testing to examine/log/render the path taken through a pHMM.
//!
//! To use this API, pick or define a "visitor" which makes decisions about
//! which transitions/emissions to make within a pHMM, and is also able to track
//! any state or perform any side-effects needed by the user. Then, call one of
//! the traverse functions, such as [`traverse_global_phmm`],
//! [`traverse_local_phmm`], and so on.
//!
//! Traversing a pHMM is complicated due to complex indexing, edge-cases at the
//! beginning and end, and the potential need to handle modules. The traversal
//! functions abstract away this complexity, allowing the visitors to focus
//! solely on the business logic.
//!
//! ## Implementing Visitors
//!
//! To implement a visitor, implement any of the traits: [`GlobalVisitor`],
//! [`LocalVisitor`], [`SemiLocalVisitor`], and/or [`DomainVisitor`]. If the
//! logic for the visitor is drastically different based on the pHMM type, it
//! might be worth having four distinct visitors, each implementing one trait.
//! Otherwise, a single visitor can implement all of the traits.
//!
//! Each visitor also has a [`finalize`] method, which consumes the visitor and
//! returns a given [`Output`]. This allows the visitor to return a resulting
//! output. If this is not needed (e.g., the visitor relies solely on
//! side-effects), then [`Output`] can be set to `()`. The [`finalize`] method
//! can also perform any error checks or other steps.
//!
//! ## Provided Visitors
//!
//! A few visitors are provided by *Zoe*, which can be composed as part of a
//! custom visitor:
//!
//! - [`GlobalAlignmentVisitor`], [`LocalAlignmentVisitor`],
//!   [`SemiLocalAlignmentVisitor`], and [`DomainAlignmentVisitor`] provide
//!   methods for traversing a pHMM given an alignment.
//! - [`SampleVisitor`] randomly samples a path through the pHMM based on the
//!   parameter weights, returning the formed sequence and alignment.
//! - [`ScoreVisitor`] is a wrapper around an existing visitor that also
//!   computes the score during traversal.
//!
//! [`GlobalAlignmentVisitor`]:
//!     crate::alignment::phmm::traverse::alignment::GlobalAlignmentVisitor
//! [`LocalAlignmentVisitor`]:
//!     crate::alignment::phmm::traverse::alignment::LocalAlignmentVisitor
//! [`SemiLocalAlignmentVisitor`]:
//!     crate::alignment::phmm::traverse::alignment::SemiLocalAlignmentVisitor
//! [`DomainAlignmentVisitor`]:
//!     crate::alignment::phmm::traverse::alignment::DomainAlignmentVisitor
//! [`SampleVisitor`]: crate::alignment::phmm::sampling::SampleVisitor
//! [`ScoreVisitor`]:
//!     crate::alignment::phmm::traverse::score_from_path::ScoreVisitor
//! [`Output`]: GlobalVisitor::Output
//! [`finalize`]: GlobalVisitor::finalize

// TODO: Document ScoreVisitor being compatible with Viterbi

use crate::{
    alignment::phmm::{
        DomainPhmm, GlobalPhmm, LocalPhmm, PhmmNumber, SemiLocalPhmm,
        components::{EmissionParams, TransitionParams},
        indexing::{DpIndex, GetMapping, GetModule, PhmmIndex},
        modules::{DomainModule, LocalModule},
        state::{PhmmState, PhmmStateOrModule},
    },
    data::ByteIndexMap,
};
use std::fmt::Display;

pub mod alignment;
pub mod score_from_path;
mod traversal_fns;
mod visitor;

pub use visitor::*;

/// A trait combining the shared methods between [`GlobalVisitor`] and
/// [`DomainVisitor`] for visiting the core pHMM and requiring the BEGIN and END
/// states to be passed through.
///
/// ## Parameters
///
/// - `V`: The type of the visitor
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
#[allow(clippy::missing_errors_doc, reason = "this trait is a wrapper around other traits")]
trait VisitCore<V, T, const S: usize> {
    /// The type of error issued by the visitor `V`.
    type Error;

    /// See [`GlobalVisitor::choose_emission`] and
    /// [`DomainVisitor::choose_emission`].
    fn choose_emission(
        &self, visitor: &mut V, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
    ) -> Result<usize, Self::Error>;

    /// See [`GlobalVisitor::choose_core_transition`] and
    /// [`DomainVisitor::choose_core_transition`].
    fn choose_core_transition(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<PhmmState, Self::Error>;

    /// See [`GlobalVisitor::choose_end_or_insert`] and
    /// [`DomainVisitor::choose_end_or_insert`].
    fn choose_end_or_insert(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<EndInsert, Self::Error>;
}

impl<V, T, const S: usize> VisitCore<V, T, S> for GlobalPhmm<T, S>
where
    V: GlobalVisitor<T, S>,
{
    type Error = V::Error;

    #[inline]
    fn choose_emission(
        &self, visitor: &mut V, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
    ) -> Result<usize, Self::Error> {
        visitor.choose_emission(layer, state, params, map, self)
    }

    #[inline]
    fn choose_core_transition(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<PhmmState, Self::Error> {
        visitor.choose_core_transition(layer, exiting, params, self)
    }

    #[inline]
    fn choose_end_or_insert(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<EndInsert, Self::Error> {
        visitor.choose_end_or_insert(layer, exiting, params, self)
    }
}

impl<V, T, const S: usize> VisitCore<V, T, S> for DomainPhmm<T, S>
where
    V: DomainVisitor<T, S>,
{
    type Error = V::Error;

    #[inline]
    fn choose_emission(
        &self, visitor: &mut V, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
    ) -> Result<usize, Self::Error> {
        visitor.choose_emission(layer, state, params, map, self)
    }

    #[inline]
    fn choose_core_transition(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<PhmmState, Self::Error> {
        visitor.choose_core_transition(layer, exiting, params, self)
    }

    #[inline]
    fn choose_end_or_insert(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<EndInsert, Self::Error> {
        visitor.choose_end_or_insert(layer, exiting, params, self)
    }
}

/// A trait combining the shared methods between [`SemiLocalVisitor`] and
/// [`LocalVisitor`] for visiting the core pHMM and optionally allowing the
/// BEGIN and END states to be skipped.
///
/// ## Parameters
///
/// - `V`: The type of the visitor
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
#[allow(clippy::missing_errors_doc, reason = "this trait is a wrapper around other traits")]
trait VisitCoreOrExit<V, T, const S: usize> {
    /// The type of error issued by the visitor `V`.
    type Error;

    /// See [`SemiLocalVisitor::choose_emission`] and
    /// [`LocalVisitor::choose_emission`].
    fn choose_emission(
        &self, visitor: &mut V, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
    ) -> Result<usize, Self::Error>;

    /// See [`SemiLocalVisitor::choose_core_transition`] and
    /// [`LocalVisitor::choose_core_transition`].
    fn choose_core_transition(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<PhmmState, Self::Error>;

    /// See [`SemiLocalVisitor::choose_end_or_insert`] and
    /// [`LocalVisitor::choose_end_or_insert`].
    fn choose_end_or_insert(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<EndInsert, Self::Error>;

    /// See [`SemiLocalVisitor::choose_core_transition_or_exit`] and
    /// [`LocalVisitor::choose_core_transition_or_exit`].
    fn choose_core_transition_or_exit(
        &self, visitor: &mut V, layer: DpIndex, params: &TransitionParams<T>, exit_param: T,
    ) -> Result<PhmmStateOrModule, Self::Error>;

    /// See [`SemiLocalVisitor::choose_end_insert_or_exit`] and
    /// [`LocalVisitor::choose_end_insert_or_exit`].
    fn choose_end_insert_or_exit(
        &self, visitor: &mut V, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
    ) -> Result<EndInsertExit, Self::Error>;

    /// See [`SemiLocalVisitor::exit_core_from_end`] and
    /// [`LocalVisitor::exit_core_from_end`].
    fn exit_core_from_end(&self, visitor: &mut V, layer: DpIndex, exit_param: T) -> Result<(), Self::Error>;
}

impl<V, T, const S: usize> VisitCoreOrExit<V, T, S> for SemiLocalPhmm<T, S>
where
    V: SemiLocalVisitor<T, S>,
{
    type Error = V::Error;

    #[inline]
    fn choose_emission(
        &self, visitor: &mut V, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
    ) -> Result<usize, Self::Error> {
        visitor.choose_emission(layer, state, params, map, self)
    }

    #[inline]
    fn choose_core_transition(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<PhmmState, Self::Error> {
        visitor.choose_core_transition(layer, exiting, params, self)
    }

    #[inline]
    fn choose_end_or_insert(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<EndInsert, Self::Error> {
        visitor.choose_end_or_insert(layer, exiting, params, self)
    }

    #[inline]
    fn choose_core_transition_or_exit(
        &self, visitor: &mut V, layer: DpIndex, params: &TransitionParams<T>, exit_param: T,
    ) -> Result<PhmmStateOrModule, Self::Error> {
        visitor.choose_core_transition_or_exit(layer, params, exit_param, self)
    }

    #[inline]
    fn choose_end_insert_or_exit(
        &self, visitor: &mut V, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
    ) -> Result<EndInsertExit, Self::Error> {
        visitor.choose_end_insert_or_exit(layer, params, exit_param, exit_from_end_param, self)
    }

    #[inline]
    fn exit_core_from_end(&self, visitor: &mut V, layer: DpIndex, exit_param: T) -> Result<(), Self::Error> {
        visitor.exit_core_from_end(layer, exit_param, self)
    }
}

impl<V, T, const S: usize> VisitCoreOrExit<V, T, S> for LocalPhmm<T, S>
where
    V: LocalVisitor<T, S>,
{
    type Error = V::Error;

    #[inline]
    fn choose_emission(
        &self, visitor: &mut V, layer: DpIndex, state: PhmmState, params: &EmissionParams<T, S>, map: &ByteIndexMap<S>,
    ) -> Result<usize, Self::Error> {
        visitor.choose_emission(layer, state, params, map, self)
    }

    #[inline]
    fn choose_core_transition(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<PhmmState, Self::Error> {
        visitor.choose_core_transition(layer, exiting, params, self)
    }

    #[inline]
    fn choose_end_or_insert(
        &self, visitor: &mut V, layer: DpIndex, exiting: PhmmState, params: &TransitionParams<T>,
    ) -> Result<EndInsert, Self::Error> {
        visitor.choose_end_or_insert(layer, exiting, params, self)
    }

    #[inline]
    fn choose_core_transition_or_exit(
        &self, visitor: &mut V, layer: DpIndex, params: &TransitionParams<T>, exit_param: T,
    ) -> Result<PhmmStateOrModule, Self::Error> {
        visitor.choose_core_transition_or_exit(layer, params, exit_param, self)
    }

    #[inline]
    fn choose_end_insert_or_exit(
        &self, visitor: &mut V, layer: DpIndex, params: &TransitionParams<T>, exit_param: T, exit_from_end_param: T,
    ) -> Result<EndInsertExit, Self::Error> {
        visitor.choose_end_insert_or_exit(layer, params, exit_param, exit_from_end_param, self)
    }

    #[inline]
    fn exit_core_from_end(&self, visitor: &mut V, layer: DpIndex, exit_param: T) -> Result<(), Self::Error> {
        visitor.exit_core_from_end(layer, exit_param, self)
    }
}

/// A trait combining the shared methods between [`DomainVisitor`] and
/// [`LocalVisitor`] for visiting each [`DomainModule`] in the pHMM.
///
/// ## Parameters
///
/// - `V`: The type of the visitor
/// - `T`: The type of the parameters in the pHMM
/// - `S`: The alphabet size of the pHMM
#[allow(clippy::missing_errors_doc, reason = "this trait is a wrapper around other traits")]
trait VisitDomainModule<V, T, const S: usize> {
    /// The type of error issued by the visitor `V`.
    type Error;

    /// See [`DomainVisitor::choose_domain_emission`] and
    /// [`LocalVisitor::choose_local_emission`].
    fn choose_domain_emission(
        &self, visitor: &mut V, params: &EmissionParams<T, S>, mapping: &ByteIndexMap<S>, loc: ModuleLocation,
    ) -> Result<usize, Self::Error>;

    /// See [`DomainVisitor::enter_module_insert`] and
    /// [`LocalVisitor::enter_module_insert`].
    fn enter_module_insert(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<bool, Self::Error>;

    /// See [`DomainVisitor::exit_module_insert`] and
    /// [`LocalVisitor::exit_module_insert`].
    fn exit_module_insert(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<bool, Self::Error>;

    /// See [`DomainVisitor::exiting_module`] and
    /// [`LocalVisitor::exiting_domain_module`].
    fn exiting_domain_module(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<(), Self::Error>;
}

impl<V, T, const S: usize> VisitDomainModule<V, T, S> for DomainPhmm<T, S>
where
    V: DomainVisitor<T, S>,
{
    type Error = V::Error;

    #[inline]
    fn choose_domain_emission(
        &self, visitor: &mut V, params: &EmissionParams<T, S>, mapping: &ByteIndexMap<S>, loc: ModuleLocation,
    ) -> Result<usize, Self::Error> {
        visitor.choose_domain_emission(params, mapping, loc, self)
    }

    #[inline]
    fn enter_module_insert(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<bool, Self::Error> {
        visitor.enter_module_insert(module, loc, self)
    }

    #[inline]
    fn exit_module_insert(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<bool, Self::Error> {
        visitor.exit_module_insert(module, loc, self)
    }

    #[inline]
    fn exiting_domain_module(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<(), Self::Error> {
        visitor.exiting_module(module, loc, self)
    }
}

impl<V, T, const S: usize> VisitDomainModule<V, T, S> for LocalPhmm<T, S>
where
    V: LocalVisitor<T, S>,
{
    type Error = V::Error;

    #[inline]
    fn choose_domain_emission(
        &self, visitor: &mut V, params: &EmissionParams<T, S>, mapping: &ByteIndexMap<S>, loc: ModuleLocation,
    ) -> Result<usize, Self::Error> {
        visitor.choose_local_emission(params, mapping, loc, self)
    }

    #[inline]
    fn enter_module_insert(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<bool, Self::Error> {
        visitor.enter_module_insert(module, loc, self)
    }

    #[inline]
    fn exit_module_insert(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<bool, Self::Error> {
        visitor.exit_module_insert(module, loc, self)
    }

    #[inline]
    fn exiting_domain_module(
        &self, visitor: &mut V, module: &DomainModule<T, S>, loc: ModuleLocation,
    ) -> Result<(), Self::Error> {
        visitor.exiting_domain_module(module, loc, self)
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum EndInsert {
    End,
    Insert,
}

impl Display for EndInsert {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EndInsert::End => write!(f, "End"),
            EndInsert::Insert => write!(f, "Insert"),
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum EndInsertExit {
    End,
    Insert,
    Exit,
}

impl Display for EndInsertExit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EndInsertExit::End => write!(f, "End"),
            EndInsertExit::Insert => write!(f, "Insert"),
            EndInsertExit::Exit => write!(f, "Exit"),
        }
    }
}

impl<T: PhmmNumber, const S: usize> DomainModule<T, S> {
    /// Lazily compute the score for skipping `inserted` residues at the
    /// beginning of the query. This should only be used for diagnostics or
    /// testing, otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    fn get_begin_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        if inserted.is_empty() {
            self.start_to_end
        } else {
            // Special casing needed in case insert_to_insert is infinite,
            // causing a NAN to appear when multiplied by 0
            let insert_to_insert = if inserted.len() > 1 {
                T::cast_from(inserted.len() - 1) * self.insert_to_insert
            } else {
                T::ZERO
            };

            self.start_to_insert
                + self.insert_to_end
                + (inserted
                    .iter()
                    .map(|x| self.background_emission[mapping.to_index(*x)])
                    .fold(T::ZERO, |acc, elem| acc + elem)
                    + insert_to_insert)
        }
    }

    /// Lazily compute the score for skipping `inserted` residues at the end of
    /// the query. This should only be used for diagnostics or testing,
    /// otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    fn get_end_score(&self, inserted: &[u8], mapping: &'static ByteIndexMap<S>) -> T {
        if inserted.is_empty() {
            self.start_to_end
        } else {
            // Special casing needed in case insert_to_insert is infinite,
            // causing a NAN to appear when multiplied by 0
            let insert_to_insert = if inserted.len() > 1 {
                T::cast_from(inserted.len() - 1) * self.insert_to_insert
            } else {
                T::ZERO
            };

            self.start_to_insert
                + self.insert_to_end
                + (inserted
                    .iter()
                    .rev()
                    .map(|x| self.background_emission[mapping.to_index(*x)])
                    .fold(T::ZERO, |acc, elem| acc + elem)
                    + insert_to_insert)
        }
    }
}

trait GetScoreDomain<T, const S: usize>:
    GetModule<Begin = DomainModule<T, S>, End = DomainModule<T, S>> + GetMapping<S>
where
    T: PhmmNumber, {
    /// Lazily compute the score for skipping `inserted` residues at the
    /// beginning of the query. This should only be used for diagnostics or
    /// testing, otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    #[inline]
    fn get_begin_score(&self, inserted: &[u8]) -> T {
        self.begin().get_begin_score(inserted, self.mapping())
    }

    /// Lazily compute the score for skipping `inserted` residues at the end of
    /// the query. This should only be used for diagnostics or testing,
    /// otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    #[inline]
    fn get_end_score(&self, inserted: &[u8]) -> T {
        self.end().get_end_score(inserted, self.mapping())
    }
}

impl<P, T, const S: usize> GetScoreDomain<T, S> for P
where
    T: PhmmNumber,
    P: GetModule<Begin = DomainModule<T, S>, End = DomainModule<T, S>> + GetMapping<S>,
{
}

trait GetScoreLocal<T, const S: usize>: GetModule<Begin = LocalModule<T, S>, End = LocalModule<T, S>> + GetMapping<S>
where
    T: PhmmNumber, {
    /// Lazily compute the score for skipping `inserted` residues at the
    /// beginning of the query. This should only be used for diagnostics or
    /// testing, otherwise [`PrecomputedDomainModule`] should be used.
    ///
    /// This is designed to give the exact same score as the precomputed
    /// version, performing all arithmetic operations in the same order so as
    /// not to change the floating point error.
    ///
    /// [`PrecomputedDomainModule`]:
    ///     crate::alignment::phmm::modules::PrecomputedDomainModule
    #[inline]
    fn get_begin_domain_score(&self, inserted: &[u8]) -> T {
        self.begin().domain_params.get_begin_score(inserted, self.mapping())
    }

    #[inline]
    fn get_begin_semilocal_score(&self, index: impl PhmmIndex) -> T {
        self.begin().semilocal_params.get_score(index)
    }

    #[inline]
    fn get_end_domain_score(&self, inserted: &[u8]) -> T {
        self.end().domain_params.get_end_score(inserted, self.mapping())
    }

    #[inline]
    fn get_end_semilocal_score(&self, index: impl PhmmIndex) -> T {
        self.end().semilocal_params.get_score(index)
    }

    #[inline]
    fn get_begin_score(&self, inserted: &[u8], index: impl PhmmIndex) -> T {
        self.get_begin_domain_score(inserted) + self.get_begin_semilocal_score(index)
    }

    #[inline]
    fn get_end_score(&self, inserted: &[u8], index: impl PhmmIndex) -> T {
        self.get_end_domain_score(inserted) + self.get_end_semilocal_score(index)
    }
}

impl<P, T, const S: usize> GetScoreLocal<T, S> for P
where
    T: PhmmNumber,
    P: GetModule<Begin = LocalModule<T, S>, End = LocalModule<T, S>> + GetMapping<S>,
{
}

/// The location of a module in a pHMM.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum ModuleLocation {
    /// The module is placed at the beginning of the pHMM.
    Begin,
    /// The module is placed at the end of the pHMM.
    End,
}
