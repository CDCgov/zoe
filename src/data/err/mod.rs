//! Error types and convenience traits for handling [`Result`].
//!
//! This module provides:
//!
//! - The error type [`ErrorWithContext`], along with the traits
//!   [`ResultWithErrorContext`] and [`WithErrorContext`], for wrapping errors
//!   with additional context while preserving the error source chain.
//! - [`OrFail`] and [`Fail`] for graceful CLI error handling with exit codes.
//!
//! ## Error Handling Philosophy in *Zoe*
//!
//! As a library, *Zoe* aims to avoid making assumptions on the style of error
//! handling chosen by users, in particular by not adopting any error handling
//! crate as a dependency.
//!
//! For specific applications, *Zoe* has enum-style error types such as
//! [`ProfileError`] or [`KmerError`], which the user can match on or display.
//! For working with files and record types, however, *Zoe* elects to use
//! [`std::io::Error`], allowing for system IO errors to be propagated and
//! function-specific error messages to be represented with
//! [`ErrorKind::InvalidData`] or [`ErrorKind::Other`].
//!
//! Similar to [`std::io::Error`], [`std::fmt::Display`] is implemented only at
//! the immediate error level. To see the full error stack when handling errors,
//! it is important to do one of the following:
//!
//! - Return a `Result<(), ErrorWithContext>` from `main.rs`
//! - Use *Zoe*'s [`OrFail`] or [`Fail`] traits
//! - Iterate through the source chain via [`Error::source`]
//! - Use an external crate like `anyhow`
//!
//! ## Error Context
//!
//! *Zoe* elects to add context to error messages by default when available,
//! such as in [`FastQReader::from_path`] which will include the path of a
//! missing/empty file. In bioinformatics, this context is very useful in
//! complex pipelines, and any runtime penalty is considered negligible compared
//! to the algorithms being run.
//!
//! This context is added using the [`WithErrorContext`] and
//! [`ResultWithErrorContext`] traits. They add context by creating a
//! [`ErrorWithContext`] struct, containing the original error (boxed) as the
//! [`Error::source`] and the context as the new top-level error, stored as a
//! [`String`].
//!
//! For flexibility, there are also the following ways to construct an
//! [`ErrorWithContext`]:
//!
//! - Conversion from [`std::io::Error`] without any additional context, using
//!   [`From::from`].
//! - Conversion from any error without any additional context, using
//!   [`WrapError`] or [`WrapResult`].
//! - Initialization with just a line of context and no source, using
//!   [`ErrorWithContext::new`]. This is conceptually similar to
//!   [`std::io::Error::other`].
//!
//! [`ErrorWithContext`]: crate::data::err::ErrorWithContext
//! [`ResultWithErrorContext`]: crate::data::err::ResultWithErrorContext
//! [`WithErrorContext`]: crate::data::err::WithErrorContext
//! [`OrFail`]: crate::data::err::OrFail
//! [`FastQReader::from_path`]: crate::prelude::FastQReader::from_path
//! [`Error::source`]: std::error::Error::source
//! [`ProfileError`]: crate::alignment::ProfileError
//! [`KmerError`]: crate::kmer::KmerError
//! [`ErrorKind::InvalidData`]: std::io::ErrorKind::InvalidData
//! [`ErrorKind::Other`]: std::io::ErrorKind::Other

mod display;
mod error_with_context;
mod handling;

pub use display::*;
pub use error_with_context::*;
pub use handling::*;

#[cfg(test)]
mod test;

/// Maximum number of errors inspected when traversing a source chain, both for
/// retrieving exit codes and for displaying the stack.
const MAX_ERROR_CHAIN_DEPTH: usize = 256;
