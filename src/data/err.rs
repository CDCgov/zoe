//! Error types and convenience traits for handling [`Result`].
//!
//! This module provides:
//!
//! - The error type [`ErrorWithContext`], along with the traits
//!   [`ResultWithErrorContext`] and [`WithErrorContext`], for wrapping errors
//!   with additional context while preserving the error source chain.
//! - [`GetCode`], [`OrFail`], and [`Fail`] for graceful CLI error handling with
//!   exit codes.
//!
//! ## Error Handling Philosophy
//!
//! *Zoe* often wraps [`std::io::Error`] for I/O and parsing operations rather
//! than defining new error types. Domain-specific errors (e.g., alignment,
//! distance) use enums in their respective modules. [`std::fmt::Display`] is
//! implemented only at the immediate error level—callers and loggers must
//! iterate through the source chain via [`std::error::Error::source`], use
//! *Zoe*'s [`OrFail`] or [`Fail`], or an external crates like `anyhow`.
//!
//! *Zoe* elects to add context by default when available, such as in
//! [`FastQReader::from_path`] which will include the path of a missing/empty
//! file. In bioinformatics, this context is very useful in complex pipelines,
//! and any runtime penalty is considered negligible compared to the algorithms
//! being run.
//!
//! This context is added using the [`WithErrorContext`] and
//! [`ResultWithErrorContext`] traits. They add context by creating a
//! [`ErrorWithContext`] struct, containing the original error (boxed) as the
//! [`Error::source`] and the context as the new top-level error, stored as a
//! [`String`]. When used in conjuction with an error handling library such as
//! `anyhow` or *Zoe*'s [`OrFail`] or [`Fail`] traits, the stack of errors can
//! be displayed in an application.
//!
//! [`ErrorWithContext`] can also be constructed directly without a source
//! error. In applications that are avoiding dependencies such as `anyhow` and
//! do not want to use [`std::io::Error::other`], [`ErrorWithContext::new`] is a
//! viable option.
//!
//! [`ErrorWithContext`]: crate::data::err::ErrorWithContext
//! [`ResultWithErrorContext`]: crate::data::err::ResultWithErrorContext
//! [`WithErrorContext`]: crate::data::err::WithErrorContext
//! [`GetCode`]: crate::data::err::GetCode
//! [`OrFail`]: crate::data::err::OrFail
//! [`FastQReader::from_path`]: crate::prelude::FastQReader::from_path
//! [`Error::source`]: std::error::Error::source

use std::{
    error::Error,
    fmt::{Debug, Display, Write},
    hint::cold_path,
    path::Path,
};

/// A macro for unwrapping a [`Result`] and propagating any error as a
/// `Some(Err(e))`.
///
/// This is especially useful for fallible iterators, where results need to be
/// wrapped in [`Some`].
#[macro_export]
macro_rules! unwrap_or_return_some_err {
    ($expression:expr) => {
        match $expression {
            Ok(v) => v,
            Err(e) => return Some(Err(e)),
        }
    };
}

/// Trait for specifying getting exit codes originating from IO errors.
///
/// Implementing this trait allows the error type to work with [`OrFail`]. If
/// this is being implemented for a top-level error containing a nested error,
/// one should manually implement [`get_code`] to retrieve the underlying code
/// for the nested error, rather than using the blanket implementation.
///
/// [`get_code`]: GetCode::get_code
pub trait GetCode {
    /// Retrieves the exit code associated with a given error.
    ///
    /// ## Validity
    ///
    /// If this method is manually implemented, then it must recursively call
    /// [`get_code`] on [`Error::source`]. Any other behavior or logic is not
    /// guaranteed to be consistent within *Zoe*, and may or may not be
    /// correctly applied when using wrapped errors ([`ErrorWithContext`]).
    ///
    /// The blanket implementation returns `1`. We also implement on
    /// [`std::io::Error`] to return the [`raw_os_error`] if available.
    ///
    /// [`raw_os_error`]: std::io::Error::raw_os_error
    /// [`get_code`]: GetCode::get_code
    #[inline]
    #[must_use]
    fn get_code(&self) -> i32 {
        1
    }
}

impl GetCode for std::io::Error {
    #[inline]
    fn get_code(&self) -> i32 {
        if let Some(code) = self.raw_os_error() {
            return code;
        }

        let mut source = self.source();
        while let Some(err) = source {
            if let Some(e) = err.downcast_ref::<std::io::Error>()
                && let Some(code) = e.raw_os_error()
            {
                return code;
            }
            source = err.source();
        }

        1
    }
}

/// A trait for providing more graceful error reporting and aborting.
///
/// A status code is provided by [`GetCode`], and any context available in
/// [`Error::source`] is displayed.
///
/// <div class="warning note">
///
/// **Note**
///
/// To get full utility out of this trait, custom top-level errors should
/// manually implement [`std::error::Error::source`] and get whatever field or
/// variants contains the nested errors. In addition, [`GetCode`] should
/// likewise be implemented manually to retrieve the underlying codes for nested
/// [`std::io::Error`].
///
/// </div>
pub trait OrFail<T> {
    /// Unwraps the result, writing the error and any information in
    /// [`Error::source`] to stderr.
    fn unwrap_or_fail(self) -> T;

    /// Unwraps the result, writing the provided message, the error, and any
    /// information in [`Error::source`] to stderr.
    fn unwrap_or_die(self, msg: &str) -> T;
}

impl<T, E> OrFail<T> for Result<T, E>
where
    E: GetCode + Display + Error + 'static,
{
    fn unwrap_or_fail(self) -> T {
        match self {
            Ok(result) => result,
            Err(e) => {
                cold_path();
                e.fail()
            }
        }
    }

    fn unwrap_or_die(self, msg: &str) -> T {
        match self {
            Ok(result) => result,
            Err(e) => {
                cold_path();
                e.die(msg)
            }
        }
    }
}

/// A trait for providing more graceful error reporting and aborting. For
/// similar methods on [`Result`], see [`OrFail`].
///
/// A status code is provided by [`GetCode`], and any context available in
/// [`Error::source`] is displayed.
///
/// <div class="warning note">
///
/// **Note**
///
/// To get full utility out of this trait, custom top-level errors should
/// manually implement [`std::error::Error::source`] and get whatever field or
/// variants contains the nested errors. In addition, [`GetCode`] should
/// likewise be implemented manually to retrieve the underlying codes for nested
/// [`std::io::Error`].
///
/// </div>
pub trait Fail {
    /// Exits the program, writing the error and any information in
    /// [`Error::source`] to stderr.
    fn fail(self) -> !;

    /// Exits the program, writing the provided message, the error, and any
    /// information in [`Error::source`] to stderr.
    fn die(self, msg: &str) -> !;
}

impl<E> Fail for E
where
    E: GetCode + Display + Error + 'static,
{
    fn fail(self) -> ! {
        if let Ok(bin) = std::env::current_exe() {
            eprintln!("Error in {b}", b = bin.display());
        } else {
            eprintln!("Error in program");
        }

        print_stack(&self);
        std::process::exit(self.get_code());
    }

    fn die(self, msg: &str) -> ! {
        if let Ok(bin) = std::env::current_exe() {
            eprintln!("Error in {b}: {msg}", b = bin.display());
        } else {
            eprintln!("Error: {msg}");
        }

        print_stack(&self);
        std::process::exit(self.get_code());
    }
}

/// An error type supporting context and a backtrace.
///
/// Specifically, this error can hold up to three things:
///
/// 1. An optional source error message, which this error wraps. Using
///    [`unwrap_or_fail`] or [`unwrap_or_die`] cause the source error to be
///    shown in the backtrace. This source is accessible via [`Error::source`].
/// 2. A line of context describing the error. This appears as one item in the
///    [`OrFail`] backtrace.
/// 3. Any subitems (additional indented lines with more information that appear
///    below the line of context). This is useful for including the values of
///    variables or other useful information.
///
/// This can be converted to [`std::io::Error`] with [`Into`]. Hence, in
/// functions returning [`std::io::Result`], the `?` operator can be used after
/// adding context.
///
/// [`with_subitem`]: WithSubitem::with_subitem
/// [`unwrap_or_fail`]: OrFail::unwrap_or_fail
/// [`unwrap_or_die`]: OrFail::unwrap_or_die
#[must_use]
#[derive(Debug)]
pub struct ErrorWithContext {
    /// The inner representation. Using a single fat pointer is better than
    /// storing the description and source fields directly, since it minimizes
    /// the size of [`ErrorWithContext`] and hence the size of `Result<T,
    /// ErrorWithContext>`.
    repr: Box<ErrorWithContextRepr>,
}

impl ErrorWithContext {
    /// Constructs a new [`ErrorWithContext`] with the given description,
    /// without a source error or any subitems.
    ///
    /// The `description` may be anything implementing `Into<String>`. Passing
    /// an owned `String` avoids an extra allocation.
    pub fn new(description: impl Into<String>) -> Self {
        ErrorWithContext {
            repr: Box::new(ErrorWithContextRepr {
                description: description.into(),
                subitem:     None,
                source:      None,
            }),
        }
    }
}

/// The inner representation for an [`ErrorWithContext`]. This is wrapped in a
/// [`Box`] in [`ErrorWithContext`] to reduce the memory of `Result<T,
/// ErrorWithContext>` in the `Ok` case.
#[derive(Debug)]
struct ErrorWithContextRepr {
    /// The context that was added to the error.
    description: String,

    /// Any subitems attached to the error, separated by new lines.
    subitem: Option<String>,

    /// The source error which the context is added to.
    source: Option<Box<dyn Error + Send + Sync>>,
}

impl Display for ErrorWithContextRepr {
    /// Displays the description of the error as well as any subitems.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.description)?;
        if let Some(subitem) = &self.subitem {
            write!(
                f,
                "\n| {}",
                IndentWrapper {
                    val:    subitem,
                    indent: "| ",
                }
            )?;
        }
        Ok(())
    }
}

impl std::fmt::Display for ErrorWithContext {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.repr)
    }
}

impl Error for ErrorWithContext {
    #[inline]
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.repr.source {
            Some(source) => Some(source.as_ref()),
            None => None,
        }
    }
}

impl From<ErrorWithContext> for std::io::Error {
    #[inline]
    fn from(e: ErrorWithContext) -> Self {
        std::io::Error::other(e)
    }
}

/// An extension trait for [`Error`] allowing additional context to be added via
/// a [`ErrorWithContext`].
pub trait WithErrorContext {
    /// Wraps the error in an [`ErrorWithContext`] with the given description.
    ///
    /// The `description` may be anything implementing `Into<String>`. Passing
    /// an owned `String` avoids an extra allocation.
    fn with_context(self, description: impl Into<String>) -> ErrorWithContext;

    /// Wraps the error in an [`ErrorWithContext`] by adding type context.
    fn with_type_context<T>(self) -> ErrorWithContext;

    /// Wraps the error in an [`ErrorWithContext`] by adding path context.
    ///
    /// The context will be formatted as `msg: 'path'`. The `msg` may be
    /// anything implementing [`Display`].
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> ErrorWithContext;

    /// Deprecated, use [`with_path_context`] instead.
    ///
    /// [`with_path_context`]: WithErrorContext::with_path_context
    #[deprecated(
        since = "0.0.27",
        note = "please use `with_path_context` instead. This function will be removed in v0.0.29"
    )]
    fn with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> ErrorWithContext;
}

impl<E: Error + Send + Sync + 'static> WithErrorContext for E {
    // Do not inline, since this is cold code
    fn with_context(self, description: impl Into<String>) -> ErrorWithContext {
        ErrorWithContext {
            repr: Box::new(ErrorWithContextRepr {
                description: description.into(),
                subitem:     None,
                source:      Some(Box::new(self)),
            }),
        }
    }

    // Do not inline, since this is cold code
    fn with_type_context<T>(self) -> ErrorWithContext {
        let name = std::any::type_name::<T>();
        let description = format!(
            "Failure in {}",
            name.split('<').next().unwrap_or(name).rsplit("::").next().unwrap_or(name)
        );

        ErrorWithContext {
            repr: Box::new(ErrorWithContextRepr {
                description,
                subitem: None,
                source: Some(Box::new(self)),
            }),
        }
    }

    // Do not inline, since this is cold code
    fn with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> ErrorWithContext {
        self.with_path_context(msg, file)
    }

    // Do not inline, since this is cold code
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> ErrorWithContext {
        Self::with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}

/// An extension trait for [`ErrorWithContext`] allowing an indented subitem to
/// be added to the error (without adding a new error to the backtrace).
pub trait WithSubitem {
    /// Adds a subitem with the given `message` to the error without adding a
    /// new error to the backtrace.
    ///
    /// [`WithErrorContext::with_context`] creates a new entry in the backtrace
    /// (displayed using `→`), whereas this method adds a message that is
    /// indented beneath the error and indicated using `|`. For example:
    ///
    /// ```text
    /// Error in /path/to/binary
    ///   → Failed to load reads from file: input.fastq
    ///   → Failed to deinterleave records due to mismatching headers
    ///     | Header 1: SIM:1:FCX:1:15:6329:1045 1:N:0:2
    ///     | Header 2: SIM:1:FCX:1:15:2345:1001 2:N:0:2
    ///   → x_pos fields did not agree!
    /// ```
    fn with_subitem(self, message: impl Into<String>) -> ErrorWithContext;
}

impl WithSubitem for ErrorWithContext {
    fn with_subitem(mut self, message: impl Into<String>) -> ErrorWithContext {
        let subitem = &mut self.repr.subitem;
        let message = message.into();
        if let Some(subitem) = subitem {
            subitem.push('\n');
            subitem.push_str(&message);
        } else {
            *subitem = Some(message);
        }
        self
    }
}

/// An extension trait for [`Result`] allowing additional context to be added to
/// an [`Err`] variant via a [`ErrorWithContext`].
///
/// The methods are similar to [`WithErrorContext`], but are implemented for
/// results.
pub trait ResultWithErrorContext {
    /// The type of the [`Ok`] variant in the result.
    type Ok;

    /// Wraps the [`Err`] variant in an [`ErrorWithContext`] with the given
    /// description.
    ///
    /// The `description` may be anything implementing `Into<String>`. Passing
    /// an owned `String` avoids an extra allocation.
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_context(self, description: impl Into<String>) -> Result<Self::Ok, ErrorWithContext>;

    /// Wraps the [`Err`] variant in an [`ErrorWithContext`] by adding type
    /// context.
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_type_context<T>(self) -> Result<Self::Ok, ErrorWithContext>;

    /// Wraps the [`Err`] variant in an [`ErrorWithContext`] by adding path
    /// context.
    ///
    /// The context will be formatted as `msg: 'path'`. The `msg` may be
    /// anything implementing [`Display`].
    ///
    /// ## Errors
    ///
    /// Propagates errors in `self`, with the added context.
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> Result<Self::Ok, ErrorWithContext>;

    /// Deprecated, use [`with_path_context`] instead.
    ///
    /// [`with_path_context`]: ResultWithErrorContext::with_path_context
    #[allow(clippy::missing_errors_doc)]
    #[deprecated(
        since = "0.0.27",
        note = "please use `with_path_context` instead. This function will be removed in v0.0.29"
    )]
    fn with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> Result<Self::Ok, ErrorWithContext>;
}

impl<Ok, E: WithErrorContext> ResultWithErrorContext for Result<Ok, E> {
    type Ok = Ok;

    #[inline]
    fn with_context(self, description: impl Into<String>) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| {
            cold_path();
            e.with_context(description)
        })
    }

    #[inline]
    fn with_type_context<T>(self) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| {
            cold_path();
            e.with_type_context::<T>()
        })
    }

    #[inline]
    fn with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> Result<Ok, ErrorWithContext> {
        self.with_path_context(msg, file)
    }

    #[inline]
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> Result<Ok, ErrorWithContext> {
        self.map_err(|e| {
            cold_path();
            e.with_path_context(msg, file)
        })
    }
}

/// An extension trait for `Result<T, WithErrorContext>` allowing information to
/// be attached to an [`Err`] variant.
///
/// The methods are similar to [`WithSubitem`], but are implemented for results.
pub trait ResultWithSubitem {
    /// Adds a subitem with the given `message` to an [`Err`] variant without
    /// adding a new error to the backtrace.
    ///
    /// [`ResultWithErrorContext::with_context`] creates a new entry in the
    /// backtrace (displayed using `→`), whereas this method adds a message that
    /// is indented beneath the error and indicated using `|`. For example:
    ///
    /// ```text
    /// Error in /path/to/binary
    ///   → Failed to load reads from file: input.fastq
    ///   → Failed to deinterleave records due to mismatching headers
    ///     | Header 1: SIM:1:FCX:1:15:6329:1045 1:N:0:2
    ///     | Header 2: SIM:1:FCX:1:15:2345:1001 2:N:0:2
    ///   → x_pos fields did not agree!
    /// ```
    #[must_use]
    fn with_subitem(self, message: impl Into<String>) -> Self;
}

impl<T> ResultWithSubitem for Result<T, ErrorWithContext> {
    #[inline]
    fn with_subitem(self, message: impl Into<String>) -> Self {
        self.map_err(|e| {
            cold_path();
            e.with_subitem(message)
        })
    }
}

/// A wrapper around [`std::fmt::Formatter`] which automatically indents all new
/// lines with a specified string.
///
/// This is a helper struct for the formatting used in [`fail`] and [`die`].
///
/// [`fail`]: Fail::fail
/// [`die`]: Fail::die
struct IndentFormatter<'a, 'b> {
    formatter: &'a mut std::fmt::Formatter<'b>,
    indent:    &'static str,
}

impl Write for IndentFormatter<'_, '_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        let mut parts = s.split('\n');
        let Some(first_part) = parts.next() else { return Ok(()) };
        self.formatter.write_str(first_part)?;

        for part in parts {
            self.formatter.write_char('\n')?;
            self.formatter.write_str(self.indent)?;
            self.formatter.write_str(part)?;
        }

        Ok(())
    }

    fn write_char(&mut self, c: char) -> std::fmt::Result {
        if c == '\n' {
            self.formatter.write_char('\n')?;
            self.formatter.write_str(self.indent)
        } else {
            self.formatter.write_char(c)
        }
    }
}

/// A wrapper type altering the implementation of [`Display`], such that any new
/// lines are automatically indented with a specified string.
///
/// This is a helper struct for the formatting used in [`fail`] and [`die`].
///
/// [`fail`]: Fail::fail
/// [`die`]: Fail::die
struct IndentWrapper<T> {
    /// The value to display.
    val:    T,
    /// The string to use when indenting lines after the first.
    indent: &'static str,
}

impl<T: Display> Display for IndentWrapper<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            IndentFormatter {
                formatter: f,
                indent:    self.indent,
            },
            "{}",
            self.val
        )
    }
}

/// Prints an error and its backtrace using [`Error::source`].
///
/// This encapsulates the shared logic between [`unwrap_or_die`] and
/// [`unwrap_or_fail`]. Dynamic errors are used to prevent monomorphization on
/// the cold path.
///
/// [`unwrap_or_die`]: OrFail::unwrap_or_die
/// [`unwrap_or_fail`]: OrFail::unwrap_or_fail
#[cold]
fn print_stack(err: &(dyn Error + 'static)) {
    // Wrap the error in Some so that we don't have to write the same logic
    // twice
    let mut maybe_err = Some(err);

    while let Some(err) = maybe_err {
        eprintln!(
            "  → {err}",
            err = IndentWrapper {
                val:    err,
                indent: "    ",
            }
        );

        maybe_err = err.source();
    }
}
