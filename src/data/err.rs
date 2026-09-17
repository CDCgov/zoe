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
//! [`ErrorWithContext`] can also be constructed directly without a source
//! error. In applications that are avoiding dependencies such as `anyhow` and
//! do not want to use [`std::io::Error::other`], [`ErrorWithContext::new`] is a
//! viable option.
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

use std::{
    error::Error,
    fmt::{Debug, Display, Write},
    hint::cold_path,
    path::Path,
};

/// Maximum number of errors inspected when traversing a source chain, both for
/// retrieving exit codes and for displaying the stack.
const MAX_ERROR_CHAIN_DEPTH: usize = 256;

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

/// Finds the first raw OS error code in an error's chain.
///
/// Unlike [`Error::source`], this descends into transparent sources and IO
/// payloads. These can hold a raw OS error code which is not visible by the
/// wrapping error.
fn io_exit_code(error: &(dyn Error + 'static)) -> i32 {
    let mut current = Some(error);

    for _ in 0..MAX_ERROR_CHAIN_DEPTH {
        let Some(error) = current else {
            break;
        };

        // Extract the raw_os_error if present
        if let Some(io_error) = error.downcast_ref::<std::io::Error>()
            && let Some(code) = io_error.raw_os_error()
        {
            return code;
        }

        // Get any transparent errors, IO payloads, or the source error
        current = if let Some(context) = error.downcast_ref::<ErrorWithContext>() {
            context.repr.source.as_ref().map(ErrorSource::as_error)
        } else if let Some(payload) = error.downcast_ref::<std::io::Error>().and_then(std::io::Error::get_ref) {
            Some(payload)
        } else {
            error.source()
        };
    }

    1
}

/// Writes an error stack to stderr and exits with `code`.
fn exit_with_error(error: &(dyn Error + 'static), msg: Option<&str>, code: i32) -> ! {
    match (std::env::current_exe(), msg) {
        (Ok(bin), Some(msg)) => eprintln!("Error in {bin}: {msg}", bin = bin.display()),
        (Ok(bin), None) => eprintln!("Error in {bin}", bin = bin.display()),
        (Err(_), Some(msg)) => eprintln!("Error: {msg}"),
        (Err(_), None) => eprintln!("Error in program"),
    }

    eprint!("{}", error.display_stack());
    std::process::exit(code);
}

/// A trait for providing more graceful error reporting and aborting.
pub trait OrFail<T> {
    /// Unwraps the result, writing the error and any information in
    /// [`Error::source`] to stderr.
    ///
    /// A raw OS error code is used as the exit code when one is available in
    /// the source chain; otherwise, the process exits with code `1`.
    fn unwrap_or_fail(self) -> T;

    /// Unwraps the result, writing the provided message, the error, and any
    /// information in [`Error::source`] to stderr.
    ///
    /// A raw OS error code is used as the exit code when one is available in
    /// the source chain; otherwise, the process exits with code `1`.
    fn unwrap_or_die(self, msg: &str) -> T;

    /// Unwraps the result, writing the error and any information in
    /// [`Error::source`] to stderr and exiting with `code` on failure.
    fn unwrap_or_exit(self, code: i32) -> T;
}

impl<T, E> OrFail<T> for Result<T, E>
where
    E: Error + 'static,
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

    fn unwrap_or_exit(self, code: i32) -> T {
        match self {
            Ok(result) => result,
            Err(error) => {
                cold_path();
                exit_with_error(&error, None, code)
            }
        }
    }
}

/// A trait for providing more graceful error reporting and aborting. For
/// similar methods on [`Result`], see [`OrFail`].
///
/// A raw OS error code is used when one is available in the source chain;
/// otherwise, the process exits with code `1`. Any context available in
/// [`Error::source`] is displayed.
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
    E: Error + 'static,
{
    #[cold]
    fn fail(self) -> ! {
        let code = io_exit_code(&self);
        exit_with_error(&self, None, code)
    }

    #[cold]
    fn die(self, msg: &str) -> ! {
        let code = io_exit_code(&self);
        exit_with_error(&self, Some(msg), code)
    }
}

/// An error type supporting context and a backtrace.
///
/// Specifically, this error can hold up to three things:
///
/// 1. A line of context describing the error. This appears as one item in the
///    [`OrFail`] backtrace.
/// 2. Any subitems (additional indented lines with more information that appear
///    below the line of context). This is useful for including the values of
///    variables or other useful information.
/// 3. An optional source error, which this error wraps. Using
///    [`unwrap_or_fail`] or [`unwrap_or_die`] cause the source error to be
///    shown in the backtrace. This source is accessible via [`Error::source`].
///
/// This can be converted to [`std::io::Error`] with [`Into`]. Hence, in
/// functions returning [`std::io::Result`], the `?` operator can be used after
/// adding context.
///
/// A [`std::io::Error`] can also be converted into a [`ErrorWithContext`]
/// without adding additional context via [`From`]. The implementation ensures
/// that the error does not appear twice in the public source chain. This can be
/// useful to return `Result<(), ErrorWithContext>` from the `main` function,
/// which will automatically format any errors and the source chain upon
/// failure.
///
/// [`with_subitem`]: WithSubitem::with_subitem
/// [`unwrap_or_fail`]: OrFail::unwrap_or_fail
/// [`unwrap_or_die`]: OrFail::unwrap_or_die
#[must_use]
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

    /// Constructs a new [`ErrorWithContext`] with the given description and
    /// path, without a source error or any subitems.
    ///
    /// The `description` may be anything implementing `Display`. An allocation
    /// will occur for the resulting error message.
    pub fn new_with_path(description: impl Display, path: impl AsRef<Path>) -> Self {
        let description = format!("{description}: '{path}'", path = path.as_ref().display());

        ErrorWithContext {
            repr: Box::new(ErrorWithContextRepr {
                description,
                subitem: None,
                source: None,
            }),
        }
    }
}

impl Debug for ErrorWithContext {
    /// Formats the value using the debug formatter.
    ///
    /// By default, this will print the binary name (if available) followed by a
    /// formatted backtrace of the error. If using the alternate display with
    /// `{:#}`, a traditional debug format is used.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if f.alternate() {
            f.debug_struct("ErrorWithContext").field("repr", &self.repr).finish()
        } else {
            if let Ok(bin) = std::env::current_exe() {
                writeln!(f, "Error in {b}", b = bin.display())?;
            } else {
                writeln!(f, "Error in program")?;
            }

            writeln!(f, "{}", self.display_stack())
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
    source: Option<ErrorSource>,
}

/// The source error contained within an [`ErrorWithContext`].
#[derive(Debug)]
enum ErrorSource {
    /// An ordinary source which appears in the displayed error stack and is
    /// returned by [`Error::source`].
    Chained(Box<dyn Error + Send + Sync>),

    /// A source error which should be skipped in the source stack.
    /// [`Error::source`] skips this source, returning the source of the
    /// contained error instead.
    ///
    /// Errors cannot be directly removed from the stack because
    /// [`Error::source`] does not include the [`Send`] and [`Sync`] trait
    /// bounds, so this variant allows them to be skipped when traversing the
    /// stack.
    Transparent(Box<dyn Error + Send + Sync>),
}

impl ErrorSource {
    /// Returns the contained error regardless of whether it is transparent or
    /// not. This is used in [`io_exit_code`], which checks codes for all
    /// errors.
    fn as_error(&self) -> &(dyn Error + 'static) {
        match self {
            Self::Chained(source) | Self::Transparent(source) => source.as_ref(),
        }
    }
}

impl Display for ErrorWithContextRepr {
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

impl Display for ErrorWithContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.repr)
    }
}

impl Error for ErrorWithContext {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.repr.source {
            Some(ErrorSource::Chained(source)) => Some(source.as_ref()),
            Some(ErrorSource::Transparent(source)) => source.source(),
            None => None,
        }
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
}

impl<E: Error + Send + Sync + 'static> WithErrorContext for E {
    // Do not inline, since this is cold code
    fn with_context(self, description: impl Into<String>) -> ErrorWithContext {
        ErrorWithContext {
            repr: Box::new(ErrorWithContextRepr {
                description: description.into(),
                subitem:     None,
                source:      Some(ErrorSource::Chained(Box::new(self))),
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

        Self::with_context(self, description)
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

/// An extension trait for an [`Error`] enabling it to be displayed alongside
/// its error stack using [`Error::source`].
///
/// If aborting, consider using [`fail`] or [`die`]. Otherwise, this can be
/// helpful in displaying warnings.
///
/// [`fail`]: Fail::fail
/// [`die`]: Fail::die
pub trait DisplayErrStack {
    /// Returns a displayable representation of the error alongside its error
    /// stack using [`Error::source`].
    ///
    /// This using `→` and two spaces of indent before each item, and includes a
    /// newline at the end.
    fn display_stack(&self) -> ErrStackDisplay<'_>;
}

impl<E> DisplayErrStack for E
where
    E: Error + 'static,
{
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

impl DisplayErrStack for dyn Error + 'static {
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

impl DisplayErrStack for dyn Error + Send + 'static {
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

impl DisplayErrStack for dyn Error + Sync + 'static {
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

impl DisplayErrStack for dyn Error + Send + Sync + 'static {
    fn display_stack(&self) -> ErrStackDisplay<'_> {
        ErrStackDisplay(self)
    }
}

/// A display wrapper around an error that shows the error and its sources (with
/// [`Error::source`]) in a list, using `→` and two spaces of indent before each
/// item.
///
/// This includes a newline at the end.
pub struct ErrStackDisplay<'a>(&'a (dyn Error + 'static));

impl Display for ErrStackDisplay<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut maybe_err = Some(self.0);

        for _ in 0..MAX_ERROR_CHAIN_DEPTH {
            let Some(err) = maybe_err else {
                return Ok(());
            };

            writeln!(
                f,
                "  → {err}",
                err = IndentWrapper {
                    val:    err,
                    indent: "    ",
                }
            )?;

            maybe_err = err.source();
        }

        if maybe_err.is_some() {
            writeln!(f, "  → [error source chain truncated]")?;
        }

        Ok(())
    }
}

impl From<ErrorWithContext> for std::io::Error {
    #[inline]
    fn from(e: ErrorWithContext) -> Self {
        std::io::Error::other(e)
    }
}

impl From<std::io::Error> for ErrorWithContext {
    fn from(error: std::io::Error) -> Self {
        // This downcast avoids an extra allocation when doing a round-trip
        // conversion, and prevents the subitems from being squashed into the
        // description.
        let error = match error.downcast::<Self>() {
            Ok(error) => return error,
            Err(error) => error,
        };

        ErrorWithContext {
            repr: Box::new(ErrorWithContextRepr {
                description: error.to_string(),
                subitem:     None,
                source:      Some(ErrorSource::Transparent(Box::new(error))),
            }),
        }
    }
}

impl From<std::convert::Infallible> for ErrorWithContext {
    fn from(error: std::convert::Infallible) -> Self {
        match error {}
    }
}
